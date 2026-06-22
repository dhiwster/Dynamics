"""
solver.py — Simulation backend for the DQD gate sequence.

Imports the compiled Hamiltonian tools from circuit.py and provides two
distinct simulation modes:

  Mode A (run_state)    : time-domain state evolution via qt.mesolve
  Mode B (run_unitary)  : full unitary propagator via qt.propagator

Post-processing helpers extract DQD-only observables from the full
photon ⊗ DQD Hilbert space via partial trace and subspace projection.
"""

from __future__ import annotations

import numpy as np
import qutip as qt

from .system import DQDsystem
from .circuit import (
    DQDSequenceCompiler,
    iswap_H,
    xrot_H,
    zrot_H,
)


# ---------------------------------------------------------------------------
# State evolution
# ---------------------------------------------------------------------------

def run_state(
    H,
    psi0: qt.Qobj,
    tlist,
    c_ops=None,
    options=None,
) -> qt.solver.Result:
    """
    Evolve an initial state under H using qt.mesolve.

    Parameters
    ----------
    H       : Qobj or list — time-dependent Hamiltonian in QuTiP format,
              e.g. [H_drift, [H_ctrl, f(t)]] from DQDSequenceCompiler.compile()
    psi0    : initial state (ket or density matrix)
    tlist   : array-like of time points (ns)
    c_ops   : list of collapse operators; empty list → unitary evolution
    options : dict of solver options forwarded to qt.mesolve

    Returns
    -------
    result : qt.solver.Result  (.states holds the state at each tlist point)
    """
    c_ops = c_ops or []
    kwargs = {} if options is None else {"options": options}
    return qt.mesolve(H, psi0, tlist, c_ops=c_ops, **kwargs)


def dqd_populations(
    result: qt.solver.Result,
    sys: DQDsystem,
    state_index: int = -1,
    basis: str = 'spin_charge',
):
    """
    Post-process a mesolve result to extract DQD populations by partially
    tracing out the photon (cavity) subspace.

    Parameters
    ----------
    result      : qt.solver.Result from run_state
    sys         : DQDsystem instance
    state_index : which time-step to analyse (default: -1 → final state)
    basis       : 'energy'      — project onto eigenstates of sys.H_dqd,
                                  sorted by energy
                  'spin_charge' — diagonal of rho_dqd in the spin ⊗ charge
                                  computational basis

    Returns
    -------
    basis='energy'      → (p_ground, p_excited1) : float 2-tuple
    basis='spin_charge' → np.ndarray of shape (16,), basis order:
                          |↑R↑R⟩, |↑R↑L⟩, |↑R↓R⟩, |↑R↓L⟩,
                          |↑L↑R⟩, |↑L↑L⟩, |↑L↓R⟩, |↑L↓L⟩,
                          |↓R↑R⟩, |↓R↑L⟩, |↓R↓R⟩, |↓R↓L⟩,
                          |↓L↑R⟩, |↓L↑L⟩, |↓L↓R⟩, |↓L↓L⟩
    """

    state = result.states[state_index]
    rho_full = qt.ket2dm(state) if state.type == 'ket' else state

    # Partial trace: subsystem 0 = photon; keep subsystems 1–4 (DQD1 ⊗ DQD2)
    # Full dims: [[PHOTON_MAX, 2, 2, 2, 2], [PHOTON_MAX, 2, 2, 2, 2]]
    rho_dqd = rho_full.ptrace([1, 2, 3, 4])

    if basis == 'energy':
        _, eigvecs = sys.H_dqd.eigenstates()
        gs = eigvecs[0]
        e1 = eigvecs[1]
        p_ground   = float(qt.expect(gs * gs.dag(), rho_dqd).real)
        p_excited1 = float(qt.expect(e1 * e1.dag(), rho_dqd).real)
        return p_ground, p_excited1

    elif basis == 'spin_charge':
        return rho_dqd.full().diagonal().real.copy()

    else:
        raise ValueError(f"basis must be 'energy' or 'spin_charge', got {basis!r}")


# ---------------------------------------------------------------------------
# Unitary evolution
# ---------------------------------------------------------------------------

def run_unitary(H, tlist, options=None) -> qt.Qobj:
    """
    Compute the time-evolution operator U(T) over the full pulse duration.

    qt.propagator integrates dU/dt = −i H U from t=0 to t=T=tlist[-1].
    Passing the full tlist gives the ODE solver fine-grained time steps,
    which is important for accuracy with rapidly oscillating time-dependent
    Hamiltonians (e.g. the EDSR drive).

    Parameters
    ----------
    H       : Qobj or list — Hamiltonian in QuTiP format,
              e.g. [H_drift, [H_ctrl, f(t)]] from DQDSequenceCompiler.compile()
    tlist   : array-like of time points (ns); propagator at tlist[-1] is returned
    options : dict of solver options forwarded to qt.propagator

    Returns
    -------
    U : qt.Qobj — full unitary at t = tlist[-1]
        dims = [[PHOTON_MAX, 2, 2, 2, 2], [PHOTON_MAX, 2, 2, 2, 2]]
    """
    kwargs = {} if options is None else {"options": options}
    result = qt.propagator(H, tlist, **kwargs)
    return result[-1] if isinstance(result, list) else result


def vacuum_subspace_unitary(U: qt.Qobj, sys=None) -> np.ndarray:
    """
    Extract the zero-photon (vacuum) block of the propagator U.

    In the QuTiP tensor ordering photon ⊗ DQD, the state |n=0, k⟩ maps to
    matrix index k (k = 0 … 15), so the vacuum block is the top-left 16×16
    sub-matrix of U expressed in the computational basis.

    Parameters
    ----------
    U   : qt.Qobj — full propagator from run_unitary
    sys : DQDsystem instance (accepted for API consistency; not used here)

    Returns
    -------
    U_vac : np.ndarray, shape (16, 16), complex
        Matrix elements ⟨0, i | U | 0, j⟩ for all 16 DQD basis states i, j.
    """
    dqd_dim = 16   # 4 (DQD1: spin ⊗ charge) × 4 (DQD2: spin ⊗ charge)
    return U.full()[:dqd_dim, :dqd_dim]


def dqd_qubit_unitary(U: qt.Qobj, sys: DQDsystem) -> np.ndarray:
    """
    Project the propagator onto the 2D qubit subspace in the vacuum sector.

    The qubit subspace is spanned by:
        |0 photons, DQD ground state⟩  and  |0 photons, DQD first excited state⟩

    DQD eigenstates are obtained by diagonalising sys.H_dqd.

    Parameters
    ----------
    U   : qt.Qobj — full propagator from run_unitary
    sys : DQDsystem instance

    Returns
    -------
    U_q : np.ndarray, shape (2, 2), complex
        U_q[i, j] = ⟨0, vᵢ | U | 0, vⱼ⟩
        where v₀ = DQD ground state, v₁ = DQD first excited state.
    """
    U_vac = vacuum_subspace_unitary(U)

    _, eigvecs = sys.H_dqd.eigenstates()
    v0 = eigvecs[0].full().flatten()
    v1 = eigvecs[1].full().flatten()

    return np.array([
        [v0.conj() @ U_vac @ v0,  v0.conj() @ U_vac @ v1],
        [v1.conj() @ U_vac @ v0,  v1.conj() @ U_vac @ v1],
    ])


def two_qubit_unitary(U_full: qt.Qobj, sys: DQDsystem) -> np.ndarray:
    """
    Project the propagator onto the two-qubit subspace in the vacuum sector.

    The two-qubit basis is built from tensor products of individual DQD
    ground/excited eigenstates at ε=0:
        |00⟩ = |g,g⟩,  |01⟩ = |g,e⟩,  |10⟩ = |e,g⟩,  |11⟩ = |e,e⟩

    Parameters
    ----------
    U_full : qt.Qobj — full propagator (photon ⊗ DQD1 ⊗ DQD2) from run_unitary
    sys    : DQDsystem instance

    Returns
    -------
    U_4x4 : np.ndarray, shape (4, 4), complex
        ⟨qi | U_vac | qj⟩ in the zero-photon subspace.
    """
    U_vac = vacuum_subspace_unitary(U_full)   # 16 × 16

    H_single = sys.tc * sys.tx + sys.Bz / 2 * sys.sz + sys.bx / 2 * sys.sx * sys.tz
    _, evecs = H_single.eigenstates(sort='low')
    g = evecs[0].full().flatten()   # 4-dim DQD ground state
    e = evecs[1].full().flatten()   # 4-dim DQD first excited state

    q_basis = np.array([
        np.kron(g, g),   # |00⟩
        np.kron(g, e),   # |01⟩
        np.kron(e, g),   # |10⟩
        np.kron(e, e),   # |11⟩
    ])

    return q_basis.conj() @ U_vac @ q_basis.T
