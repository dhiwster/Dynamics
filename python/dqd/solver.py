"""DQD-specific post-processing helpers for simulated evolution."""

from __future__ import annotations

import numpy as np
import qutip as qt

from .system import DQDsystem


def dqd_populations(
    result: qt.solver.Result,
    sys: DQDsystem,
    state_index: int = -1,
    basis: str = "spin_charge",
):
    """
    Extract DQD populations after tracing out the photon subspace.

    ``basis="energy"`` projects onto the lowest two eigenstates of
    ``sys.H_dqd``. ``basis="spin_charge"`` returns the diagonal of the DQD
    reduced density matrix in the computational basis.
    """

    state = result.states[state_index]
    rho_full = qt.ket2dm(state) if state.type == "ket" else state
    rho_dqd = rho_full.ptrace([1, 2, 3, 4])

    if basis == "energy":
        _, eigvecs = sys.H_dqd.eigenstates()
        gs = eigvecs[0]
        e1 = eigvecs[1]
        p_ground = float(qt.expect(gs * gs.dag(), rho_dqd).real)
        p_excited1 = float(qt.expect(e1 * e1.dag(), rho_dqd).real)
        return p_ground, p_excited1

    if basis == "spin_charge":
        return rho_dqd.full().diagonal().real.copy()

    raise ValueError(f"basis must be 'energy' or 'spin_charge', got {basis!r}")


def vacuum_subspace_unitary(U: qt.Qobj, sys=None) -> np.ndarray:
    """Extract the zero-photon block of the propagator."""

    del sys
    dqd_dim = 16
    return U.full()[:dqd_dim, :dqd_dim]


def dqd_qubit_unitary(U: qt.Qobj, sys: DQDsystem) -> np.ndarray:
    """Project the propagator onto the DQD qubit subspace in the vacuum sector."""

    U_vac = vacuum_subspace_unitary(U)

    _, eigvecs = sys.H_dqd.eigenstates()
    v0 = eigvecs[0].full().flatten()
    v1 = eigvecs[1].full().flatten()

    return np.array(
        [
            [v0.conj() @ U_vac @ v0, v0.conj() @ U_vac @ v1],
            [v1.conj() @ U_vac @ v0, v1.conj() @ U_vac @ v1],
        ]
    )


def two_qubit_unitary(U_full: qt.Qobj, sys: DQDsystem) -> np.ndarray:
    """Project the propagator onto the two-qubit vacuum-sector subspace."""

    U_vac = vacuum_subspace_unitary(U_full)

    H_single = sys.tc * sys.tx + sys.Bz / 2 * sys.sz + sys.bx / 2 * sys.sx * sys.tz
    _, evecs = H_single.eigenstates(sort="low")
    g = evecs[0].full().flatten()
    e = evecs[1].full().flatten()

    q_basis = np.array(
        [
            np.kron(g, g),
            np.kron(g, e),
            np.kron(e, g),
            np.kron(e, e),
        ]
    )

    return q_basis.conj() @ U_vac @ q_basis.T


__all__ = [
    "dqd_populations",
    "dqd_qubit_unitary",
    "two_qubit_unitary",
    "vacuum_subspace_unitary",
]
