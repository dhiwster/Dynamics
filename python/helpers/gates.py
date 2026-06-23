from __future__ import annotations

import numpy as np
import qutip as qt

from .constants import RAD_PER_NS_PER_MICROELECTRONVOLT
from .pulses import Pulse, constant
from .system import DQDsystem

_MHz2GHz = 1e-3


# ---------------------------------------------------------------------------
# Internal operator helpers
# ---------------------------------------------------------------------------

def _H0_ns(dqd: DQDsystem) -> qt.Qobj:
    return dqd.H_static * _MHz2GHz


def _eps1_op(dqd: DQDsystem) -> qt.Qobj:
    return qt.tensor(qt.qeye(dqd.photon_max), dqd.tz1) * (
        0.5 * RAD_PER_NS_PER_MICROELECTRONVOLT
    )


def _eps2_op(dqd: DQDsystem) -> qt.Qobj:
    return qt.tensor(qt.qeye(dqd.photon_max), dqd.tz2) * (
        0.5 * RAD_PER_NS_PER_MICROELECTRONVOLT
    )


def _edsr1_op(dqd: DQDsystem) -> qt.Qobj:
    return dqd.H_edsr1_amplitude * RAD_PER_NS_PER_MICROELECTRONVOLT


def _edsr2_op(dqd: DQDsystem) -> qt.Qobj:
    return dqd.H_edsr2_amplitude * RAD_PER_NS_PER_MICROELECTRONVOLT


# ---------------------------------------------------------------------------
# Calibration helpers
# ---------------------------------------------------------------------------

def zrot_delta_freq_GHz(dqd: DQDsystem, eps_amp: float) -> float:
    """Qubit eigenfrequency shift ΔEσ/ħ [GHz] caused by a DC ε detuning."""
    tc, bx, Bz = dqd.tc, dqd.bx, dqd.Bz
    tx, sx, sz, tz = dqd.tx, dqd.sx, dqd.sz, dqd.tz

    H_base = tc * tx + Bz / 2 * sz + bx / 2 * sx * tz
    H_det  = H_base + eps_amp / 2 * tz

    e_base = np.sort(H_base.eigenenergies())
    e_det  = np.sort(H_det.eigenenergies())

    return abs((e_det[1] - e_det[0]) - (e_base[1] - e_base[0])) * RAD_PER_NS_PER_MICROELECTRONVOLT


def zrot_time_ns(dqd: DQDsystem, eps_amp: float, angle: float) -> float:
    """Pulse duration [ns] for Rz(angle) via a DC ε detuning."""
    delta_f = zrot_delta_freq_GHz(dqd, eps_amp)
    if delta_f == 0:
        raise ValueError("eps_amp=0 gives zero qubit frequency shift")
    return abs(angle) / (2 * np.pi * delta_f)


# ---------------------------------------------------------------------------
# Gate Hamiltonian builders
# ---------------------------------------------------------------------------

def iswap_H(dqd: DQDsystem) -> list:
    """iSWAP Hamiltonian at ε=0. Gate duration: dqd.iSWAP_gate_time() * 1e3 ns."""
    return [_H0_ns(dqd)]


def xrot_H(
    dqd: DQDsystem,
    target: int,
    t_start: float,
    t_end: float,
    eps_idle: float | None = None,
) -> list:
    """EDSR X-rotation Hamiltonian for target DQD (1 or 2)."""
    if target not in (1, 2):
        raise ValueError("target must be 1 or 2")
    eps_idle   = dqd.epsilon_idle if eps_idle is None else eps_idle
    Esigma_GHz = dqd.Esigma * RAD_PER_NS_PER_MICROELECTRONVOLT
    H0         = _H0_ns(dqd)
    H_drv, H_idl = (_edsr1_op(dqd), _eps2_op(dqd)) if target == 1 else (_edsr2_op(dqd), _eps1_op(dqd))

    def _edsr(t, **kwargs):
        return np.cos(Esigma_GHz * t) if t_start <= t < t_end else 0.0

    def _idle(t, **kwargs):
        return eps_idle if not (t_start <= t < t_end) else 0.0

    return [H0, [H_drv, _edsr], [H_idl, _idle]]


def zrot_H(
    dqd: DQDsystem,
    target: int,
    t_start: float,
    t_end: float,
    eps_amp: float,
    eps_idle: float | None = None,
) -> list:
    """Z-rotation Hamiltonian via ε detuning pulse on target DQD (1 or 2)."""
    if target not in (1, 2):
        raise ValueError("target must be 1 or 2")
    eps_idle = dqd.epsilon_idle if eps_idle is None else eps_idle
    H0       = _H0_ns(dqd)
    H_tgt, H_idl = (_eps1_op(dqd), _eps2_op(dqd)) if target == 1 else (_eps2_op(dqd), _eps1_op(dqd))

    def _eps_pulse(t, **kwargs):
        return eps_amp if t_start <= t < t_end else 0.0

    def _idle_full(t, **kwargs):
        return eps_idle if not (t_start <= t < t_end) else 0.0

    return [H0, [H_tgt, _eps_pulse], [H_idl, _idle_full]]


# ---------------------------------------------------------------------------
# Eigenstate preparation
# ---------------------------------------------------------------------------

def vacuum_eigenstates(dqd: DQDsystem) -> list:
    """Ground and first excited eigenstates of the single-DQD H₀ (no photon)."""
    H0_single = dqd.tc * dqd.tx + dqd.Bz / 2 * dqd.sz + dqd.bx / 2 * dqd.sx * dqd.tz
    _, states = H0_single.eigenstates(sort='low')
    return list(states)


def initial_full_state(
    dqd: DQDsystem,
    dqd1_eigenstate_idx: int = 0,
    dqd2_eigenstate_idx: int = 1,
) -> qt.Qobj:
    """Build |ψ₀⟩ = |n=0⟩_photon ⊗ |φ_i⟩_DQD1 ⊗ |φ_j⟩_DQD2."""
    states = vacuum_eigenstates(dqd)
    return qt.tensor(
        qt.basis(dqd.photon_max, 0),
        states[dqd1_eigenstate_idx],
        states[dqd2_eigenstate_idx],
    )


# ---------------------------------------------------------------------------
# High-level Hamiltonian assembly (from builders.py)
# ---------------------------------------------------------------------------

def _delta(base: float, pulse: Pulse) -> Pulse:
    def f(t: float, args: dict | None = None, **kwargs) -> float:
        if kwargs:
            args = {**(args or {}), **kwargs}
        return pulse(t, args) - base
    return f


def build_H(
    dqd: DQDsystem,
    tc1: Pulse | None = None,
    tc2: Pulse | None = None,
    eps1: Pulse | None = None,
    eps2: Pulse | None = None,
    same_for_both: bool = True,
) -> list:
    """Build a DQD QuTiP Hamiltonian list from absolute control pulses."""
    tc1 = tc1 or constant(dqd.tc)
    eps1 = eps1 or constant(dqd.epsilon)

    if same_for_both:
        tc2 = tc2 or tc1
        eps2 = eps2 or eps1
    else:
        tc2 = tc2 or constant(dqd.tc)
        eps2 = eps2 or constant(dqd.epsilon)

    return [
        dqd.H_static,
        [dqd.H_tc1_op, _delta(dqd.tc, tc1)],
        [dqd.H_tc2_op, _delta(dqd.tc, tc2)],
        [dqd.H_eps1_op, eps1],
        [dqd.H_eps2_op, eps2],
    ]
