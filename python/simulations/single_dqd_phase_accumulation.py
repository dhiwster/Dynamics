"""
Single-DQD phase-accumulation draft simulation.

This script implements the 4 x 4 model from:

    DQD Individual Simulation with Phase Accumulation.md

Basis order:

    |up,R>, |up,L>, |down,R>, |down,L>

Hamiltonian, in micro-eV:

    H(t) = tc tau_x
         + epsilon(t) / 2 tau_z
         + Bz / 2 sigma_z
         + dBx / 2 sigma_x tau_z
         + Vac(t) cos(omega0 t + phase) tau_z

The default sequence is:

    1. EDSR Rx(pi/2) at epsilon = 0.
    2. Linear ramp to finite detuning.
    3. Hold at finite detuning.
    4. Linear ramp back to epsilon = 0.

The reported phase is the dynamic qubit phase accumulated relative to
staying at epsilon = 0:

    phi_detuning = integral [omega_ge(epsilon(t)) - omega_ge(0)] dt.

No QuTiP dependency is required; this draft uses numpy/scipy only.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp


HBAR_UEV_NS = 0.6582119569


@dataclass(frozen=True)
class DQDParams:
    """Physical parameters for one spin-charge DQD."""

    tc: float = 40.0
    dBx: float = 2.0
    Bz: float = 20.0
    Vac0: float = 1.0


@dataclass(frozen=True)
class PulseParams:
    """Pulse schedule in ns and micro-eV."""

    rx_angle: float = np.pi / 2
    epsilon_target: float = 2000.0
    ramp_time: float = 5.0
    hold_time: float = 20.0
    points_per_carrier_cycle: int = 24
    max_step_fraction_of_period: float = 1.0 / 20.0


@dataclass(frozen=True)
class Operators:
    """Single-DQD spin and charge operators."""

    sx: np.ndarray
    sy: np.ndarray
    sz: np.ndarray
    tx: np.ndarray
    ty: np.ndarray
    tz: np.ndarray


def pauli_operators() -> Operators:
    """Return Pauli operators in spin tensor charge ordering."""

    ident = np.eye(2, dtype=complex)
    sx_1q = np.array([[0, 1], [1, 0]], dtype=complex)
    sy_1q = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz_1q = np.array([[1, 0], [0, -1]], dtype=complex)

    return Operators(
        sx=np.kron(sx_1q, ident),
        sy=np.kron(sy_1q, ident),
        sz=np.kron(sz_1q, ident),
        tx=np.kron(ident, sx_1q),
        ty=np.kron(ident, sy_1q),
        tz=np.kron(ident, sz_1q),
    )


OPS = pauli_operators()


def static_hamiltonian(params: DQDParams, epsilon: float) -> np.ndarray:
    """Microscopic single-DQD Hamiltonian at fixed detuning, in micro-eV."""

    return (
        params.tc * OPS.tx
        + 0.5 * epsilon * OPS.tz
        + 0.5 * params.Bz * OPS.sz
        + 0.5 * params.dBx * (OPS.sx @ OPS.tz)
    )


def eigensystem(params: DQDParams, epsilon: float) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted eigenvalues and eigenvectors for H(epsilon)."""

    energies, vectors = np.linalg.eigh(static_hamiltonian(params, epsilon))
    order = np.argsort(energies)
    return energies[order], vectors[:, order]


def qubit_splitting(params: DQDParams, epsilon: float) -> float:
    """Ground-to-first-excited qubit splitting in angular frequency, rad/ns."""

    energies, _ = eigensystem(params, epsilon)
    return float((energies[1] - energies[0]) / HBAR_UEV_NS)


def edsr_matrix_element(params: DQDParams) -> float:
    """Magnitude of <g|tau_z|e> at epsilon = 0."""

    _, vectors = eigensystem(params, epsilon=0.0)
    g = vectors[:, 0]
    e = vectors[:, 1]
    return float(abs(np.vdot(g, OPS.tz @ e)))


def rx_duration(params: DQDParams, angle: float) -> float:
    """Resonant EDSR pulse duration for Rx(angle), in ns."""

    matrix_element = edsr_matrix_element(params)
    if matrix_element <= 0:
        raise ValueError("EDSR matrix element is zero; cannot calibrate Rx pulse")
    return abs(angle) * HBAR_UEV_NS / (params.Vac0 * matrix_element)


def make_epsilon_schedule(t_rx: float, pulse: PulseParams) -> Callable[[float], float]:
    """Piecewise epsilon(t): zero during Rx, ramp, hold, ramp back."""

    t_ramp_up_end = t_rx + pulse.ramp_time
    t_hold_end = t_ramp_up_end + pulse.hold_time
    t_ramp_down_end = t_hold_end + pulse.ramp_time

    def epsilon(t: float) -> float:
        if t < t_rx:
            return 0.0
        if t < t_ramp_up_end:
            return pulse.epsilon_target * (t - t_rx) / pulse.ramp_time
        if t < t_hold_end:
            return pulse.epsilon_target
        if t < t_ramp_down_end:
            return pulse.epsilon_target * (1.0 - (t - t_hold_end) / pulse.ramp_time)
        return 0.0

    return epsilon


def make_drive_envelope(t_rx: float, amplitude: float) -> Callable[[float], float]:
    """EDSR amplitude is on only during the initial resonant Rx pulse."""

    def envelope(t: float) -> float:
        return amplitude if 0.0 <= t < t_rx else 0.0

    return envelope


def lab_hamiltonian(
    params: DQDParams,
    epsilon_of_t: Callable[[float], float],
    drive_of_t: Callable[[float], float],
    carrier_omega: float,
    drive_phase: float,
) -> Callable[[float], np.ndarray]:
    """Build H(t) in micro-eV."""

    def hamiltonian(t: float) -> np.ndarray:
        eps = epsilon_of_t(t)
        drive = drive_of_t(t) * np.cos(carrier_omega * t + drive_phase)
        return static_hamiltonian(params, eps) + drive * OPS.tz

    return hamiltonian


def evolve_state(
    psi0: np.ndarray,
    hamiltonian: Callable[[float], np.ndarray],
    t_eval: np.ndarray,
    max_step: float,
) -> np.ndarray:
    """Integrate Schrodinger evolution and return states as columns."""

    def rhs(t: float, psi: np.ndarray) -> np.ndarray:
        return -1j * (hamiltonian(t) @ psi) / HBAR_UEV_NS

    result = solve_ivp(
        rhs,
        (float(t_eval[0]), float(t_eval[-1])),
        psi0,
        t_eval=t_eval,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
        max_step=max_step,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.y


def build_time_grid(
    params: DQDParams,
    pulse: PulseParams,
    t_rx: float,
    total_time: float,
) -> tuple[np.ndarray, float]:
    """Use the zero-detuning carrier period to set a stable lab-frame grid."""

    carrier_period = 2 * np.pi / qubit_splitting(params, epsilon=0.0)
    dt = carrier_period / pulse.points_per_carrier_cycle
    n_points = int(np.ceil(total_time / dt)) + 1
    t_eval = np.linspace(0.0, total_time, n_points)
    max_step = carrier_period * pulse.max_step_fraction_of_period
    return t_eval, max_step


def detuning_phase(
    params: DQDParams,
    epsilon_of_t: Callable[[float], float],
    t_eval: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute dynamic phase relative to epsilon = 0 along the schedule."""

    omega0 = qubit_splitting(params, epsilon=0.0)
    delta_omega = np.array(
        [qubit_splitting(params, epsilon_of_t(float(t))) - omega0 for t in t_eval]
    )
    phase = cumulative_trapezoid(delta_omega, t_eval, initial=0.0)
    return delta_omega, phase


def analyze_final_state(params: DQDParams, psi_final: np.ndarray) -> dict[str, complex | float]:
    """Project the final state onto the epsilon=0 qubit subspace."""

    _, vectors = eigensystem(params, epsilon=0.0)
    g = vectors[:, 0]
    e = vectors[:, 1]
    amp_g = np.vdot(g, psi_final)
    amp_e = np.vdot(e, psi_final)
    pop_g = float(abs(amp_g) ** 2)
    pop_e = float(abs(amp_e) ** 2)
    leakage = max(0.0, 1.0 - pop_g - pop_e)
    relative_phase = float(np.angle(amp_e) - np.angle(amp_g))
    relative_phase = float(np.angle(np.exp(1j * relative_phase)))
    return {
        "amp_g": amp_g,
        "amp_e": amp_e,
        "pop_g": pop_g,
        "pop_e": pop_e,
        "leakage": leakage,
        "relative_phase": relative_phase,
    }


def run_simulation(params: DQDParams, pulse: PulseParams, drive_phase: float = 0.0) -> dict:
    """Run the default phase-accumulation sequence."""

    t_rx = rx_duration(params, pulse.rx_angle)
    total_time = t_rx + 2 * pulse.ramp_time + pulse.hold_time
    t_eval, max_step = build_time_grid(params, pulse, t_rx, total_time)

    epsilon_of_t = make_epsilon_schedule(t_rx, pulse)
    drive_of_t = make_drive_envelope(t_rx, params.Vac0)
    carrier_omega = qubit_splitting(params, epsilon=0.0)
    hamiltonian = lab_hamiltonian(
        params=params,
        epsilon_of_t=epsilon_of_t,
        drive_of_t=drive_of_t,
        carrier_omega=carrier_omega,
        drive_phase=drive_phase,
    )

    _, vectors = eigensystem(params, epsilon=0.0)
    psi0 = vectors[:, 0]
    states = evolve_state(psi0, hamiltonian, t_eval, max_step=max_step)
    delta_omega, phase = detuning_phase(params, epsilon_of_t, t_eval)
    final = analyze_final_state(params, states[:, -1])

    return {
        "params": params,
        "pulse": pulse,
        "t_rx": t_rx,
        "total_time": total_time,
        "carrier_omega": carrier_omega,
        "carrier_period": 2 * np.pi / carrier_omega,
        "edsr_matrix_element": edsr_matrix_element(params),
        "t_eval": t_eval,
        "states": states,
        "epsilon": np.array([epsilon_of_t(float(t)) for t in t_eval]),
        "delta_omega": delta_omega,
        "phase": phase,
        "final": final,
    }


def print_summary(result: dict) -> None:
    """Print a compact numerical summary."""

    params: DQDParams = result["params"]
    pulse: PulseParams = result["pulse"]
    final = result["final"]
    phase = result["phase"][-1]

    print("Single-DQD phase-accumulation draft")
    print("-----------------------------------")
    print(f"tc={params.tc:g} ueV, dBx={params.dBx:g} ueV, Bz={params.Bz:g} ueV")
    print(f"Vac0={params.Vac0:g} ueV, epsilon_target={pulse.epsilon_target:g} ueV")
    print(f"<g|tau_z|e> at eps=0: {result['edsr_matrix_element']:.8f}")
    print(
        f"carrier omega_ge(0): {result['carrier_omega']:.6f} rad/ns "
        f"(period {result['carrier_period']:.6f} ns)"
    )
    print(f"Rx(pi/2) duration: {result['t_rx']:.6f} ns")
    print(f"ramp/hold/ramp: {pulse.ramp_time:g} / {pulse.hold_time:g} / {pulse.ramp_time:g} ns")
    print(f"total time: {result['total_time']:.6f} ns")
    print()
    print(f"detuning phase: {phase:.9f} rad ({np.degrees(phase):.6f} deg)")
    print(f"detuning phase mod 2pi: {np.angle(np.exp(1j * phase)):.9f} rad")
    print(f"final Pg={final['pop_g']:.9f}, Pe={final['pop_e']:.9f}, leakage={final['leakage']:.3e}")
    print(f"final relative qubit phase: {final['relative_phase']:.9f} rad")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tc", type=float, default=DQDParams.tc)
    parser.add_argument("--dBx", type=float, default=DQDParams.dBx)
    parser.add_argument("--Bz", type=float, default=DQDParams.Bz)
    parser.add_argument("--Vac0", type=float, default=DQDParams.Vac0)
    parser.add_argument("--epsilon-target", type=float, default=PulseParams.epsilon_target)
    parser.add_argument("--ramp-time", type=float, default=PulseParams.ramp_time)
    parser.add_argument("--hold-time", type=float, default=PulseParams.hold_time)
    parser.add_argument("--points-per-carrier-cycle", type=int, default=PulseParams.points_per_carrier_cycle)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = DQDParams(tc=args.tc, dBx=args.dBx, Bz=args.Bz, Vac0=args.Vac0)
    pulse = PulseParams(
        epsilon_target=args.epsilon_target,
        ramp_time=args.ramp_time,
        hold_time=args.hold_time,
        points_per_carrier_cycle=args.points_per_carrier_cycle,
    )
    result = run_simulation(params, pulse)
    print_summary(result)


if __name__ == "__main__":
    main()
