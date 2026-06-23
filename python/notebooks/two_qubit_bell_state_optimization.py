from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import qutip as qt
from qutip import basis, qeye, tensor

repo_root = Path(__file__).resolve()
for candidate in [repo_root, *repo_root.parents]:
    if (candidate / "python").exists() and (candidate / "pyproject.toml").exists():
        repo_root = candidate
        break

python_root = repo_root / "python"
if str(python_root) not in sys.path:
    sys.path.insert(0, str(python_root))

from helpers import DQDsystem  # noqa: E402
from helpers import RAD_PER_US_PER_MICROELECTRONVOLT, TWO_PI  # noqa: E402


@dataclass(frozen=True)
class WorkflowConfig:
    tc_ghz: float = 10.0
    bx_uev: float = 3.0
    bz_uev: float = 35.0
    wc_mhz: float = 5.5e3
    gc_mhz: float = 50.0
    photon_max: int = 10
    epsilon_prepare_ghz: float = 80.0
    omega_drive_mhz: float = 100.0
    epsilon_ramp_cycles: float = 20.0
    ramp_points_per_charge_period: int = 80


@dataclass(frozen=True)
class OptimizationConfig:
    prep_scan_ns: float = 10.0
    prep_points: int = 10001
    prep_window_ns: float = 0.4
    prep_stride_ns: float = 0.025
    gate_max_fraction: float = 0.35
    gate_points: int = 241
    refine_points: int = 41
    coarse_top_k: int = 3
    score_metric: str = "readout_rr_weighted_phase_optimized_bell_overlap"


class ConstantHamiltonianEvolution:
    def __init__(self, hamiltonian: qt.Qobj):
        evals, evecs = hamiltonian.eigenstates()
        self.evals = np.array(evals, dtype=float)
        self.basis_matrix = np.column_stack([vec.full().ravel() for vec in evecs])

    def evolve_one(self, psi0: qt.Qobj, t_us: float) -> qt.Qobj:
        psi0_vec = psi0.full().ravel()
        coeffs = self.basis_matrix.conj().T @ psi0_vec
        state_vec = self.basis_matrix @ (coeffs * np.exp(-1j * self.evals * t_us))
        return qt.Qobj(state_vec[:, None], dims=psi0.dims)

    def evolve_many(self, psi0: qt.Qobj, tlist_us: np.ndarray) -> list[qt.Qobj]:
        psi0_vec = psi0.full().ravel()
        coeffs = self.basis_matrix.conj().T @ psi0_vec
        return [
            qt.Qobj(
                self.basis_matrix @ (coeffs * np.exp(-1j * self.evals * float(t_us)))[:, None],
                dims=psi0.dims,
            )
            for t_us in tlist_us
        ]


def ghz_to_uev(ghz: float) -> float:
    return TWO_PI * ghz * 1_000.0 / RAD_PER_US_PER_MICROELECTRONVOLT


def lab_frame_product_state(photon_max: int, spin1: str = "up", spin2: str = "down") -> qt.Qobj:
    spin_index = {"up": 0, "down": 1}
    psi1 = basis([2, 2], [spin_index[spin1], 0])
    psi2 = basis([2, 2], [spin_index[spin2], 0])
    return tensor(basis(photon_max, 0), psi1, psi2)


def spin_diagonalization_data(dqd: DQDsystem) -> dict[str, float]:
    theta = np.arctan2(dqd.bx, dqd.Bz)
    return {
        "theta": float(theta),
        "spin_field": float(np.hypot(dqd.Bz, dqd.bx)),
        "cos": float(np.cos(theta)),
        "sin": float(np.sin(theta)),
    }


def single_spin_diagonalization_unitary(dqd: DQDsystem) -> qt.Qobj:
    theta = spin_diagonalization_data(dqd)["theta"]
    return (0.5j * theta * dqd.sy * dqd.tz).expm()


def full_spin_diagonalization_unitary(dqd: DQDsystem) -> qt.Qobj:
    u_spin = single_spin_diagonalization_unitary(dqd)
    return tensor(qeye(dqd.photon_max), u_spin, u_spin)


def spin_diagonalized_single_hamiltonian(dqd: DQDsystem, epsilon_uev: float) -> qt.Qobj:
    data = spin_diagonalization_data(dqd)
    return (
        dqd.tc * (data["cos"] * dqd.tx - data["sin"] * dqd.sy * dqd.ty)
        + 0.5 * epsilon_uev * dqd.tz
        + 0.5 * data["spin_field"] * dqd.sz
    )


def spin_basis_qubit_resonance_rad_per_us(dqd: DQDsystem, epsilon_uev: float) -> float:
    energies = np.sort(spin_diagonalized_single_hamiltonian(dqd, epsilon_uev).eigenenergies())
    return float((energies[1] - energies[0]) * RAD_PER_US_PER_MICROELECTRONVOLT)


def transform_full_operator_to_spin_basis(op: qt.Qobj, dqd: DQDsystem) -> qt.Qobj:
    u_full = full_spin_diagonalization_unitary(dqd)
    return (u_full * op * u_full.dag()).tidyup(atol=1e-10)


def spin_populations(states: list[qt.Qobj], spin_index: int) -> tuple[np.ndarray, np.ndarray]:
    proj_up = basis(2, 0) * basis(2, 0).dag()
    proj_down = basis(2, 1) * basis(2, 1).dag()
    rho_spin = [qt.ket2dm(state).ptrace(spin_index) for state in states]
    p_up = np.array([qt.expect(proj_up, rho) for rho in rho_spin], dtype=float)
    p_down = np.array([qt.expect(proj_down, rho) for rho in rho_spin], dtype=float)
    return p_up, p_down


def mixed_state_concurrence(rho: qt.Qobj) -> float:
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    yy = np.kron(sigma_y, sigma_y)
    rho_mat = rho.full()
    r_mat = rho_mat @ yy @ rho_mat.conj() @ yy
    eigvals = np.linalg.eigvals(r_mat)
    lambdas = np.sort(np.sqrt(np.maximum(np.real_if_close(eigvals), 0.0)))[::-1]
    return float(max(0.0, lambdas[0] - lambdas[1] - lambdas[2] - lambdas[3]))


def spin_basis_two_spin_density_matrix(state: qt.Qobj, dqd: DQDsystem) -> qt.Qobj:
    return qt.ket2dm(full_spin_diagonalization_unitary(dqd) * state).ptrace([1, 3])


def spin_basis_rr_spin_density_matrix(state: qt.Qobj, dqd: DQDsystem) -> tuple[qt.Qobj, float]:
    psi_spin_basis = full_spin_diagonalization_unitary(dqd) * state
    rho_spin_basis = qt.ket2dm(psi_spin_basis)
    proj_r = basis(2, 0) * basis(2, 0).dag()
    proj_rr = tensor(qeye(dqd.photon_max), qeye(2), proj_r, qeye(2), proj_r)
    rho_rr_full = proj_rr * rho_spin_basis * proj_rr
    rr_weight = float(rho_rr_full.tr().real)
    rho_rr_spin = rho_rr_full.ptrace([1, 3])
    if rr_weight > 1e-15:
        rho_rr_spin = rho_rr_spin / rr_weight
    return rho_rr_spin, rr_weight


def bell_readout_metrics(state: qt.Qobj, dqd: DQDsystem, bell_target: np.ndarray) -> dict[str, float]:
    bell_ket = qt.Qobj(bell_target[:, None], dims=[[2, 2], [1, 1]])
    spin_rho = spin_basis_two_spin_density_matrix(state, dqd)
    rr_spin_rho, rr_weight = spin_basis_rr_spin_density_matrix(state, dqd)
    spin_bell_overlap = float(qt.expect(bell_ket * bell_ket.dag(), spin_rho).real)
    rr_bell_overlap = float(qt.expect(bell_ket * bell_ket.dag(), rr_spin_rho).real)
    spin_phase_optimized = phase_optimized_odd_bell_overlap(spin_rho)
    rr_phase_optimized = phase_optimized_odd_bell_overlap(rr_spin_rho)
    spin_phase = phase_optimized_odd_bell_phase(spin_rho)
    rr_phase = phase_optimized_odd_bell_phase(rr_spin_rho)
    return {
        "spin_bell_overlap": spin_bell_overlap,
        "spin_phase_optimized_bell_overlap": spin_phase_optimized,
        "spin_phase_optimized_bell_phase_rad": spin_phase,
        "spin_concurrence": mixed_state_concurrence(spin_rho),
        "rr_spin_weight": rr_weight,
        "rr_spin_bell_overlap": rr_bell_overlap,
        "rr_spin_phase_optimized_bell_overlap": rr_phase_optimized,
        "rr_spin_phase_optimized_bell_phase_rad": rr_phase,
        "rr_spin_weighted_bell_overlap": rr_weight * rr_bell_overlap,
        "rr_spin_weighted_phase_optimized_bell_overlap": rr_weight * rr_phase_optimized,
        "rr_spin_concurrence": mixed_state_concurrence(rr_spin_rho) if rr_weight > 1e-15 else 0.0,
    }


def phase_optimized_odd_bell_overlap(rho: qt.Qobj) -> float:
    """Max overlap with (|+-> + exp(i phi)|-+>)/sqrt(2)."""
    mat = rho.full()
    return float((0.5 * (mat[1, 1].real + mat[2, 2].real) + abs(mat[1, 2])).real)


def phase_optimized_odd_bell_phase(rho: qt.Qobj) -> float:
    """Phase phi for (|+-> + exp(i phi)|-+>)/sqrt(2)."""
    return float(-np.angle(rho.full()[1, 2]))


def ramp_detuning(
    dqd: DQDsystem,
    psi0: qt.Qobj,
    epsilon_start_uev: float,
    epsilon_stop_uev: float,
    ramp_time_us: float,
    ramp_samples: int,
    options: dict[str, Any],
) -> qt.Qobj:
    def ramp_coeff(t, epsilon_start_uev, epsilon_stop_uev, ramp_time_us):
        if ramp_time_us <= 0.0:
            return epsilon_stop_uev
        s = min(max(t / ramp_time_us, 0.0), 1.0)
        return epsilon_start_uev + (epsilon_stop_uev - epsilon_start_uev) * s

    tlist_us = np.linspace(0.0, ramp_time_us, ramp_samples)
    result = qt.mesolve(
        [dqd.H_static, [dqd.H_eps1_op, ramp_coeff], [dqd.H_eps2_op, ramp_coeff]],
        psi0,
        tlist_us,
        c_ops=[],
        e_ops=[],
        args={
            "epsilon_start_uev": epsilon_start_uev,
            "epsilon_stop_uev": epsilon_stop_uev,
            "ramp_time_us": ramp_time_us,
        },
        options=options,
    )
    return result.states[-1]


def build_workflow(config: WorkflowConfig) -> dict[str, Any]:
    tc = ghz_to_uev(config.tc_ghz)
    dqd = DQDsystem(
        tc=tc,
        bx=config.bx_uev,
        Bz=config.bz_uev,
        Vac0=0.0,
        wc=TWO_PI * config.wc_mhz,
        gc=TWO_PI * config.gc_mhz,
        photon_max=config.photon_max,
        epsilon_idle=0.0,
    )
    epsilon_prepare_uev = ghz_to_uev(config.epsilon_prepare_ghz)
    charge_period_ns = TWO_PI / (2.0 * tc * (RAD_PER_US_PER_MICROELECTRONVOLT / 1000.0))
    epsilon_ramp_ns = config.epsilon_ramp_cycles * charge_period_ns
    omega_drive = spin_basis_qubit_resonance_rad_per_us(dqd, epsilon_prepare_uev)
    return {
        "dqd": dqd,
        "epsilon_prepare_uev": epsilon_prepare_uev,
        "omega_drive": omega_drive,
        "omega_drive_ghz": omega_drive / (TWO_PI * 1000.0),
        "omega_drive_mhz": TWO_PI * config.omega_drive_mhz,
        "epsilon_ramp_us": epsilon_ramp_ns / 1000.0,
        "ramp_samples": int(np.ceil(config.epsilon_ramp_cycles * config.ramp_points_per_charge_period)) + 1,
        "gate_time_us": float(dqd.iSWAP_gate_time()),
    }


def run_prep_scan(
    workflow: dict[str, Any],
    opt: OptimizationConfig,
    solver_options: dict[str, Any],
) -> dict[str, Any]:
    dqd = workflow["dqd"]
    u_full_spin = full_spin_diagonalization_unitary(dqd)
    psi0_lab = lab_frame_product_state(dqd.photon_max, spin1="down", spin2="down")
    psi0_spin_basis = u_full_spin * psi0_lab
    h_drive_lab = tensor(qeye(dqd.photon_max), dqd.sx1)
    h_static_lab = (
        dqd.H_static
        + workflow["epsilon_prepare_uev"] * dqd.H_eps1_op
        + workflow["epsilon_prepare_uev"] * dqd.H_eps2_op
    )
    h_drive = transform_full_operator_to_spin_basis(h_drive_lab, dqd)
    h_static = transform_full_operator_to_spin_basis(h_static_lab, dqd)

    def drive_coeff(t, omega_drive_mhz, omega_drive):
        return omega_drive_mhz * np.cos(omega_drive * t)

    prep_scan_ns = np.linspace(0.0, opt.prep_scan_ns, opt.prep_points)
    prep_scan_us = prep_scan_ns / 1000.0
    result = qt.mesolve(
        [h_static, [h_drive, drive_coeff]],
        psi0_spin_basis,
        prep_scan_us,
        c_ops=[],
        e_ops=[],
        args={
            "omega_drive_mhz": workflow["omega_drive_mhz"],
            "omega_drive": workflow["omega_drive"],
        },
        options=solver_options,
    )
    p_up, _ = spin_populations(result.states, spin_index=1)
    return {
        "states_spin_basis": result.states,
        "times_ns": prep_scan_ns,
        "p_up": p_up,
        "u_full_spin": u_full_spin,
    }


def prep_candidate_indices(prep: dict[str, Any], opt: OptimizationConfig) -> np.ndarray:
    center_idx = int(np.argmax(prep["p_up"]))
    times_ns = prep["times_ns"]
    center_ns = float(times_ns[center_idx])
    low = center_ns - opt.prep_window_ns
    high = center_ns + opt.prep_window_ns
    stride = max(1, int(round(opt.prep_stride_ns / float(times_ns[1] - times_ns[0]))))
    indices = np.where((times_ns >= low) & (times_ns <= high))[0]
    return indices[::stride]


def score_candidate(candidate: dict[str, Any], metric: str) -> float:
    if metric == "readout_rr_weighted_bell_overlap":
        return candidate["readout_rr_spin_weighted_bell_overlap"]
    if metric == "readout_rr_weighted_phase_optimized_bell_overlap":
        return candidate["readout_rr_spin_weighted_phase_optimized_bell_overlap"]
    if metric == "readout_rr_bell_overlap":
        return candidate["readout_rr_spin_bell_overlap"]
    if metric == "readout_rr_phase_optimized_bell_overlap":
        return candidate["readout_rr_spin_phase_optimized_bell_overlap"]
    if metric == "readout_spin_bell_overlap":
        return candidate["readout_spin_bell_overlap"]
    if metric == "readout_spin_phase_optimized_bell_overlap":
        return candidate["readout_spin_phase_optimized_bell_overlap"]
    if metric == "gate_rr_weighted_bell_overlap":
        return candidate["gate_rr_spin_weighted_bell_overlap"]
    if metric == "gate_rr_weighted_phase_optimized_bell_overlap":
        return candidate["gate_rr_spin_weighted_phase_optimized_bell_overlap"]
    if metric == "gate_rr_bell_overlap":
        return candidate["gate_rr_spin_bell_overlap"]
    if metric == "gate_rr_phase_optimized_bell_overlap":
        return candidate["gate_rr_spin_phase_optimized_bell_overlap"]
    if metric == "gate_spin_bell_overlap":
        return candidate["gate_spin_bell_overlap"]
    if metric == "gate_spin_phase_optimized_bell_overlap":
        return candidate["gate_spin_phase_optimized_bell_overlap"]
    raise ValueError(f"unknown score metric: {metric}")


def coarse_gate_score(item: dict[str, Any], metric: str) -> float:
    metric_to_gate_key = {
        "readout_rr_weighted_bell_overlap": "rr_spin_weighted_bell_overlap",
        "readout_rr_weighted_phase_optimized_bell_overlap": (
            "rr_spin_weighted_phase_optimized_bell_overlap"
        ),
        "readout_rr_bell_overlap": "rr_spin_bell_overlap",
        "readout_rr_phase_optimized_bell_overlap": "rr_spin_phase_optimized_bell_overlap",
        "readout_spin_bell_overlap": "spin_bell_overlap",
        "readout_spin_phase_optimized_bell_overlap": "spin_phase_optimized_bell_overlap",
        "gate_rr_weighted_bell_overlap": "rr_spin_weighted_bell_overlap",
        "gate_rr_weighted_phase_optimized_bell_overlap": (
            "rr_spin_weighted_phase_optimized_bell_overlap"
        ),
        "gate_rr_bell_overlap": "rr_spin_bell_overlap",
        "gate_rr_phase_optimized_bell_overlap": "rr_spin_phase_optimized_bell_overlap",
        "gate_spin_bell_overlap": "spin_bell_overlap",
        "gate_spin_phase_optimized_bell_overlap": "spin_phase_optimized_bell_overlap",
    }
    return item["gate_metrics"][metric_to_gate_key[metric]]


def optimize_bell_state(
    workflow: dict[str, Any],
    opt: OptimizationConfig,
    solver_options: dict[str, Any],
) -> dict[str, Any]:
    dqd = workflow["dqd"]
    bell_target = np.array([0.0, 1.0, 1.0j, 0.0], dtype=complex) / np.sqrt(2.0)
    prep = run_prep_scan(workflow, opt, solver_options)
    gate_evolver = ConstantHamiltonianEvolution(dqd.H_static)
    prep_indices = prep_candidate_indices(prep, opt)
    coarse_gate_tlist_us = np.linspace(
        0.0,
        opt.gate_max_fraction * workflow["gate_time_us"],
        opt.gate_points,
    )

    coarse: list[dict[str, Any]] = []
    for n, prep_idx in enumerate(prep_indices, start=1):
        prep_state_lab = prep["u_full_spin"].dag() * prep["states_spin_basis"][int(prep_idx)]
        ramped_state = ramp_detuning(
            dqd,
            prep_state_lab,
            workflow["epsilon_prepare_uev"],
            0.0,
            workflow["epsilon_ramp_us"],
            workflow["ramp_samples"],
            solver_options,
        )
        gate_states = gate_evolver.evolve_many(ramped_state, coarse_gate_tlist_us)
        for gate_idx, gate_state in enumerate(gate_states):
            gate_metrics = bell_readout_metrics(gate_state, dqd, bell_target)
            coarse.append(
                {
                    "prep_idx": int(prep_idx),
                    "prep_time_ns": float(prep["times_ns"][int(prep_idx)]),
                    "gate_hold_us": float(coarse_gate_tlist_us[gate_idx]),
                    "ramped_state": ramped_state,
                    "gate_metrics": gate_metrics,
                }
            )
        print(f"coarse prep {n}/{len(prep_indices)}: {prep['times_ns'][int(prep_idx)]:.4f} ns")

    coarse.sort(key=lambda item: coarse_gate_score(item, opt.score_metric), reverse=True)
    coarse_top = coarse[: max(1, opt.coarse_top_k)]

    gate_step_us = float(coarse_gate_tlist_us[1] - coarse_gate_tlist_us[0])
    refined: list[dict[str, Any]] = []
    for item in coarse_top:
        low = max(0.0, item["gate_hold_us"] - gate_step_us)
        high = min(opt.gate_max_fraction * workflow["gate_time_us"], item["gate_hold_us"] + gate_step_us)
        refine_tlist_us = np.linspace(low, high, opt.refine_points)
        for gate_hold_us in refine_tlist_us:
            gate_state = gate_evolver.evolve_one(item["ramped_state"], float(gate_hold_us))
            readout_state = ramp_detuning(
                dqd,
                gate_state,
                0.0,
                workflow["epsilon_prepare_uev"],
                workflow["epsilon_ramp_us"],
                workflow["ramp_samples"],
                solver_options,
            )
            gate_metrics = bell_readout_metrics(gate_state, dqd, bell_target)
            readout_metrics = bell_readout_metrics(readout_state, dqd, bell_target)
            candidate = {
                "prep_idx": item["prep_idx"],
                "prep_time_ns": item["prep_time_ns"],
                "gate_hold_us": float(gate_hold_us),
                "gate_spin_bell_overlap": gate_metrics["spin_bell_overlap"],
                "gate_spin_phase_optimized_bell_overlap": gate_metrics["spin_phase_optimized_bell_overlap"],
                "gate_spin_phase_optimized_bell_phase_rad": gate_metrics[
                    "spin_phase_optimized_bell_phase_rad"
                ],
                "gate_spin_concurrence": gate_metrics["spin_concurrence"],
                "gate_rr_spin_weight": gate_metrics["rr_spin_weight"],
                "gate_rr_spin_bell_overlap": gate_metrics["rr_spin_bell_overlap"],
                "gate_rr_spin_phase_optimized_bell_overlap": gate_metrics[
                    "rr_spin_phase_optimized_bell_overlap"
                ],
                "gate_rr_spin_phase_optimized_bell_phase_rad": gate_metrics[
                    "rr_spin_phase_optimized_bell_phase_rad"
                ],
                "gate_rr_spin_weighted_bell_overlap": gate_metrics["rr_spin_weighted_bell_overlap"],
                "gate_rr_spin_weighted_phase_optimized_bell_overlap": gate_metrics[
                    "rr_spin_weighted_phase_optimized_bell_overlap"
                ],
                "gate_rr_spin_concurrence": gate_metrics["rr_spin_concurrence"],
                "readout_spin_bell_overlap": readout_metrics["spin_bell_overlap"],
                "readout_spin_phase_optimized_bell_overlap": readout_metrics[
                    "spin_phase_optimized_bell_overlap"
                ],
                "readout_spin_phase_optimized_bell_phase_rad": readout_metrics[
                    "spin_phase_optimized_bell_phase_rad"
                ],
                "readout_spin_concurrence": readout_metrics["spin_concurrence"],
                "readout_rr_spin_weight": readout_metrics["rr_spin_weight"],
                "readout_rr_spin_bell_overlap": readout_metrics["rr_spin_bell_overlap"],
                "readout_rr_spin_phase_optimized_bell_overlap": readout_metrics[
                    "rr_spin_phase_optimized_bell_overlap"
                ],
                "readout_rr_spin_phase_optimized_bell_phase_rad": readout_metrics[
                    "rr_spin_phase_optimized_bell_phase_rad"
                ],
                "readout_rr_spin_weighted_bell_overlap": readout_metrics["rr_spin_weighted_bell_overlap"],
                "readout_rr_spin_weighted_phase_optimized_bell_overlap": readout_metrics[
                    "rr_spin_weighted_phase_optimized_bell_overlap"
                ],
                "readout_rr_spin_concurrence": readout_metrics["rr_spin_concurrence"],
            }
            candidate["score"] = score_candidate(candidate, opt.score_metric)
            refined.append(candidate)

    best = max(refined, key=lambda item: item["score"])
    return {
        "workflow": {
            "epsilon_prepare_uev": workflow["epsilon_prepare_uev"],
            "omega_drive_rad_per_us": workflow["omega_drive"],
            "omega_drive_ghz": workflow["omega_drive_ghz"],
            "epsilon_ramp_us": workflow["epsilon_ramp_us"],
            "ramp_samples": workflow["ramp_samples"],
            "gate_time_us": workflow["gate_time_us"],
        },
        "optimization": asdict(opt),
        "best": best,
        "coarse_top": [
            {
                "prep_time_ns": item["prep_time_ns"],
                "gate_hold_us": item["gate_hold_us"],
                **{f"gate_{k}": v for k, v in item["gate_metrics"].items()},
            }
            for item in coarse_top
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-metric", default=OptimizationConfig.score_metric, choices=[
        "readout_rr_weighted_bell_overlap",
        "readout_rr_weighted_phase_optimized_bell_overlap",
        "readout_rr_bell_overlap",
        "readout_rr_phase_optimized_bell_overlap",
        "readout_spin_bell_overlap",
        "readout_spin_phase_optimized_bell_overlap",
        "gate_rr_weighted_bell_overlap",
        "gate_rr_weighted_phase_optimized_bell_overlap",
        "gate_rr_bell_overlap",
        "gate_rr_phase_optimized_bell_overlap",
        "gate_spin_bell_overlap",
        "gate_spin_phase_optimized_bell_overlap",
    ])
    parser.add_argument("--prep-scan-ns", type=float, default=OptimizationConfig.prep_scan_ns)
    parser.add_argument("--prep-points", type=int, default=OptimizationConfig.prep_points)
    parser.add_argument("--prep-window-ns", type=float, default=OptimizationConfig.prep_window_ns)
    parser.add_argument("--prep-stride-ns", type=float, default=OptimizationConfig.prep_stride_ns)
    parser.add_argument("--gate-max-fraction", type=float, default=OptimizationConfig.gate_max_fraction)
    parser.add_argument("--gate-points", type=int, default=OptimizationConfig.gate_points)
    parser.add_argument("--refine-points", type=int, default=OptimizationConfig.refine_points)
    parser.add_argument("--coarse-top-k", type=int, default=OptimizationConfig.coarse_top_k)
    parser.add_argument("--output-json", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workflow_config = WorkflowConfig()
    opt_config = OptimizationConfig(
        prep_scan_ns=args.prep_scan_ns,
        prep_points=args.prep_points,
        prep_window_ns=args.prep_window_ns,
        prep_stride_ns=args.prep_stride_ns,
        gate_max_fraction=args.gate_max_fraction,
        gate_points=args.gate_points,
        refine_points=args.refine_points,
        coarse_top_k=args.coarse_top_k,
        score_metric=args.score_metric,
    )
    solver_options = {"nsteps": 200000, "atol": 1e-9, "rtol": 1e-9, "store_states": True}
    workflow = build_workflow(workflow_config)
    result = optimize_bell_state(workflow, opt_config, solver_options)

    print(json.dumps(result["workflow"], indent=2))
    print("best")
    print(json.dumps(result["best"], indent=2))
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(f"wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
