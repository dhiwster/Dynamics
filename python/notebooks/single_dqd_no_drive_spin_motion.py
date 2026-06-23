"""
Standalone 4x4 single-DQD verification of no-drive spin motion.

This script isolates the effect seen in the resonator workflow notebook:
the DQD1 spin can oscillate even when the explicit drive amplitude is zero.

Basis order:

    |up,R>, |up,L>, |down,R>, |down,L>

Hamiltonian, in micro-eV:

    H = tc tau_x
      + epsilon / 2 tau_z
      + Bz / 2 sigma_z
      + bx / 2 sigma_x tau_z

There is intentionally no cavity, no second DQD, and no time-dependent drive.
The only dynamics come from free evolution under this 4x4 Hamiltonian.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from helpers import (
    RAD_PER_NS_PER_MICROELECTRONVOLT,
    TWO_PI,
    single_dqd_eigensystem,
    single_dqd_hamiltonian,
    single_dqd_numpy_operators,
)


OPS = single_dqd_numpy_operators()


@dataclass(frozen=True)
class SingleDQDParams:
    tc: float = 80.0
    bx: float = 20.0
    Bz: float = 24.0
    epsilon: float = 0.0


def ghz_to_uev(ghz: float) -> float:
    """Convert a frequency input in GHz to the model detuning unit, micro-eV."""
    return ghz * TWO_PI / RAD_PER_NS_PER_MICROELECTRONVOLT


def static_hamiltonian(params: SingleDQDParams) -> np.ndarray:
    return single_dqd_hamiltonian(params, params.epsilon, OPS)


def basis_state(label: str) -> np.ndarray:
    states = {
        "up_R": np.array([1.0, 0.0, 0.0, 0.0], dtype=complex),
        "up_L": np.array([0.0, 1.0, 0.0, 0.0], dtype=complex),
        "down_R": np.array([0.0, 0.0, 1.0, 0.0], dtype=complex),
        "down_L": np.array([0.0, 0.0, 0.0, 1.0], dtype=complex),
    }
    try:
        return states[label].copy()
    except KeyError as exc:
        raise ValueError(f"Unknown initial state {label!r}") from exc


def evolve_state(psi0: np.ndarray, hamiltonian_uev: np.ndarray, t_eval_ns: np.ndarray) -> np.ndarray:
    """Return states as columns, evolved under a constant 4x4 Hamiltonian."""

    generator = -1j * hamiltonian_uev * RAD_PER_NS_PER_MICROELECTRONVOLT

    def rhs(_t: float, psi: np.ndarray) -> np.ndarray:
        return generator @ psi

    result = solve_ivp(
        rhs,
        (float(t_eval_ns[0]), float(t_eval_ns[-1])),
        psi0,
        t_eval=t_eval_ns,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.y


def spin_populations(states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Spin populations after tracing out charge by direct basis summation."""
    p_up = np.abs(states[0]) ** 2 + np.abs(states[1]) ** 2
    p_down = np.abs(states[2]) ** 2 + np.abs(states[3]) ** 2
    return p_up.real, p_down.real


def charge_right_probability(states: np.ndarray) -> np.ndarray:
    return (np.abs(states[0]) ** 2 + np.abs(states[2]) ** 2).real


def eigenstate_overlaps(params: SingleDQDParams, psi0: np.ndarray) -> np.ndarray:
    _, vectors = single_dqd_eigensystem(params, params.epsilon)
    return np.abs(vectors.conj().T @ psi0) ** 2


def run_case(params: SingleDQDParams, psi0: np.ndarray, t_eval_ns: np.ndarray) -> dict[str, np.ndarray]:
    states = evolve_state(psi0, static_hamiltonian(params), t_eval_ns)
    p_up, p_down = spin_populations(states)
    p_right = charge_right_probability(states)
    overlaps = eigenstate_overlaps(params, psi0)
    return {
        "states": states,
        "p_up": p_up,
        "p_down": p_down,
        "p_right": p_right,
        "overlaps": overlaps,
    }


def make_plot(
    t_eval_ns: np.ndarray,
    full_case: dict[str, np.ndarray],
    no_bx_case: dict[str, np.ndarray],
    title_suffix: str,
) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(t_eval_ns, full_case["p_up"], label="spin up, bx on")
    axes[0].plot(t_eval_ns, full_case["p_down"], label="spin down, bx on")
    axes[0].plot(t_eval_ns, no_bx_case["p_up"], "--", label="spin up, bx = 0")
    axes[0].plot(t_eval_ns, no_bx_case["p_down"], "--", label="spin down, bx = 0")
    axes[0].set_xlabel("time (ns)")
    axes[0].set_ylabel("population")
    axes[0].set_title(f"Spin populations {title_suffix}")
    axes[0].legend()
    axes[0].grid(True)

    axes[1].plot(t_eval_ns, full_case["p_right"], label="charge-right, bx on")
    axes[1].plot(t_eval_ns, no_bx_case["p_right"], "--", label="charge-right, bx = 0")
    axes[1].set_xlabel("time (ns)")
    axes[1].set_ylabel("probability")
    axes[1].set_title(f"Charge-right occupancy {title_suffix}")
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tc", type=float, default=80.0, help="Tunnel coupling in micro-eV.")
    parser.add_argument("--bx", type=float, default=20.0, help="Transverse spin-charge coupling in micro-eV.")
    parser.add_argument("--Bz", type=float, default=24.0, help="Zeeman splitting in micro-eV.")
    parser.add_argument(
        "--epsilon-ghz",
        type=float,
        default=50.0,
        help="Detuning specified as a frequency in GHz; converted internally to micro-eV.",
    )
    parser.add_argument(
        "--initial-state",
        choices=["up_R", "up_L", "down_R", "down_L"],
        default="down_R",
        help="Bare 4x4 basis state to evolve.",
    )
    parser.add_argument("--t-end-ns", type=float, default=2.0, help="Total simulation time in ns.")
    parser.add_argument("--num-points", type=int, default=2001, help="Number of time samples.")
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="Optional path to save the plot instead of only displaying it.",
    )
    args = parser.parse_args()

    params = SingleDQDParams(
        tc=args.tc,
        bx=args.bx,
        Bz=args.Bz,
        epsilon=ghz_to_uev(args.epsilon_ghz),
    )
    params_no_bx = SingleDQDParams(
        tc=args.tc,
        bx=0.0,
        Bz=args.Bz,
        epsilon=ghz_to_uev(args.epsilon_ghz),
    )

    psi0 = basis_state(args.initial_state)
    t_eval_ns = np.linspace(0.0, args.t_end_ns, args.num_points)

    full_case = run_case(params, psi0, t_eval_ns)
    no_bx_case = run_case(params_no_bx, psi0, t_eval_ns)

    print("Initial bare state:", args.initial_state)
    print("Parameters with bx on:", params)
    print("Parameters with bx off:", params_no_bx)
    print("Initial-state overlaps with finite-detuning eigenstates (bx on):", full_case["overlaps"])
    print("Initial-state overlaps with finite-detuning eigenstates (bx off):", no_bx_case["overlaps"])
    print(
        "Spin-up range with bx on:",
        float(full_case["p_up"].min()),
        float(full_case["p_up"].max()),
    )
    print(
        "Spin-up range with bx off:",
        float(no_bx_case["p_up"].min()),
        float(no_bx_case["p_up"].max()),
    )

    title_suffix = f"from |{args.initial_state.replace('_', ', ')}>"
    fig = make_plot(t_eval_ns, full_case, no_bx_case, title_suffix)
    if args.save is not None:
        fig.savefig(args.save, dpi=150, bbox_inches="tight")
        print("Saved plot to:", args.save)

    if "agg" not in matplotlib.get_backend().lower():
        plt.show()


if __name__ == "__main__":
    main()
