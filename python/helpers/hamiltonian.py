from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .constants import RAD_PER_NS_PER_MICROELECTRONVOLT
from .operators import (
    SingleDQDOperatorsNumpy,
    SingleDQDOperatorsQutip,
    single_dqd_numpy_operators,
    single_dqd_qutip_operators,
)


@dataclass(frozen=True)
class SingleDQDSpectrum:
    phi_p: float
    phi_m: float
    phi_bar: float
    e2: float
    e3: float
    Esigma: float
    Etau: float


def _resolve_coefficients(model: Any) -> tuple[float, float, float]:
    tc = float(model.tc)
    Bz = float(model.Bz)
    if hasattr(model, "bx"):
        bx = float(model.bx)
    elif hasattr(model, "dBx"):
        bx = float(model.dBx)
    else:
        raise AttributeError("single-DQD model must define either 'bx' or 'dBx'")
    return tc, bx, Bz


def _matmul(left: Any, right: Any) -> Any:
    try:
        return left @ right
    except TypeError:
        return left * right


def single_dqd_hamiltonian(model: Any, epsilon: float, ops: Any) -> Any:
    tc, bx, Bz = _resolve_coefficients(model)
    return (
        tc * ops.tx
        + 0.5 * epsilon * ops.tz
        + 0.5 * Bz * ops.sz
        + 0.5 * bx * _matmul(ops.sx, ops.tz)
    )


def single_dqd_numpy_hamiltonian(
    model: Any,
    epsilon: float,
    ops: SingleDQDOperatorsNumpy | None = None,
) -> np.ndarray:
    return single_dqd_hamiltonian(model, epsilon, ops or single_dqd_numpy_operators())


def single_dqd_qutip_hamiltonian(
    model: Any,
    epsilon: float,
    ops: SingleDQDOperatorsQutip | None = None,
) -> Any:
    return single_dqd_hamiltonian(model, epsilon, ops or single_dqd_qutip_operators())


def single_dqd_spectrum(model: Any) -> SingleDQDSpectrum:
    tc, bx, Bz = _resolve_coefficients(model)

    phi_p = np.arctan2(bx, 2 * tc + Bz)
    if phi_p < 0:
        phi_p += np.pi

    phi_m = np.arctan2(bx, 2 * tc - Bz)
    if phi_m < 0:
        phi_m += np.pi

    e2 = 0.5 * np.sqrt((2 * tc - Bz) ** 2 + bx**2)
    e3 = 0.5 * np.sqrt((2 * tc + Bz) ** 2 + bx**2)
    return SingleDQDSpectrum(
        phi_p=phi_p,
        phi_m=phi_m,
        phi_bar=0.5 * (phi_p + phi_m),
        e2=e2,
        e3=e3,
        Esigma=e3 - e2,
        Etau=e3 + e2,
    )


def single_dqd_eigensystem(model: Any, epsilon: float) -> tuple[np.ndarray, np.ndarray]:
    energies, vectors = np.linalg.eigh(single_dqd_numpy_hamiltonian(model, epsilon))
    order = np.argsort(energies)
    return energies[order], vectors[:, order]


def single_dqd_qubit_splitting(model: Any, epsilon: float) -> float:
    energies, _ = single_dqd_eigensystem(model, epsilon)
    return float((energies[1] - energies[0]) * RAD_PER_NS_PER_MICROELECTRONVOLT)


def single_dqd_tau_z_matrix_element(model: Any, epsilon: float = 0.0) -> float:
    _, vectors = single_dqd_eigensystem(model, epsilon)
    g = vectors[:, 0]
    e = vectors[:, 1]
    tz = single_dqd_numpy_operators().tz
    return float(abs(np.vdot(g, tz @ e)))
