from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np

TWO_PI = 2 * np.pi
HBAR_UEV_NS = 0.6582119569
HBAR_UEV_US = HBAR_UEV_NS / 1000.0


@dataclass(frozen=True)
class SingleDQDOperatorsNumpy:
    sx: np.ndarray
    sy: np.ndarray
    sz: np.ndarray
    tx: np.ndarray
    ty: np.ndarray
    tz: np.ndarray


@dataclass(frozen=True)
class SingleDQDOperatorsQutip:
    sx: Any
    sy: Any
    sz: Any
    sm: Any
    sp: Any
    tx: Any
    ty: Any
    tz: Any


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


@lru_cache(maxsize=1)
def single_dqd_numpy_operators() -> SingleDQDOperatorsNumpy:
    ident = np.eye(2, dtype=complex)
    sx_1q = np.array([[0, 1], [1, 0]], dtype=complex)
    sy_1q = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz_1q = np.array([[1, 0], [0, -1]], dtype=complex)

    return SingleDQDOperatorsNumpy(
        sx=np.kron(sx_1q, ident),
        sy=np.kron(sy_1q, ident),
        sz=np.kron(sz_1q, ident),
        tx=np.kron(ident, sx_1q),
        ty=np.kron(ident, sy_1q),
        tz=np.kron(ident, sz_1q),
    )


@lru_cache(maxsize=1)
def single_dqd_qutip_operators() -> SingleDQDOperatorsQutip:
    from qutip import qeye, sigmam, sigmap, sigmax, sigmay, sigmaz, tensor

    return SingleDQDOperatorsQutip(
        sx=tensor(sigmax(), qeye(2)),
        sy=tensor(sigmay(), qeye(2)),
        sz=tensor(sigmaz(), qeye(2)),
        sm=tensor(sigmam(), qeye(2)),
        sp=tensor(sigmap(), qeye(2)),
        tx=tensor(qeye(2), sigmax()),
        ty=tensor(qeye(2), sigmay()),
        tz=tensor(qeye(2), sigmaz()),
    )


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

    e2 = 0.5 * np.sqrt((2 * tc - Bz) ** 2 + bx ** 2)
    e3 = 0.5 * np.sqrt((2 * tc + Bz) ** 2 + bx ** 2)
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
    return float((energies[1] - energies[0]) / HBAR_UEV_NS)


def single_dqd_tau_z_matrix_element(model: Any, epsilon: float = 0.0) -> float:
    _, vectors = single_dqd_eigensystem(model, epsilon)
    g = vectors[:, 0]
    e = vectors[:, 1]
    tz = single_dqd_numpy_operators().tz
    return float(abs(np.vdot(g, tz @ e)))
