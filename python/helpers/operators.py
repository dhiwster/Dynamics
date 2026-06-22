from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np


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
