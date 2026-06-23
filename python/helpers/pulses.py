from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

Pulse = Callable[[float, dict | None], float]


def square(t_start: float, t_stop: float, value: float = 1.0) -> Pulse:
    """Square pulse: `value` on [t_start, t_stop), 0 otherwise."""

    def f(t: float, args: dict | None = None) -> float:
        return value * (
            np.heaviside(t - t_start, 1.0) - np.heaviside(t - t_stop, 1.0)
        )

    return f


def constant(value: float = 0.0) -> Pulse:
    """Time-uniform coefficient."""

    def f(t: float, args: dict | None = None) -> float:
        return value

    return f


def piecewise(segments: Sequence[tuple[float, float, float]]) -> Pulse:
    """Piecewise-constant pulse from (t_start, t_stop, value) segments."""

    def f(t: float, args: dict | None = None) -> float:
        return sum(
            v * (np.heaviside(t - t0, 1.0) - np.heaviside(t - t1, 1.0))
            for t0, t1, v in segments
        )

    return f


__all__ = ["Pulse", "square", "constant", "piecewise"]
