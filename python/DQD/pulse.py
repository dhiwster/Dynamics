from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

Pulse = Callable[[float, dict | None], float]


# ---------------------------------------------------------------------------
# Pulse shapes
# ---------------------------------------------------------------------------

def square(t_start: float, t_stop: float, value: float = 1.0) -> Pulse:
    """
    Square pulse: returns `value` during [t_start, t_stop), 0 otherwise.

    Parameters
    ----------
    t_start, t_stop : pulse window (same units as the simulation time)
    value           : amplitude of the pulse (e.g. tc in μeV, or epsilon in μeV)

    Returns
    -------
    f(t, args=None) : callable for use as a qutip time-dependent coefficient
    """
    def f(t: float, args: dict | None = None) -> float:
        return value * (np.heaviside(t - t_start, 1.0)
                        - np.heaviside(t - t_stop,  1.0))
    return f


def constant(value: float = 0.0) -> Pulse:
    """Time-uniform coefficient (e.g. tc held at a fixed value throughout)."""
    def f(t: float, args: dict | None = None) -> float:
        return value
    return f


def piecewise(segments: Sequence[tuple[float, float, float]]) -> Pulse:
    """
    Piecewise-constant pulse built from a list of (t_start, t_stop, value) segments.
    Overlapping segments are summed.

    Example
    -------
    tc_pulse = piecewise([(0, 10, 80), (20, 30, 40)])
    """
    def f(t: float, args: dict | None = None) -> float:
        return sum(v * (np.heaviside(t - t0, 1.0) - np.heaviside(t - t1, 1.0))
                   for t0, t1, v in segments)
    return f

# ---------------------------------------------------------------------------
# Hamiltonian builder
# ---------------------------------------------------------------------------

def _delta(base: float, pulse: Pulse) -> Pulse:
    """Convert an absolute control pulse into an offset from the static value."""

    def f(t: float, args: dict | None = None) -> float:
        return pulse(t, args) - base

    return f


def build_H(
    dqd,
    tc1: Pulse | None = None,
    tc2: Pulse | None = None,
    eps1: Pulse | None = None,
    eps2: Pulse | None = None,
    same_for_both: bool = True,
):
    """
    Build a qutip time-dependent Hamiltonian list from a DQDsystem and pulse callables.

    H = [H_static,
         [H_tc1_op,  tc1(t)],
         [H_tc2_op,  tc2(t)],
         [H_eps1_op, eps1(t)],
         [H_eps2_op, eps2(t)]]

    Parameters
    ----------
    dqd            : DQDsystem instance
    tc1            : callable f(t, args) for DQD1 tunneling (μeV); uses dqd.tc if None
    tc2            : callable f(t, args) for DQD2 tunneling (μeV); defaults to tc1
                     when same_for_both=True
    eps1           : callable f(t, args) for DQD1 detuning (μeV); uses 0 if None
    eps2           : callable f(t, args) for DQD2 detuning (μeV); defaults to eps1
                     when same_for_both=True
    same_for_both  : if True, tc2 mirrors tc1 and eps2 mirrors eps1 when not provided

    Returns
    -------
    H : list in qutip mesolve/propagator format
    """
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
        [dqd.H_tc1_op,  _delta(dqd.tc, tc1)],
        [dqd.H_tc2_op,  _delta(dqd.tc, tc2)],
        [dqd.H_eps1_op, eps1],
        [dqd.H_eps2_op, eps2],
    ]
