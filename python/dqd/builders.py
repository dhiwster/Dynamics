from __future__ import annotations

from hamiltonians.pulses import Pulse, constant


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


__all__ = ["build_H"]
