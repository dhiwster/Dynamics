from __future__ import annotations

import numpy as np

from .constants import HBAR_UEV_NS
from .gates import (
    _H0_ns,
    _edsr1_op,
    _edsr2_op,
    _eps1_op,
    _eps2_op,
    zrot_delta_freq_GHz,
)
from .system import DQDsystem

hbar_ns = HBAR_UEV_NS


class DQDSequenceCompiler:
    """
    Chain gate steps with virtual Z-rotation and automatic phase tracking.

    Z-rotations are virtual (zero duration, no pulse): the rotation angle
    is added directly to that DQD's phase accumulator. Every subsequent
    X-rotation on the same DQD applies its EDSR drive with a phase offset
    equal to the accumulated value, implementing the rotation through frame
    tracking rather than a physical ε pulse.

    Phase accumulators are also updated automatically whenever a DQD sits
    idle under ε_idle during a neighbour's X-rotation:

        Δφ_idle = 2π · ΔEσ(ε_idle)/ħ · t_gate

    iSWAP holds both DQDs at ε=0, so no extra phase accumulates there.

    Parameters
    ----------
    dqd : DQDsystem
    dt  : float — time resolution of the master tlist (ns).
    """

    def __init__(self, dqd: DQDsystem, dt: float = 0.01):
        self._sys        = dqd
        self._dt         = dt
        self._steps: list = []
        self._phi_bar    = dqd.phi_bar
        self._Esigma_GHz = dqd.Esigma / hbar_ns
        self._eps_idle   = dqd.epsilon_idle
        self.phase_accumulator = [0.0, 0.0]

    @property
    def phase_accumulator_dqd1(self) -> float:
        return self.phase_accumulator[0]

    @property
    def phase_accumulator_dqd2(self) -> float:
        return self.phase_accumulator[1]

    def reset_phases(self) -> 'DQDSequenceCompiler':
        self.phase_accumulator = [0.0, 0.0]
        return self

    def add_iswap(self) -> 'DQDSequenceCompiler':
        duration = self._sys.iSWAP_gate_time() * 1e3
        self._steps.append({'type': 'iswap', 'duration': duration})
        return self

    def add_xrot(self, target: int, angle: float) -> 'DQDSequenceCompiler':
        if target not in (1, 2):
            raise ValueError("target must be 1 or 2")

        duration  = abs(angle) * hbar_ns / (self._sys.Vac0 * np.sin(self._phi_bar))
        phi_drive = self.phase_accumulator[target - 1]

        idle = 3 - target
        delta_f = zrot_delta_freq_GHz(self._sys, self._eps_idle)
        self.phase_accumulator[idle - 1] += 2 * np.pi * delta_f * duration

        self._steps.append({
            'type':      'xrot',
            'target':    target,
            'angle':     angle,
            'duration':  duration,
            'phi_drive': phi_drive,
        })
        return self

    def add_zrot(self, target: int, angle: float) -> 'DQDSequenceCompiler':
        if target not in (1, 2):
            raise ValueError("target must be 1 or 2")
        self.phase_accumulator[target - 1] += angle
        return self

    @property
    def step_durations(self) -> list:
        return [s['duration'] for s in self._steps]

    @property
    def total_time(self) -> float:
        return sum(s['duration'] for s in self._steps)

    def __len__(self) -> int:
        return len(self._steps)

    def __repr__(self) -> str:
        phi1_deg = np.degrees(self.phase_accumulator[0])
        phi2_deg = np.degrees(self.phase_accumulator[1])
        lines = [
            f"DQDSequenceCompiler  ({len(self)} steps, {self.total_time:.4f} ns total)",
            f"  Phase accumulators:  DQD1 = {phi1_deg:.3f}°   DQD2 = {phi2_deg:.3f}°",
        ]
        offsets = [0.0]
        for s in self._steps:
            offsets.append(offsets[-1] + s['duration'])
        for i, s in enumerate(self._steps):
            t0, t1 = offsets[i], offsets[i + 1]
            if s['type'] == 'iswap':
                desc = "iSWAP  (ε=0 both DQDs)"
            else:
                desc = (f"Rx({np.degrees(s['angle']):.1f}°)  target={s['target']}"
                        f"  φ_drive={np.degrees(s['phi_drive']):.2f}°"
                        f"  Vac0={self._sys.Vac0} μeV")
            lines.append(f"  [{t0:.3f} – {t1:.3f} ns]  {desc}")
        return "\n".join(lines)

    def compile(self) -> tuple:
        """
        Compile the gate sequence to a QuTiP time-dependent Hamiltonian.

        Returns
        -------
        H_td  : list compatible with qt.mesolve and qt.propagator
        tlist : np.ndarray — master time axis in ns
        """
        if not self._steps:
            raise ValueError("No physical gates added — call add_iswap / add_xrot first")

        dqd        = self._sys
        eps_idle   = self._eps_idle
        Esigma_GHz = self._Esigma_GHz

        offsets = [0.0]
        for s in self._steps:
            offsets.append(offsets[-1] + s['duration'])
        T = offsets[-1]

        tlist   = np.arange(0.0, T + self._dt, self._dt)
        H_drift = _H0_ns(dqd)
        terms: list = []

        for i, step in enumerate(self._steps):
            t0    = offsets[i]
            t1    = offsets[i + 1]
            gtype = step['type']

            if gtype == 'iswap':
                pass

            elif gtype == 'xrot':
                tgt       = step['target']
                phi_drive = step['phi_drive']
                H_drv = _edsr1_op(dqd) if tgt == 1 else _edsr2_op(dqd)
                H_idl = _eps2_op(dqd)  if tgt == 1 else _eps1_op(dqd)

                def _edsr(t, args=None, _t0=t0, _t1=t1, _w=Esigma_GHz, _phi=phi_drive):
                    return np.cos(_w * t + _phi) if _t0 <= t < _t1 else 0.0

                def _idle_x(t, args=None, _t0=t0, _t1=t1, _e=eps_idle):
                    return _e if _t0 <= t < _t1 else 0.0

                terms.append([H_drv, _edsr])
                terms.append([H_idl, _idle_x])

        return [H_drift] + terms, tlist
