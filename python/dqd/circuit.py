"""
circuit.py — DQD gate Hamiltonian builders and sequence compiler.

Units
-----
  Time       : nanoseconds (ns)
  Energy     : μeV
  Frequency  : GHz = rad/ns
  Conversion : E[μeV] / ħ[μeV·ns] → ω[GHz]

Tunable control parameter (ONLY dynamic knob):
  ε (epsilon) — detuning, enters H as (ε/2) τ_z

Fixed (embedded in DQDsystem, never pulsed):
  Bz — Zeeman splitting
  bx — micro-magnet / Rashba SOI
  tc — tunneling coupling
"""

from __future__ import annotations

import numpy as np
import qutip as qt

from .system import DQDsystem

# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

hbar_ns  = 0.6582119569   # ħ in μeV·ns
_MHz2GHz = 1e-3

# ---------------------------------------------------------------------------
# Internal operator helpers  (all take a DQDsystem instance)
# ---------------------------------------------------------------------------

def _H0_ns(dqd: DQDsystem) -> qt.Qobj:
    """Full static Hamiltonian at ε=0 in GHz."""
    return dqd.H_static * _MHz2GHz


def _eps1_op(dqd: DQDsystem) -> qt.Qobj:
    """DQD1 detuning operator [GHz/μeV]: τ_z1 / (2ħ)."""
    return qt.tensor(qt.qeye(dqd.photon_max), dqd.tz1) / (2 * hbar_ns)


def _eps2_op(dqd: DQDsystem) -> qt.Qobj:
    """DQD2 detuning operator [GHz/μeV]: τ_z2 / (2ħ)."""
    return qt.tensor(qt.qeye(dqd.photon_max), dqd.tz2) / (2 * hbar_ns)


def _edsr1_op(dqd: DQDsystem) -> qt.Qobj:
    """DQD1 EDSR drive operator [GHz]: Vac0 · τ_z1 / ħ."""
    return dqd.H_edsr1_amplitude / hbar_ns


def _edsr2_op(dqd: DQDsystem) -> qt.Qobj:
    """DQD2 EDSR drive operator [GHz]: Vac0 · τ_z2 / ħ."""
    return dqd.H_edsr2_amplitude / hbar_ns


# ---------------------------------------------------------------------------
# Calibration helpers
# ---------------------------------------------------------------------------

def zrot_delta_freq_GHz(dqd: DQDsystem, eps_amp: float) -> float:
    """
    Qubit eigenfrequency shift ΔEσ/ħ [GHz] caused by a DC ε detuning.

    Uses the single-DQD Hamiltonian:

        H(ε) = tc·τ_x + (Bz/2)·σ_z + (bx/2)·σ_x·τ_z + (ε/2)·τ_z

    Parameters
    ----------
    dqd     : DQDsystem instance
    eps_amp : μeV, detuning amplitude

    Returns
    -------
    delta_f : GHz (always positive for eps_amp != 0)
    """
    tc, bx, Bz = dqd.tc, dqd.bx, dqd.Bz
    tx, sx, sz, tz = dqd.tx, dqd.sx, dqd.sz, dqd.tz

    H_base = tc * tx + Bz / 2 * sz + bx / 2 * sx * tz
    H_det  = H_base + eps_amp / 2 * tz

    e_base = np.sort(H_base.eigenenergies())
    e_det  = np.sort(H_det.eigenenergies())

    return abs((e_det[1] - e_det[0]) - (e_base[1] - e_base[0])) / hbar_ns


def zrot_time_ns(dqd: DQDsystem, eps_amp: float, angle: float) -> float:
    """
    Pulse duration [ns] for Rz(angle) via a DC ε detuning.

    Parameters
    ----------
    dqd     : DQDsystem instance
    eps_amp : μeV, detuning amplitude
    angle   : rad, target rotation angle

    Returns
    -------
    t_ns : ns
    """
    delta_f = zrot_delta_freq_GHz(dqd, eps_amp)
    if delta_f == 0:
        raise ValueError("eps_amp=0 gives zero qubit frequency shift")
    return abs(angle) / (2 * np.pi * delta_f)


# ---------------------------------------------------------------------------
# Gate Hamiltonian builders  (standalone; take dqd as first argument)
# ---------------------------------------------------------------------------

def iswap_H(dqd: DQDsystem) -> list:
    """
    iSWAP Hamiltonian in GHz at the ε=0 sweet spot.

    Returns [H_const] for qt.mesolve / qt.propagator.
    Gate duration: dqd.iSWAP_gate_time() * 1e3 ns.
    """
    return [_H0_ns(dqd)]


def xrot_H(
    dqd: DQDsystem,
    target: int,
    t_start: float,
    t_end: float,
    eps_idle: float | None = None,
) -> list:
    """
    EDSR X-rotation Hamiltonian for ``target`` DQD in GHz / ns.

        H = H0 + [Vac0·τ_z_target/ħ, cos(Eσ/ħ·t)] + [H_eps_idle, ε_idle]

    Parameters
    ----------
    dqd      : DQDsystem instance
    target   : 1 or 2
    t_start  : pulse start (ns)
    t_end    : pulse end   (ns)
    eps_idle : μeV, ε on the idle DQD; defaults to dqd.epsilon_idle
    """
    if target not in (1, 2):
        raise ValueError("target must be 1 or 2")
    eps_idle    = dqd.epsilon_idle if eps_idle is None else eps_idle
    Esigma_GHz  = dqd.Esigma / hbar_ns
    H0          = _H0_ns(dqd)
    H_drv, H_idl = (_edsr1_op(dqd), _eps2_op(dqd)) if target == 1 else (_edsr2_op(dqd), _eps1_op(dqd))

    def _edsr(t, args=None):
        return np.cos(Esigma_GHz * t) if t_start <= t < t_end else 0.0

    def _idle(t, args=None):
        return eps_idle if not (t_start <= t < t_end) else 0.0

    return [H0, [H_drv, _edsr], [H_idl, _idle]]


def zrot_H(
    dqd: DQDsystem,
    target: int,
    t_start: float,
    t_end: float,
    eps_amp: float,
    eps_idle: float | None = None,
) -> list:
    """
    Z-rotation Hamiltonian via ε detuning pulse on ``target`` DQD.

    Parameters
    ----------
    dqd      : DQDsystem instance
    target   : 1 or 2
    t_start  : pulse start (ns)
    t_end    : pulse end   (ns)
    eps_amp  : μeV, detuning amplitude (sign controls Rz direction)
    eps_idle : μeV, ε on the idle DQD; defaults to dqd.epsilon_idle
    """
    if target not in (1, 2):
        raise ValueError("target must be 1 or 2")
    eps_idle = dqd.epsilon_idle if eps_idle is None else eps_idle
    H0       = _H0_ns(dqd)
    H_tgt, H_idl = (_eps1_op(dqd), _eps2_op(dqd)) if target == 1 else (_eps2_op(dqd), _eps1_op(dqd))

    def _eps_pulse(t, args=None):
        return eps_amp if t_start <= t < t_end else 0.0

    def _idle_full(t, args=None):
        return eps_idle if not (t_start <= t < t_end) else 0.0

    return [H0, [H_tgt, _eps_pulse], [H_idl, _idle_full]]


# ---------------------------------------------------------------------------
# Eigenstate preparation
# ---------------------------------------------------------------------------

def vacuum_eigenstates(dqd: DQDsystem) -> list:
    """
    Ground and first excited eigenstates of the single-DQD H₀ (no photon).

    Parameters
    ----------
    dqd : DQDsystem instance

    Returns
    -------
    states : list of 4-dim Qobj kets, sorted by energy
    """
    H0_single = dqd.tc * dqd.tx + dqd.Bz / 2 * dqd.sz + dqd.bx / 2 * dqd.sx * dqd.tz
    _, states = H0_single.eigenstates(sort='low')
    return list(states)


def initial_full_state(
    dqd: DQDsystem,
    dqd1_eigenstate_idx: int = 0,
    dqd2_eigenstate_idx: int = 1,
) -> qt.Qobj:
    """
    Build |ψ₀⟩ = |n=0⟩_photon ⊗ |φ_i⟩_DQD1 ⊗ |φ_j⟩_DQD2.

    Parameters
    ----------
    dqd                  : DQDsystem instance
    dqd1_eigenstate_idx  : 0 = ground, 1 = first excited, ...
    dqd2_eigenstate_idx  : same

    Returns
    -------
    psi0 : Qobj ket, dims [[photon_max,2,2,2,2],[1,1,1,1,1]]
    """
    states = vacuum_eigenstates(dqd)
    return qt.tensor(
        qt.basis(dqd.photon_max, 0),
        states[dqd1_eigenstate_idx],
        states[dqd2_eigenstate_idx],
    )


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(dqd: DQDsystem) -> None:
    """Print key derived parameters and gate times for the given DQDsystem."""
    phi_bar      = dqd.phi_bar
    Esigma_GHz   = dqd.Esigma / hbar_ns
    g_sigma_GHz  = dqd.g_sigma * _MHz2GHz
    d_sigma_GHz  = dqd.d_sigma * _MHz2GHz
    tg_iSWAP     = dqd.iSWAP_gate_time() * 1e3
    t_x180       = np.pi * hbar_ns / (dqd.Vac0 * np.sin(phi_bar))
    r_s, _       = dqd.dispersive_ratios()

    print(f"DQD:        tc={dqd.tc} μeV,  Bz={dqd.Bz} μeV (fixed),  bx={dqd.bx} μeV (fixed)")
    print(f"Cavity:     wc={dqd.wc/2/np.pi:.0f} MHz,  gc={dqd.gc/2/np.pi:.1f} MHz")
    print(f"Qubit:      Eσ={dqd.Esigma:.4f} μeV  ({Esigma_GHz:.3f} GHz)")
    print(f"φ̄:          {np.degrees(phi_bar):.3f}°,  sin(φ̄)={np.sin(phi_bar):.6f}")
    print(f"Coupling:   g_σ={g_sigma_GHz*1e3:.4f} MHz,  Δ_σ={d_sigma_GHz*1e3:.2f} MHz,  g/Δ={r_s:.5f}")
    print(f"tg(iSWAP):  {tg_iSWAP:.4f} ns")
    print(f"t(Rx90):    {t_x180/2:.4f} ns,  t(Rx180): {t_x180:.4f} ns")
    print(f"Vac0:       {dqd.Vac0} μeV  (Rabi={dqd.Vac0/hbar_ns*np.sin(phi_bar)*1e3:.3f} MHz)")
    print(f"ε_idle:     {dqd.epsilon_idle:.1f} μeV")


# ---------------------------------------------------------------------------
# Sequence Compiler
# ---------------------------------------------------------------------------

class DQDSequenceCompiler:
    """
    Chain gate steps with virtual Z-rotation and automatic phase tracking.

    Z-rotations are virtual (zero duration, no pulse): the rotation angle
    is added directly to that DQD's phase accumulator. Every subsequent
    X-rotation on the same DQD applies its EDSR drive with a phase offset
    equal to the accumulated value, implementing the rotation through frame
    tracking rather than a physical ε pulse.

    Phase accumulators are also updated automatically whenever a DQD sits
    idle under ε_idle during a neighbour's X-rotation, because the detuning
    shifts the qubit eigenfrequency and causes phase accumulation in the
    rotating frame:

        Δφ_idle = 2π · ΔEσ(ε_idle)/ħ · t_gate

    iSWAP holds both DQDs at ε=0, so no extra phase accumulates there.

    Hardware constraints
    --------------------
    * ε is the ONLY dynamically adjustable parameter.
    * Bz, bx, tc are embedded in H_drift and never touched.

    Parameters
    ----------
    dqd : DQDsystem — provides all physical parameters and pre-built operators.
    dt  : float     — time resolution of the master tlist (ns).

    Usage
    -----
    dqd = DQDsystem(tc=100, bx=2, Bz=20, Vac0=1.0, ...)
    seq = DQDSequenceCompiler(dqd)
    seq.add_zrot(target=1, angle=np.pi / 4)   # virtual — no pulse emitted
    seq.add_xrot(target=1, angle=np.pi / 2)   # drive phase shifted by π/4
    seq.add_iswap()
    H_td, tlist = seq.compile()
    """

    def __init__(self, dqd: DQDsystem, dt: float = 0.01):
        self._sys        = dqd
        self._dt         = dt
        self._steps: list = []
        self._phi_bar    = dqd.phi_bar
        self._Esigma_GHz = dqd.Esigma / hbar_ns
        self._eps_idle   = dqd.epsilon_idle
        # One phase accumulator per DQD (index 0 = DQD1, index 1 = DQD2).
        # Updated by add_zrot (explicit virtual rotation) and by add_xrot
        # (automatic idle-detuning phase on the non-driven DQD).
        self.phase_accumulator = [0.0, 0.0]

    # ------------------------------------------------------------------
    # Phase accumulator helpers
    # ------------------------------------------------------------------

    @property
    def phase_accumulator_dqd1(self) -> float:
        """Current accumulated virtual Z phase for DQD1 (radians)."""
        return self.phase_accumulator[0]

    @property
    def phase_accumulator_dqd2(self) -> float:
        """Current accumulated virtual Z phase for DQD2 (radians)."""
        return self.phase_accumulator[1]

    def reset_phases(self) -> 'DQDSequenceCompiler':
        """Reset both phase accumulators to zero."""
        self.phase_accumulator = [0.0, 0.0]
        return self

    # ------------------------------------------------------------------
    # Gate adders  (all return self for optional chaining)
    # ------------------------------------------------------------------

    def add_iswap(self) -> 'DQDSequenceCompiler':
        """
        iSWAP: both DQDs at ε=0 for iSWAP_gate_time() ns.

        No phase accumulation — at the sweet spot the qubit frequency
        equals the drive frequency so the rotating-frame detuning is zero.
        """
        duration = self._sys.iSWAP_gate_time() * 1e3   # μs → ns
        self._steps.append({'type': 'iswap', 'duration': duration})
        return self

    def add_xrot(self, target: int, angle: float) -> 'DQDSequenceCompiler':
        """
        EDSR X-rotation on ``target`` DQD (1 or 2).

        The EDSR drive phase is offset by the target DQD's current phase
        accumulator value:

            coefficient(t) = cos(Eσ/ħ · t + φ_accum)

        After the gate the idle DQD's accumulator is incremented by the
        phase it accumulates under ε_idle:

            Δφ_idle = 2π · ΔEσ(ε_idle)/ħ · t_gate

        Pulse duration: t = |θ| · ħ / (Vac0 · sin(φ̄))

        Parameters
        ----------
        target : 1 or 2
        angle  : float, rotation angle in radians
        """
        if target not in (1, 2):
            raise ValueError("target must be 1 or 2")

        duration  = abs(angle) * hbar_ns / (self._sys.Vac0 * np.sin(self._phi_bar))
        phi_drive = self.phase_accumulator[target - 1]   # snapshot before gate

        # Idle DQD accumulates phase from ε_idle detuning during this gate
        idle = 3 - target   # 1 → 2,  2 → 1
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
        """
        Virtual Z-rotation on ``target`` DQD (1 or 2).

        No physical pulse is emitted and no time is added to the sequence.
        The angle is added to the target DQD's phase accumulator; the next
        X-rotation on that DQD will shift its EDSR drive phase accordingly,
        implementing Rz through frame tracking.

        Parameters
        ----------
        target : 1 or 2
        angle  : float, rotation angle in radians
        """
        if target not in (1, 2):
            raise ValueError("target must be 1 or 2")
        self.phase_accumulator[target - 1] += angle
        return self

    # ------------------------------------------------------------------
    # Inspection
    # ------------------------------------------------------------------

    @property
    def step_durations(self) -> list:
        """Per-step durations (ns)  — virtual Rz steps have zero duration."""
        return [s['duration'] for s in self._steps]

    @property
    def total_time(self) -> float:
        """Total physical sequence duration (ns)."""
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

    # ------------------------------------------------------------------
    # Compiler
    # ------------------------------------------------------------------

    def compile(self) -> tuple:
        """
        Compile the gate sequence to a QuTiP time-dependent Hamiltonian.

        Virtual Z-rotations do not appear in H_td; their effect is already
        baked into the ``phi_drive`` field of each xrot step recorded during
        the add_xrot call.

        Returns
        -------
        H_td  : list — [H_drift, [H_ctrl, f(t)], ...]
                  Compatible with qt.mesolve and qt.propagator.
        tlist : np.ndarray — master time axis in ns, step = dt.
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
                pass   # ε=0: no control terms needed

            elif gtype == 'xrot':
                tgt       = step['target']
                phi_drive = step['phi_drive']
                H_drv = _edsr1_op(dqd) if tgt == 1 else _edsr2_op(dqd)
                H_idl = _eps2_op(dqd)  if tgt == 1 else _eps1_op(dqd)

                # phi_drive encodes the accumulated virtual Z phase — shifts
                # the rotating frame so the drive appears phase-corrected.
                def _edsr(t, args=None, _t0=t0, _t1=t1, _w=Esigma_GHz, _phi=phi_drive):
                    return np.cos(_w * t + _phi) if _t0 <= t < _t1 else 0.0

                def _idle_x(t, args=None, _t0=t0, _t1=t1, _e=eps_idle):
                    return _e if _t0 <= t < _t1 else 0.0

                terms.append([H_drv, _edsr])
                terms.append([H_idl, _idle_x])

        return [H_drift] + terms, tlist


# ---------------------------------------------------------------------------
# Pulse-schedule visualisation
# ---------------------------------------------------------------------------

def _g_sigma_at_eps(dqd: DQDsystem, eps: float) -> float:
    """
    Effective qubit-photon coupling g_σ [MHz] at DC detuning eps [μeV].

    Computed from the off-diagonal matrix element
        g_σ(ε) = g_c · |⟨g(ε)| τ_z |e(ε)⟩|
    where |g(ε)⟩ and |e(ε)⟩ are the ground and first-excited eigenstates of
        H(ε) = t_c·τ_x + (B_z/2)·σ_z + (b_x/2)·σ_x·τ_z + (ε/2)·τ_z
    """
    H = (dqd.tc * dqd.tx + dqd.Bz / 2 * dqd.sz
         + dqd.bx / 2 * dqd.sx * dqd.tz + eps / 2 * dqd.tz)
    _, evecs = H.eigenstates(sort='low')
    g_vec = evecs[0].full().flatten()
    e_vec = evecs[1].full().flatten()
    return float(abs(g_vec.conj() @ dqd.tz.full() @ e_vec)) * dqd.gc


def plot_pulse_schedule(
    dqd: DQDsystem,
    seq: 'DQDSequenceCompiler',
    n_pts: int = 2000,
    title: str | None = None,
) -> 'matplotlib.figure.Figure':
    """
    Visualise the pulse schedule produced by a DQDSequenceCompiler.

    Outputs
    -----------
    ε1(t)  : DQD1 DC detuning [μeV]  (ε_idle during idle DQD's Rx window)
    ε2(t)  : DQD2 DC detuning [μeV]
    Vac0   : Vac₀ drive amplitude [μeV]; actual carrier is cos(Eσ/ħ·t+φ)
    g_σ(ε) : effective qubit-photon coupling; suppressed by idle detuning

    iSWAP windows are highlighted in green.

    Parameters
    ----------
    dqd   : DQDsystem
    seq   : DQDSequenceCompiler after all gates have been added
    n_pts : number of equally-spaced sample points (default 2000)
    title : optional figure suptitle

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    if not seq._steps:
        raise ValueError("Sequence is empty — add gates before plotting")

    # ── Time axis ──────────────────────────────────────────────────────────
    offsets = [0.0]
    for s in seq._steps:
        offsets.append(offsets[-1] + s['duration'])
    T    = offsets[-1]
    t    = np.linspace(0.0, T, n_pts)
    t_us = t * 1e-3          # ns → μs for x-axis

    # ── Reconstruct ε(t) and EDSR envelope(t) from step list ──────────────
    eps_idle = seq._eps_idle
    eps1 = np.zeros(n_pts)   # DQD1 DC detuning [μeV]
    eps2 = np.zeros(n_pts)   # DQD2 DC detuning [μeV]
    env1 = np.zeros(n_pts)   # DQD1 EDSR envelope [μeV]
    env2 = np.zeros(n_pts)   # DQD2 EDSR envelope [μeV]

    for i, step in enumerate(seq._steps):
        t0s, t1s = offsets[i], offsets[i + 1]
        m = (t >= t0s) & (t < t1s)
        if step['type'] == 'xrot':
            if step['target'] == 1:
                env1[m] = dqd.Vac0
                eps2[m] = eps_idle   # idle DQD gets ε_idle
            else:
                env2[m] = dqd.Vac0
                eps1[m] = eps_idle
        # iswap: ε = 0, no EDSR (default zeros)

    # ── g_σ(ε) — only {0, ε_idle} appear for DQDSequenceCompiler ─────────
    g_at_0    = _g_sigma_at_eps(dqd, 0.0)
    g_at_idle = _g_sigma_at_eps(dqd, eps_idle) if eps_idle != 0.0 else g_at_0

    g1 = np.where(np.abs(eps1) < 1e-6, g_at_0, g_at_idle)
    g2 = np.where(np.abs(eps2) < 1e-6, g_at_0, g_at_idle)

    # ── Build figure ───────────────────────────────────────────────────────
    fig, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True,
                             gridspec_kw={'hspace': 0.08})
    ax_e1, ax_e2, ax_ac, ax_g = axes

    # iSWAP window shading
    _iswap_seen = False
    for i, step in enumerate(seq._steps):
        t0_us = offsets[i] * 1e-3
        t1_us = offsets[i + 1] * 1e-3
        if step['type'] == 'iswap':
            for ax in axes:
                ax.axvspan(t0_us, t1_us, color='green', alpha=0.10, lw=0)
            _iswap_seen = True

    # Gate-boundary annotations on ax_e1
    for i, step in enumerate(seq._steps):
        if step['type'] == 'xrot':
            cx = (offsets[i] + offsets[i + 1]) / 2 * 1e-3
            ax_e1.text(cx, 1.02, f"Rx({np.degrees(step['angle']):.0f}°)\nQ{step['target']}",
                       ha='center', va='bottom', fontsize=7,
                       transform=ax_e1.get_xaxis_transform())

    # ① ε₁
    ax_e1.plot(t_us, eps1, color='C0', lw=1.5)
    ax_e1.axhline(0, color='k', lw=0.5, ls='--', alpha=0.4)
    ax_e1.set_ylabel('ε₁ (μeV)', fontsize=9)

    # ② ε₂
    ax_e2.plot(t_us, eps2, color='C1', lw=1.5)
    ax_e2.axhline(0, color='k', lw=0.5, ls='--', alpha=0.4)
    ax_e2.set_ylabel('ε₂ (μeV)', fontsize=9)

    # ③ EDSR envelopes
    ax_ac.fill_between(t_us, env1, alpha=0.65, color='C0', label='DQD1')
    ax_ac.fill_between(t_us, env2, alpha=0.65, color='C1', label='DQD2')
    ax_ac.set_ylabel(f'EDSR envelope\nVac₀ = {dqd.Vac0} μeV', fontsize=9)
    ax_ac.legend(loc='upper right', fontsize=8)

    # ④ g_σ(ε)
    ax_g.plot(t_us, g1, color='C0', lw=1.5, label='DQD1')
    ax_g.plot(t_us, g2, color='C1', lw=1.5, label='DQD2')
    ax_g.axhline(dqd.g_sigma, color='gray', lw=0.8, ls='--',
                 label=f'g_σ(ε=0) = {dqd.g_sigma:.3f} MHz')
    ax_g.set_ylabel('g_σ(ε) (MHz)', fontsize=9)
    ax_g.set_xlabel('time (μs)')
    ax_g.legend(loc='upper right', fontsize=8)

    # Combined legend: iSWAP indicator
    handles = []
    if _iswap_seen:
        handles.append(mpatches.Patch(fc='green', alpha=0.25, label='iSWAP window'))
    if handles:
        ax_e1.legend(handles=handles, loc='upper left', fontsize=8)

    fig.suptitle(title or 'DQD Pulse Schedule', fontsize=11)
    plt.tight_layout()
    return fig
