from __future__ import annotations

import numpy as np

from .constants import HBAR_UEV_NS
from .system import DQDsystem

hbar_ns = HBAR_UEV_NS
_MHz2GHz = 1e-3


def _g_sigma_at_eps(dqd: DQDsystem, eps: float) -> float:
    H = (dqd.tc * dqd.tx + dqd.Bz / 2 * dqd.sz
         + dqd.bx / 2 * dqd.sx * dqd.tz + eps / 2 * dqd.tz)
    _, evecs = H.eigenstates(sort='low')
    g_vec = evecs[0].full().flatten()
    e_vec = evecs[1].full().flatten()
    return float(abs(g_vec.conj() @ dqd.tz.full() @ e_vec)) * dqd.gc


def print_summary(dqd: DQDsystem) -> None:
    """Print key derived parameters and gate times for the given DQDsystem."""
    phi_bar     = dqd.phi_bar
    Esigma_GHz  = dqd.Esigma / hbar_ns
    g_sigma_GHz = dqd.g_sigma * _MHz2GHz
    d_sigma_GHz = dqd.d_sigma * _MHz2GHz
    tg_iSWAP    = dqd.iSWAP_gate_time() * 1e3
    t_x180      = np.pi * hbar_ns / (dqd.Vac0 * np.sin(phi_bar))
    r_s, _      = dqd.dispersive_ratios()

    print(f"DQD:        tc={dqd.tc} μeV,  Bz={dqd.Bz} μeV (fixed),  bx={dqd.bx} μeV (fixed)")
    print(f"Cavity:     wc={dqd.wc/2/np.pi:.0f} MHz,  gc={dqd.gc/2/np.pi:.1f} MHz")
    print(f"Qubit:      Eσ={dqd.Esigma:.4f} μeV  ({Esigma_GHz:.3f} GHz)")
    print(f"φ̄:          {np.degrees(phi_bar):.3f}°,  sin(φ̄)={np.sin(phi_bar):.6f}")
    print(f"Coupling:   g_σ={g_sigma_GHz*1e3:.4f} MHz,  Δ_σ={d_sigma_GHz*1e3:.2f} MHz,  g/Δ={r_s:.5f}")
    print(f"tg(iSWAP):  {tg_iSWAP:.4f} ns")
    print(f"t(Rx90):    {t_x180/2:.4f} ns,  t(Rx180): {t_x180:.4f} ns")
    print(f"Vac0:       {dqd.Vac0} μeV  (Rabi={dqd.Vac0/hbar_ns*np.sin(phi_bar)*1e3:.3f} MHz)")
    print(f"ε_idle:     {dqd.epsilon_idle:.1f} μeV")


def plot_pulse_schedule(
    dqd: DQDsystem,
    seq,
    n_pts: int = 2000,
    title: str | None = None,
) -> 'matplotlib.figure.Figure':
    """Visualise the pulse schedule produced by a DQDSequenceCompiler."""
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    if not seq._steps:
        raise ValueError("Sequence is empty — add gates before plotting")

    offsets = [0.0]
    for s in seq._steps:
        offsets.append(offsets[-1] + s['duration'])
    T    = offsets[-1]
    t    = np.linspace(0.0, T, n_pts)
    t_us = t * 1e-3

    eps_idle = seq._eps_idle
    eps1 = np.zeros(n_pts)
    eps2 = np.zeros(n_pts)
    env1 = np.zeros(n_pts)
    env2 = np.zeros(n_pts)

    for i, step in enumerate(seq._steps):
        t0s, t1s = offsets[i], offsets[i + 1]
        m = (t >= t0s) & (t < t1s)
        if step['type'] == 'xrot':
            if step['target'] == 1:
                env1[m] = dqd.Vac0
                eps2[m] = eps_idle
            else:
                env2[m] = dqd.Vac0
                eps1[m] = eps_idle

    g_at_0    = _g_sigma_at_eps(dqd, 0.0)
    g_at_idle = _g_sigma_at_eps(dqd, eps_idle) if eps_idle != 0.0 else g_at_0

    g1 = np.where(np.abs(eps1) < 1e-6, g_at_0, g_at_idle)
    g2 = np.where(np.abs(eps2) < 1e-6, g_at_0, g_at_idle)

    fig, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True,
                             gridspec_kw={'hspace': 0.08})
    ax_e1, ax_e2, ax_ac, ax_g = axes

    _iswap_seen = False
    for i, step in enumerate(seq._steps):
        t0_us = offsets[i] * 1e-3
        t1_us = offsets[i + 1] * 1e-3
        if step['type'] == 'iswap':
            for ax in axes:
                ax.axvspan(t0_us, t1_us, color='green', alpha=0.10, lw=0)
            _iswap_seen = True

    for i, step in enumerate(seq._steps):
        if step['type'] == 'xrot':
            cx = (offsets[i] + offsets[i + 1]) / 2 * 1e-3
            ax_e1.text(cx, 1.02, f"Rx({np.degrees(step['angle']):.0f}°)\nQ{step['target']}",
                       ha='center', va='bottom', fontsize=7,
                       transform=ax_e1.get_xaxis_transform())

    ax_e1.plot(t_us, eps1, color='C0', lw=1.5)
    ax_e1.axhline(0, color='k', lw=0.5, ls='--', alpha=0.4)
    ax_e1.set_ylabel('ε₁ (μeV)', fontsize=9)

    ax_e2.plot(t_us, eps2, color='C1', lw=1.5)
    ax_e2.axhline(0, color='k', lw=0.5, ls='--', alpha=0.4)
    ax_e2.set_ylabel('ε₂ (μeV)', fontsize=9)

    ax_ac.fill_between(t_us, env1, alpha=0.65, color='C0', label='DQD1')
    ax_ac.fill_between(t_us, env2, alpha=0.65, color='C1', label='DQD2')
    ax_ac.set_ylabel(f'EDSR envelope\nVac₀ = {dqd.Vac0} μeV', fontsize=9)
    ax_ac.legend(loc='upper right', fontsize=8)

    ax_g.plot(t_us, g1, color='C0', lw=1.5, label='DQD1')
    ax_g.plot(t_us, g2, color='C1', lw=1.5, label='DQD2')
    ax_g.axhline(dqd.g_sigma, color='gray', lw=0.8, ls='--',
                 label=f'g_σ(ε=0) = {dqd.g_sigma:.3f} MHz')
    ax_g.set_ylabel('g_σ(ε) (MHz)', fontsize=9)
    ax_g.set_xlabel('time (μs)')
    ax_g.legend(loc='upper right', fontsize=8)

    handles = []
    if _iswap_seen:
        handles.append(mpatches.Patch(fc='green', alpha=0.25, label='iSWAP window'))
    if handles:
        ax_e1.legend(handles=handles, loc='upper left', fontsize=8)

    fig.suptitle(title or 'DQD Pulse Schedule', fontsize=11)
    plt.tight_layout()
    return fig
