import numpy as np
import qutip as qt
from qutip import (tensor, sigmax, sigmay, sigmaz, sigmam, sigmap,
                   qeye, basis, destroy)

Pi2  = 2 * np.pi
hbar = 6.582119569e-4  # ℏ in μeV·ns (so energies in μeV → frequencies in MHz)


class DQDsystem:
    """
    Two double quantum dots (DQDs) coupled via a microwave resonator.

    **Basis** — spin ⊗ charge: |↑,R⟩  |↑,L⟩  |↓,R⟩  |↓,L⟩

    - σ_z : ↑ = +1 (index 0),  ↓ = −1 (index 1)
    - τ_z : R = +1 (index 0),  L = −1 (index 1)

    **Single-DQD Hamiltonian** (symmetric, ε = 0):

        H_DQD = t_c·τ_x + (B_z/2)·σ_z + (b_x/2)·σ_x·τ_z

    **Full system Hamiltonian** (photon ⊗ DQD₁ ⊗ DQD₂, units: MHz):

        H = H_DQD1/ħ + H_DQD2/ħ + ω_c·a†a
          + g_c·τ_z1·(a + a†) + g_c·τ_z2·(a + a†)

    Reference: https://arxiv.org/abs/1902.07649
    """

    def __init__(self, tc, bx, Bz, Vac0, wc=Pi2 * 5e3, gc=Pi2 * 50, photon_max=10,
                 epsilon_idle=0, tc_idle=None):
        """
        Parameters
        ----------
        tc        : tunneling coupling (μeV)
        bx        : transverse magnetic field gradient (μeV)
        Bz        : Zeeman splitting (μeV), positive → ↑ has higher energy
        wc        : cavity angular frequency (MHz)
        gc        : bare spin-photon coupling (MHz)
        photon_max: Fock-space truncation for the resonator
        epsilon   : detuning / dot asymmetry (μeV), 0 = symmetric DQD
        """
        self.tc         = tc
        self.bx         = bx
        self.Bz         = Bz
        self.Vac0       = Vac0
        self.wc         = wc
        self.gc         = gc
        self.photon_max = photon_max
        self.epsilon_idle = epsilon_idle
        self.epsilon      = 0.0          # sweet-spot default; only ε_idle is pulsed
        self.tc_idle      = tc if tc_idle is None else tc_idle

        self._build_single_dqd_operators()
        self._build_two_dqd_operators()
        self._compute_eigenbasis()
        self._build_hamiltonians()

    # ------------------------------------------------------------------
    # Operator construction
    # ------------------------------------------------------------------

    def _build_single_dqd_operators(self):
        """Spin and charge Pauli operators for a single DQD (4×4)."""
        self.sx = tensor(sigmax(), qeye(2))
        self.sy = tensor(sigmay(), qeye(2))
        self.sz = tensor(sigmaz(), qeye(2))
        self.sm = tensor(sigmam(), qeye(2))
        self.sp = tensor(sigmap(), qeye(2))

        self.tx = tensor(qeye(2), sigmax())
        self.ty = tensor(qeye(2), sigmay())
        self.tz = tensor(qeye(2), sigmaz())

    def _build_two_dqd_operators(self):
        """Operators embedded in the two-DQD space (16×16) and full system."""
        # DQD1: acts on first 4-dim factor, identity on second
        self.sx1 = tensor(self.sx, qeye([2, 2]))
        self.sy1 = tensor(self.sy, qeye([2, 2]))
        self.sz1 = tensor(self.sz, qeye([2, 2]))
        self.sm1 = tensor(self.sm, qeye([2, 2]))
        self.sp1 = tensor(self.sp, qeye([2, 2]))
        self.tx1 = tensor(self.tx, qeye([2, 2]))
        self.ty1 = tensor(self.ty, qeye([2, 2]))
        self.tz1 = tensor(self.tz, qeye([2, 2]))

        # DQD2: identity on first 4-dim factor, acts on second
        self.sx2 = tensor(qeye([2, 2]), self.sx)
        self.sy2 = tensor(qeye([2, 2]), self.sy)
        self.sz2 = tensor(qeye([2, 2]), self.sz)
        self.sm2 = tensor(qeye([2, 2]), self.sm)
        self.sp2 = tensor(qeye([2, 2]), self.sp)
        self.tx2 = tensor(qeye([2, 2]), self.tx)
        self.ty2 = tensor(qeye([2, 2]), self.ty)
        self.tz2 = tensor(qeye([2, 2]), self.tz)

        # Photon operator in full space: photon ⊗ DQD1 ⊗ DQD2
        self.a = tensor(destroy(self.photon_max), qeye([2, 2, 2, 2]))

    # ------------------------------------------------------------------
    # Eigenbasis and derived parameters
    # ------------------------------------------------------------------

    def _compute_eigenbasis(self):
        """
        Diagonalisation angles and effective coupling parameters.

        Energy eigenvalues of H_DQD (at ε=0):
            ±e2 = ±½ √[(2t_c − B_z)² + b_x²]   (lower doublet)
            ±e3 = ±½ √[(2t_c + B_z)² + b_x²]   (upper doublet)

        Qubit transition:  E_σ = e3 − e2
        Charge transition: E_τ = e3 + e2

        Rotation angles (see paper Appendix):
            φ₊ = arctan(b_x / (2t_c + B_z))
            φ₋ = arctan(b_x / (2t_c − B_z))
            φ̄  = (φ₊ + φ₋) / 2

        Effective couplings in the dispersive eigenbasis:
            g_σ = g_c sin(φ̄)    (qubit coupling)
            g_τ = g_c cos(φ̄)    (charge coupling)
        """
        tc, bx, Bz = self.tc, self.bx, self.Bz

        phi_p = np.arctan2(bx, 2 * tc + Bz)
        if phi_p < 0:
            phi_p += np.pi

        phi_m = np.arctan2(bx, 2 * tc - Bz)
        if phi_m < 0:
            phi_m += np.pi

        self.phi_p   = phi_p
        self.phi_m   = phi_m
        self.phi_bar = 0.5 * (phi_p + phi_m)

        self.e2 = 0.5 * np.sqrt((2 * tc - Bz) ** 2 + bx ** 2)
        self.e3 = 0.5 * np.sqrt((2 * tc + Bz) ** 2 + bx ** 2)

        self.Esigma = self.e3 - self.e2  # μeV
        self.Etau   = self.e3 + self.e2  # μeV

        # Detunings from the cavity (MHz)
        self.d_sigma = self.Esigma / hbar - self.wc
        self.d_tau   = self.Etau   / hbar - self.wc

        # Effective dispersive coupling strengths (MHz)
        self.g_sigma = self.gc * np.sin(self.phi_bar)
        self.g_tau   = self.gc * np.cos(self.phi_bar)

    # ------------------------------------------------------------------
    # Hamiltonian
    # ------------------------------------------------------------------

    def _build_hamiltonians(self):
        r"""
        Construct all Hamiltonian terms in the lab frame.
        All energies in MHz (i.e. DQD terms divided by hbar).

        H = H_DQD1 + H_DQD2 + H_photon + H_int1 + H_int2

        The DQD terms are split so that tc and epsilon can be made
        time-dependent via the operator properties H_tc1_op etc.:
            H_DQDi = H_DQDi_static + tc*H_tc{i}_op + epsilon*H_eps{i}_op
        """
        tc, bx, Bz, eps = self.tc, self.bx, self.Bz, self.epsilon
        wc, gc = self.wc, self.gc
        N = self.photon_max
        a = self.a

        # Static parts
        H_static1 = tensor(qeye(N),
                               tc * self.tx1 + Bz / 2 * self.sz1
                                + bx / 2 * self.sx1 * self.tz1) / hbar
        H_static2 = tensor(qeye(N),
                               tc * self.tx2 + Bz / 2 * self.sz2
                               + bx / 2 * self.sx2 * self.tz2) / hbar
        H_photon = wc * a.dag() * a
        H_int1   = gc * tensor(qeye(N), self.tz1) * (a + a.dag())
        H_int2   = gc * tensor(qeye(N), self.tz2) * (a + a.dag())
        self.H_static = H_static1 + H_static2 + H_photon + H_int1 + H_int2
        self.H_int1    = H_int1
        self.H_int2    = H_int2

        # Dynamic control operators. H_static already includes tc and epsilon=0,
        # so time-dependent tc controls should be supplied as offsets from tc.
        self.H_tc1_op = tensor(qeye(N), self.tx1) / hbar
        self.H_tc2_op = tensor(qeye(N), self.tx2) / hbar
        self.H_eps1_op = tensor(qeye(N), self.tz1) / (2 * hbar)
        self.H_eps2_op = tensor(qeye(N), self.tz2) / (2 * hbar)

        # Detuned parts (dynamic)
        self.H_detuned1  = tensor(qeye(N), (self.tc_idle - self.tc) * self.tx1 + self.epsilon_idle / 2 * self.tz1) / hbar
        self.H_detuned2  = tensor(qeye(N), (self.tc_idle - self.tc) * self.tx2 + self.epsilon_idle / 2 * self.tz2) / hbar
        
        # EDSR terms (dynamic)
        self.H_edsr1_amplitude = tensor(qeye(N), self.Vac0 * self.tz1)
        self.H_edsr2_amplitude = tensor(qeye(N), self.Vac0 * self.tz2)


    # ------------------------------------------------------------------
    # Diagonalisation
    # ------------------------------------------------------------------

    def diagonalization_unitary(self):
        """
        Diagonalisation unitaries for a single DQD.

        Returns
        -------
        U_flop : rotation to spin-charge hybrid eigenbasis
                 U_flop = R_sy_tx(1/2(φ₊−φ₋)) · R_sx_ty(φbar)
        U_ob   : rotation to charge eigenbasis  R_ty(−π/2)

        Full single-DQD diagonalisation: U = U_flop · U_ob
        """
        Rty    = lambda theta: tensor(qeye(2), (-0.5j * sigmay() * theta).expm())
        Rsxty  = (-0.5j * self.sx * self.ty * self.phi_bar).expm()
        Rsytx  = (-0.5j * self.sy * self.tx * 0.5 * (self.phi_p - self.phi_m)).expm()
        U_flop = Rsytx * Rsxty
        U_ob   = Rty(-np.pi / 2)
        return U_flop, U_ob

    def diagonalized_H_DQD(self):
        """Return the single-DQD Hamiltonian in the diagonal eigenbasis."""
        H_single = self.tc * self.tx + self.Bz / 2 * self.sz + self.bx / 2 * self.sx * self.tz
        U_flop, U_ob = self.diagonalization_unitary()
        U = U_flop * U_ob
        return (U * H_single * U.dag()).tidyup(atol=1e-10)

    # ------------------------------------------------------------------
    # State preparation
    # ------------------------------------------------------------------

    def initialize_state(self, dqd1_spin='up', dqd2_spin='down'):
        """
        Build the initial state in the lab frame.

        Both DQDs start in the right dot.  The state is constructed in
        the computational basis and then rotated into the lab frame via
        the inverse of the diagonalisation unitary.

        Parameters
        ----------
        dqd1_spin : 'up'  → |↑,R⟩   'down' → |↓,R⟩
        dqd2_spin : 'up'  → |↑,R⟩   'down' → |↓,R⟩

        Returns
        -------
        psi0 : ket in full Hilbert space (photon ⊗ DQD1 ⊗ DQD2)
        """
        s    = {'up': 0, 'down': 1}
        psi1 = basis([2, 2], [s[dqd1_spin], 0])
        psi2 = basis([2, 2], [s[dqd2_spin], 0])
        psi0 = tensor(basis(self.photon_max, 0), psi1, psi2)

        U_flop, U_ob = self.diagonalization_unitary()
        U_full    = tensor(qeye(self.photon_max), U_flop, U_flop)
        U_ob_full = tensor(qeye(self.photon_max), U_ob,   U_ob)
        return U_full * U_ob_full * psi0

    # ------------------------------------------------------------------
    # Dispersive regime diagnostics
    # ------------------------------------------------------------------

    def iSWAP_gate_time(self):
        """
        iSWAP gate time in microseconds () from the dispersive approximation.

        t_g = (π/2) |Δ_σ| / g_σ²    if Δ_σ < 0
            = (3π/2)|Δ_σ| / g_σ²   if Δ_σ > 0
        """
        sign = 1 if self.d_sigma < 0 else 3
        return sign * np.pi / 2 * abs(self.d_sigma) / self.g_sigma ** 2

    def dispersive_ratios(self):
        """
        Return (r_σ, r_τ) = (|g_σ/Δ_σ|, |g_τ/Δ_τ|).
        Dispersive regime requires |g_τ/Δ_τ| ≪ |g_σ/Δ_σ| ≪ 1)
        """
        r_sigma = abs(self.g_sigma / self.d_sigma)
        r_tau   = abs(self.g_tau   / self.d_tau)
        return r_sigma, r_tau
