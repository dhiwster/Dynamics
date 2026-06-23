import numpy as np
from qutip import basis, destroy, qeye, sigmay, tensor

from hamiltonians.models import (
    RAD_PER_US_PER_MICROELECTRONVOLT,
    TWO_PI,
    SingleDQDOperatorsQutip,
    single_dqd_qutip_hamiltonian,
    single_dqd_qutip_operators,
    single_dqd_spectrum,
)

Pi2 = TWO_PI
RAD_PER_US_PER_MICROEV = RAD_PER_US_PER_MICROELECTRONVOLT


class DQDsystem:
    """
    Two double quantum dots coupled through a microwave resonator.

    Single-DQD basis ordering: |up,R>, |up,L>, |down,R>, |down,L>.
    """

    def __init__(
        self,
        tc,
        bx,
        Bz,
        Vac0,
        wc=Pi2 * 5e3,
        gc=Pi2 * 50,
        photon_max=10,
        epsilon_idle=0,
        tc_idle=None,
    ):
        self.tc = tc
        self.bx = bx
        self.Bz = Bz
        self.Vac0 = Vac0
        self.wc = wc
        self.gc = gc
        self.photon_max = photon_max
        self.epsilon_idle = epsilon_idle
        self.epsilon = 0.0
        self.tc_idle = tc if tc_idle is None else tc_idle

        self._build_single_dqd_operators()
        self._build_two_dqd_operators()
        self._compute_eigenbasis()
        self._build_hamiltonians()

    # ------------------------------------------------------------------
    # Operator construction
    # ------------------------------------------------------------------

    def _build_single_dqd_operators(self):
        """Spin and charge Pauli operators for one DQD."""
        ops = single_dqd_qutip_operators()
        self.sx = ops.sx
        self.sy = ops.sy
        self.sz = ops.sz
        self.sm = ops.sm
        self.sp = ops.sp
        self.tx = ops.tx
        self.ty = ops.ty
        self.tz = ops.tz

    def _build_two_dqd_operators(self):
        """Embed single-DQD operators in the two-DQD Hilbert space."""
        self.sx1 = tensor(self.sx, qeye([2, 2]))
        self.sy1 = tensor(self.sy, qeye([2, 2]))
        self.sz1 = tensor(self.sz, qeye([2, 2]))
        self.sm1 = tensor(self.sm, qeye([2, 2]))
        self.sp1 = tensor(self.sp, qeye([2, 2]))
        self.tx1 = tensor(self.tx, qeye([2, 2]))
        self.ty1 = tensor(self.ty, qeye([2, 2]))
        self.tz1 = tensor(self.tz, qeye([2, 2]))

        self.sx2 = tensor(qeye([2, 2]), self.sx)
        self.sy2 = tensor(qeye([2, 2]), self.sy)
        self.sz2 = tensor(qeye([2, 2]), self.sz)
        self.sm2 = tensor(qeye([2, 2]), self.sm)
        self.sp2 = tensor(qeye([2, 2]), self.sp)
        self.tx2 = tensor(qeye([2, 2]), self.tx)
        self.ty2 = tensor(qeye([2, 2]), self.ty)
        self.tz2 = tensor(qeye([2, 2]), self.tz)

        self.ops1 = SingleDQDOperatorsQutip(
            sx=self.sx1,
            sy=self.sy1,
            sz=self.sz1,
            sm=self.sm1,
            sp=self.sp1,
            tx=self.tx1,
            ty=self.ty1,
            tz=self.tz1,
        )
        self.ops2 = SingleDQDOperatorsQutip(
            sx=self.sx2,
            sy=self.sy2,
            sz=self.sz2,
            sm=self.sm2,
            sp=self.sp2,
            tx=self.tx2,
            ty=self.ty2,
            tz=self.tz2,
        )

        self.a = tensor(destroy(self.photon_max), qeye([2, 2, 2, 2]))

    # ------------------------------------------------------------------
    # Eigenbasis and derived parameters
    # ------------------------------------------------------------------

    def _compute_eigenbasis(self):
        """Cache the single-DQD hybridization data used throughout the model."""
        spectrum = single_dqd_spectrum(self)
        self.phi_p = spectrum.phi_p
        self.phi_m = spectrum.phi_m
        self.phi_bar = spectrum.phi_bar
        self.e2 = spectrum.e2
        self.e3 = spectrum.e3
        self.Esigma = spectrum.Esigma
        self.Etau = spectrum.Etau

        self.d_sigma = self.Esigma * RAD_PER_US_PER_MICROEV - self.wc
        self.d_tau = self.Etau * RAD_PER_US_PER_MICROEV - self.wc
        self.g_sigma = self.gc * np.sin(self.phi_bar)
        self.g_tau = self.gc * np.cos(self.phi_bar)

    # ------------------------------------------------------------------
    # Hamiltonian
    # ------------------------------------------------------------------

    def _build_hamiltonians(self):
        """Construct the full photon + two-DQD Hamiltonian in MHz units."""
        wc, gc = self.wc, self.gc
        N = self.photon_max
        a = self.a

        H_static1 = (
            tensor(qeye(N), single_dqd_qutip_hamiltonian(self, self.epsilon, self.ops1))
            * RAD_PER_US_PER_MICROEV
        )
        H_static2 = (
            tensor(qeye(N), single_dqd_qutip_hamiltonian(self, self.epsilon, self.ops2))
            * RAD_PER_US_PER_MICROEV
        )
        H_photon = wc * a.dag() * a
        H_int1 = gc * tensor(qeye(N), self.tz1) * (a + a.dag())
        H_int2 = gc * tensor(qeye(N), self.tz2) * (a + a.dag())
        self.H_static = H_static1 + H_static2 + H_photon + H_int1 + H_int2
        self.H_int1 = H_int1
        self.H_int2 = H_int2

        self.H_tc1_op = tensor(qeye(N), self.tx1) * RAD_PER_US_PER_MICROEV
        self.H_tc2_op = tensor(qeye(N), self.tx2) * RAD_PER_US_PER_MICROEV
        self.H_eps1_op = tensor(qeye(N), self.tz1) * (0.5 * RAD_PER_US_PER_MICROEV)
        self.H_eps2_op = tensor(qeye(N), self.tz2) * (0.5 * RAD_PER_US_PER_MICROEV)

        self.H_detuned1 = tensor(
            qeye(N),
            (self.tc_idle - self.tc) * self.tx1 + self.epsilon_idle / 2 * self.tz1,
        ) * RAD_PER_US_PER_MICROEV
        self.H_detuned2 = tensor(
            qeye(N),
            (self.tc_idle - self.tc) * self.tx2 + self.epsilon_idle / 2 * self.tz2,
        ) * RAD_PER_US_PER_MICROEV

        self.H_edsr1_amplitude = tensor(qeye(N), self.Vac0 * self.tz1)
        self.H_edsr2_amplitude = tensor(qeye(N), self.Vac0 * self.tz2)

    # ------------------------------------------------------------------
    # Diagonalisation
    # ------------------------------------------------------------------

    def diagonalization_unitary(self):
        """
        Return the analytic single-DQD diagonalisation unitaries.

        U_flop rotates to the spin-charge hybrid basis and U_ob rotates to the
        charge eigenbasis.
        """

        Rty = lambda theta: tensor(qeye(2), (-0.5j * sigmay() * theta).expm())
        Rsxty = (-0.5j * self.sx * self.ty * self.phi_bar).expm()
        Rsytx = (-0.5j * self.sy * self.tx * 0.5 * (self.phi_p - self.phi_m)).expm()
        U_flop = Rsytx * Rsxty
        U_ob = Rty(-np.pi / 2)
        return U_flop, U_ob

    def diagonalized_H_DQD(self):
        """Return the single-DQD Hamiltonian in the diagonal eigenbasis."""
        H_single = single_dqd_qutip_hamiltonian(self, epsilon=0.0)
        U_flop, U_ob = self.diagonalization_unitary()
        U = U_flop * U_ob
        return (U * H_single * U.dag()).tidyup(atol=1e-10)

    # ------------------------------------------------------------------
    # State preparation
    # ------------------------------------------------------------------

    def initialize_state(self, dqd1_spin="up", dqd2_spin="down"):
        """Build the default lab-frame initial state."""
        s = {"up": 0, "down": 1}
        psi1 = basis([2, 2], [s[dqd1_spin], 0])
        psi2 = basis([2, 2], [s[dqd2_spin], 0])
        psi0 = tensor(basis(self.photon_max, 0), psi1, psi2)

        U_flop, U_ob = self.diagonalization_unitary()
        U_full = tensor(qeye(self.photon_max), U_flop, U_flop)
        U_ob_full = tensor(qeye(self.photon_max), U_ob, U_ob)
        return U_full * U_ob_full * psi0

    # ------------------------------------------------------------------
    # Dispersive regime diagnostics
    # ------------------------------------------------------------------

    def iSWAP_gate_time(self):
        """Return the dispersive iSWAP gate time in microseconds."""
        sign = 1 if self.d_sigma < 0 else 3
        return sign * np.pi / 2 * abs(self.d_sigma) / self.g_sigma**2

    def dispersive_ratios(self):
        """Return (|g_sigma / d_sigma|, |g_tau / d_tau|)."""
        r_sigma = abs(self.g_sigma / self.d_sigma)
        r_tau = abs(self.g_tau / self.d_tau)
        return r_sigma, r_tau
