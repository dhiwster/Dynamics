# Spin-Spin Coupling Derivation

$$
\text{NCAlgebra available: True}
$$


## Single-site eigenvector rotation for the charge-spin blocks

Using the same style of 2x2 eigenvector rotation as in the IO theory notebook.

$$
h_{\mathrm{block},-} = \begin{pmatrix} (B_z - 2 t_c)/2 & b_x/2 \\ b_x/2 & (-B_z + 2 t_c)/2 \end{pmatrix}
$$

$$
h_{\mathrm{block},+} = \begin{pmatrix} (-B_z - 2 t_c)/2 & b_x/2 \\ b_x/2 & (B_z + 2 t_c)/2 \end{pmatrix}
$$

$$
\theta_{-} = \arctan\left(b_x/(2 t_c - B_z)\right)
$$

$$
\theta_{+} = \arctan\left(b_x/(2 t_c + B_z)\right)
$$

$$
\bar{\theta} = (\theta_{+} + \theta_{-})/2
$$

$$
E_{-} = \sqrt{b_x^2 + (B_z - 2 t_c)^2}/2
$$

$$
E_{+} = \sqrt{b_x^2 + (B_z + 2 t_c)^2}/2
$$

$$
U_{-} = \begin{pmatrix} \cos\left(\theta_{-}/2\right) & -\sin\left(\theta_{-}/2\right) \\ \sin\left(\theta_{-}/2\right) & \cos\left(\theta_{-}/2\right) \end{pmatrix}
$$

$$
U_{+} = \begin{pmatrix} \cos\left(\theta_{+}/2\right) & -\sin\left(\theta_{+}/2\right) \\ \sin\left(\theta_{+}/2\right) & \cos\left(\theta_{+}/2\right) \end{pmatrix}
$$

$$
U_{-} h_{-} U_{-}^{\dagger} = \begin{pmatrix} -\sqrt{b_x^2 + (B_z - 2 t_c)^2}/2 & 0 \\ 0 & \sqrt{b_x^2 + (B_z - 2 t_c)^2}/2 \end{pmatrix}
$$

$$
U_{+} h_{+} U_{+}^{\dagger} = \begin{pmatrix} -\sqrt{b_x^2 + (B_z + 2 t_c)^2}/2 & 0 \\ 0 & \sqrt{b_x^2 + (B_z + 2 t_c)^2}/2 \end{pmatrix}
$$

$$
U_{-} \sigma_x U_{-}^{\dagger} = \begin{pmatrix} -\sin\left(\theta_{-}\right) & \cos\left(\theta_{-}\right) \\ \cos\left(\theta_{-}\right) & \sin\left(\theta_{-}\right) \end{pmatrix}
$$

$$
U_{+} \sigma_x U_{+}^{\dagger} = \begin{pmatrix} -\sin\left(\theta_{+}\right) & \cos\left(\theta_{+}\right) \\ \cos\left(\theta_{+}\right) & \sin\left(\theta_{+}\right) \end{pmatrix}
$$

The helper eigenRotvector2x2 is kept as the implementation route; the matrices above are the compact angle-form equivalents.

At the operator level, this first basis change is what turns the bare dipole coupling into the paper structure -g_{\tau} tau_x + g_{\sigma} tau_z \sigma_x.

Paper mapping: g_{\tau} = gc \cos\left(\bar{\theta}\right), g_{\sigma} = gc \sin\left(\bar{\theta}\right).

This is the basis-change step that produces the later tau_x and tau_z \sigma_x structure.


## Full symbolic model before lower-charge projection

$$
H_{0,\mathrm{full}} = (E_{\sigma,1} \sigma_z^{(1)})/2 + (E_{\sigma,2} \sigma_z^{(2)})/2 + (E_{\tau,1} \tau_z^{(1)})/2 + (E_{\tau,2} \tau_z^{(2)})/2 + \omega_c a^{\dagger}\,a
$$

$$
H_{I,\mathrm{full}} = (-(g_{\tau,1} \tau_x^{(1)}) - g_{\tau,2} \tau_x^{(2)} + g_{\sigma,1} \tau_z^{(1)}\,\\sigma_x^{(1)} + g_{\sigma,2} \tau_z^{(2)}\,\\sigma_x^{(2)})\,(a + a^{\dagger})
$$


$$
Lower-charge projection: tau_z^{1} = tau_z^{2} = -1, tau_x^{i} -> 0
$$

$$
H_{0,\mathrm{proj}} = -1/2 E_{\tau,1} - E_{\tau,2}/2 + (\omega_{q1} \sigma_z^{(1)})/2 + (\omega_{q2} \sigma_z^{(2)})/2 + \omega_c a^{\dagger}\,a
$$

$$
H_{I,\mathrm{proj}} = -(g_{\sigma,1} (\sigma_-^{(1)} + \sigma_+^{(1)})\,(a + a^{\dagger})) - g_{\sigma,2} (\sigma_-^{(2)} + \sigma_+^{(2)})\,(a + a^{\dagger})
$$


## Rotating-frame transformation and rotating-wave approximation

Following the reference notebook pattern: H_rot = U H U^{\dagger} + i (dU/dt) U^{\dagger}.

$$
U_rot(t) = \exp\left(i t (\omega_c a^{\dagger}\,a + (\omega_{q1}/2) \sigma_z^{(1)} + (\omega_{q2}/2) \sigma_z^{(2)})\right)
$$

Detuning convention: \delta_i = \omega_{qi} - \omega_c.

$$
H_{I,\mathrm{rot}}(t) = g_1 (a\,\sigma_-^{(1)}/E^(I (\omega_c + \omega_{q1}) t) + e^{I \delta_1 t} a\,\sigma_+^{(1)} + a^{\dagger}\,\sigma_-^{(1)}/e^{I \delta_1 t} + E^(I (\omega_c + \omega_{q1}) t) a^{\dagger}\,\sigma_+^{(1)}) + g_2 (a\,\sigma_-^{(2)}/E^(I (\omega_c + \omega_{q2}) t) + e^{I \delta_2 t} a\,\sigma_+^{(2)} + a^{\dagger}\,\sigma_-^{(2)}/e^{I \delta_2 t} + E^(I (\omega_c + \omega_{q2}) t) a^{\dagger}\,\sigma_+^{(2)})
$$

RWA drops the fast terms exp[+- i (\omega_c + \omega_{qi}) t] a sm_i and a^{\dagger} sp_i.

$$
H_{I,\mathrm{RWA}}(t) = g_1 (e^{I \delta_1 t} a\,\sigma_+^{(1)} + a^{\dagger}\,\sigma_-^{(1)}/e^{I \delta_1 t}) + g_2 (e^{I \delta_2 t} a\,\sigma_+^{(2)} + a^{\dagger}\,\sigma_-^{(2)}/e^{I \delta_2 t})
$$

For the static dispersive SW step we strip the slow phases and use H_{I,\mathrm{RWA}} = g_1 a\,\sigma_+^{(1)} + g_2 a\,\sigma_+^{(2)} + g_1 a^{\dagger}\,\sigma_-^{(1)} + g_2 a^{\dagger}\,\sigma_-^{(2)}


## Symbolic Schrieffer-Wolff generator for the static RWA Hamiltonian

Using the standard dispersive SW choice from the supplementary derivation.

$$
H_{I,\mathrm{RWA}} = g_1 a\,\sigma_+^{(1)} + g_2 a\,\sigma_+^{(2)} + g_1 a^{\dagger}\,\sigma_-^{(1)} + g_2 a^{\dagger}\,\sigma_-^{(2)}
$$

$$
S = (g_1 a\,\sigma_+^{(1)})/\delta_1 + (g_2 a\,\sigma_+^{(2)})/\delta_2 - (g_1 a^{\dagger}\,\sigma_-^{(1)})/\delta_1 - (g_2 a^{\dagger}\,\sigma_-^{(2)})/\delta_2
$$

SW check ([S,H0] + H_{I,\mathrm{RWA}} == 0): True


## Second-order empty-cavity dispersive Hamiltonian

Evaluating 1/2 [S, H_{I,\mathrm{RWA}}] with same-site spin algebra and cross-site commutativity.

$$
H_{2,\mathrm{vac}} = g_1^2/(2 \delta_1) + g_2^2/(2 \delta_2) + (g_1^2 \sigma_z^{(1)})/(2 \delta_1) + (g_2^2 \sigma_z^{(2)})/(2 \delta_2) + (g_1 g_2 \sigma_-^{(1)}\,\sigma_+^{(2)})/(2 \delta_1) + (g_1 g_2 \sigma_-^{(1)}\,\sigma_+^{(2)})/(2 \delta_2) + (g_1 g_2 \sigma_+^{(1)}\,\sigma_-^{(2)})/(2 \delta_1) + (g_1 g_2 \sigma_+^{(1)}\,\sigma_-^{(2)})/(2 \delta_2)
$$

$$
H_{d,\mathrm{rot}} = (g_1 g_2 \sigma_-^{(1)}\,\sigma_+^{(2)})/(2 \delta_1) + (g_1 g_2 \sigma_-^{(1)}\,\sigma_+^{(2)})/(2 \delta_2) + (g_1 g_2 \sigma_+^{(1)}\,\sigma_-^{(2)})/(2 \delta_1) + (g_1 g_2 \sigma_+^{(1)}\,\sigma_-^{(2)})/(2 \delta_2)
$$

$$
J_{12} = ((\delta_1^{-1} + \delta_2^{-1}) g_1 g_2)/2
$$


## First-order transformed cavity annihilation operator

Using a -> a + [S,a] to first order, then projecting to the empty cavity.

$$
a_{\mathrm{eff},\mathrm{vac}} = (g_1 \sigma_-^{(1)})/\delta_1 + (g_2 \sigma_-^{(2)})/\delta_2
$$

$$
\text{With paper convention } \mathcal{D}[c]\,\rho = 2 c \rho c^{\dagger} - c^{\dagger} c\,\rho - \rho\, c^{\dagger} c
$$

the Purcell term is (kappa/2) D[(g_1 \sigma_-^{(1)})/\delta_1 + (g_2 \sigma_-^{(2)})/\delta_2].


## Identical-DQD specialization matching the paper

$$
H_{d,identical} = (g_{\sigma}^2 \sigma_-^{(1)}\,\sigma_+^{(2)})/\Delta_1 + (g_{\sigma}^2 \sigma_+^{(1)}\,\sigma_-^{(2)})/\Delta_1
$$

$$
J_{identical} = g_{\sigma}^2/\Delta_1
$$

$$
a_{\mathrm{eff},\mathrm{vac},identical} = (g_{\sigma} \sigma_-^{(1)})/\Delta_1 + (g_{\sigma} \sigma_-^{(2)})/\Delta_1
$$

$$
Paper Eq. (2): Hd = (g_{\sigma}^2/\Delta_1) (\sigma_+^{(1)} \, \sigma_-^{(2)} + \sigma_-^{(1)} \, \sigma_+^{(2)}).
$$
