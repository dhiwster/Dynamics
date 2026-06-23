Source workflow notebook: `C:\Users\wishuang\Codex\DynamicSimulations\python\notebooks\two_qubit_gate_with_resonator.ipynb`

Published workspace copy: `C:\Users\wishuang\Codex\DynamicSimulations\references\published\two_qubit_gate_with_resonator_workflow.md`

Gate-validation notebook: `C:\Users\wishuang\Codex\DynamicSimulations\python\notebooks\cnot_from_iswaps.ipynb`

Reference notebook: `C:\Users\wishuang\Codex\DynamicSimulations\python\notebooks\DQDlabframe.ipynb`

# Two-qubit gate with resonator workflow

This note follows the code paths that are currently implemented.
The notebook is now arranged as one linear workflow rather than separate disconnected blocks.

## Current workflow notebook parameters

- The two DQDs are still symmetric.
- The current parameter block uses:
  - `tc = ghz_to_uev(10.0)`
  - `bx = 3 ueV`
  - `Bz = 35 ueV`
  - `wc = 2 pi * 5.5e3`
  - `gc = 2 pi * 50`
  - `photon_max = 10`
- Preparation detuning is set to `epsilon_1 = epsilon_2 = 80 GHz` and converted into the Hamiltonian detuning coefficient in `ueV`.
- The prep drive is applied only to DQD1.
- `omega_drive = Bz / hbar`.
- `Omega_drive = 2 pi * 100 MHz`.
- The notebook now uses explicit linear detuning ramps both into and out of the two-qubit gate.
- The current ramp heuristic is:
  - charge period from `tc`: about `0.05 ns`
  - ramp cycles: `20`
  - ramp time: about `1.0 ns`
  - ramp samples: `1601`
- The dispersive gate estimate is still:
  - `t_gate ~ 455.55 us`

## Linear notebook flow

### 1. Parameters and model

The notebook first builds one symmetric `DQDsystem` object and reports the derived scales used later:

- prep detuning in `ueV`
- drive frequencies
- charge-period estimate from `tc`
- ramp duration and sample count
- dispersive gate estimate
- dispersive ratios

### 2. Calibrate the finite-detuning preparation pulse

The prep starts from the bare lab-frame state:

- photon vacuum
- DQD1 in `|down, R>`
- DQD2 in `|down, R>`

The prep Hamiltonian contains:

- the full static photon plus two-DQD Hamiltonian
- finite detuning on both DQDs
- a time-dependent DQD1-only drive

The notebook chooses the shortest prep duration that brings the reduced DQD1 spin populations closest to 50/50.
For the current parameter set that crossing occurs at about `2.379 ns`, with `p_up ~= 0.5008` and `p_down ~= 0.4992`.

In the true bare basis just before the ramp, the largest workflow-prep populations are approximately:

- `Rup,Rdown`: `0.4850`
- `Rdown,Rdown`: `0.4799`
- `Ldown,Rdown`: `0.0102`
- `Rup,Ldown`: `0.0078`

### 3. Run the main workflow sequence

The main workflow now runs in one direct sequence:

1. finite-detuning prep
2. adiabatic ramp from `epsilon_prepare` to `0`
3. gate evolution at `epsilon = 0`
4. adiabatic ramp back to finite detuning for readout

The main section also plots:

- the schedule for detuning and drive amplitude
- gate-stage spin, spin-pair, and charge-right observables
- readout-ramp spin, spin-pair, and charge-right observables
- a 16-state bare-basis bar chart after the whole sequence

The readout summary is reported from the post-ramp finite-detuning state, while also keeping the gate-end values at `epsilon = 0` for comparison.

### 4. Bell-state analysis with the same workflow prep

The Bell analysis stays consistent with the same workflow prep model instead of switching to the gate-helper single-qubit pulse used in `cnot_from_iswaps.ipynb`.
It:

- searches near the first pi-like maximum of the same finite-detuning DQD1 prep pulse
- ramps `epsilon` linearly to `0`
- scans the exchange hold time with a coarse grid and then a local refined grid
- scores the result using the Bell overlap inside the projected zero-photon qubit subspace, multiplied by the weight that remains in that subspace
- ramps the selected gate-point state back to finite detuning before the final readout-spin comparison

For the current symmetric workflow parameters with the ramps included, the best Bell-like point found by that scan is approximately:

- prep time: `4.969 ns`
- gate hold: `95.36 us`
- projected Bell overlap: `0.8960`
- projected concurrence: `0.8762`
- qubit-subspace weight: `3.64e-4`
- full-state Bell weight: `3.26e-4`
- gate-point bare-spin Bell overlap after tracing out charge and cavity: `0.0326`
- gate-point bare-spin concurrence after tracing out charge and cavity: `0.8838`
- readout-point bare-spin Bell overlap after the ramp back to finite detuning: `0.0344`
- readout-point bare-spin concurrence after the ramp back to finite detuning: `0.8819`

So once the finite-detuning to zero-detuning transition is treated as an explicit adiabatic ramp, the state that reaches the gate point is no longer concentrated in the effective zero-photon qubit manifold used by the projected Bell analysis.

### 5. Plot all 16 true bare-basis populations

The Bell section plots all 16 zero-photon DQD populations across the selected prep, ramp, gate, and readout-ramp stages in the true bare basis:

- `Rup`
- `Lup`
- `Rdown`
- `Ldown`

for each DQD.

It also plots three dedicated 16-state bar charts:

- right after the single-qubit gate and before the ramp
- after the ramp and before the two-qubit gate
- after the readout ramp back to finite detuning

### 6. Trace out charge and compare the spin Bell state

At the selected Bell-like point, the notebook traces out charge and cavity after the readout ramp and plots the reduced two-spin state in the `uu`, `ud`, `du`, `dd` basis, including both:

- spin populations
- the magnitude of the reduced density matrix compared with the ideal Bell-state target

## Separate gate-based entanglement validation

A separate gate-level validation lives in `cnot_from_iswaps.ipynb`.
That notebook uses the reusable single-qubit and two-qubit gate helpers directly.

Its Bell-state demo is:

- start from `|gg>` in the vacuum sector
- apply a calibrated `Rx(pi)` on Q2
- hold the exchange interaction for a calibrated fraction of the iSWAP time

That route still gives a cleaner effective Bell state because the single-qubit preparation there is already aligned with the gate-level qubit basis and does not include the finite-detuning workflow ramp.

## Current status

- `two_qubit_gate_with_resonator.ipynb` now runs as one linear workflow.
- A published workspace copy exists at `references/published/two_qubit_gate_with_resonator_workflow.md`.
- The notebook uses `80 GHz` prep detuning.
- The notebook includes explicit linear ramps both into and out of the `epsilon = 0` gate point.
- The Bell section includes true bare-basis 16-state plots over the full prep-plus-ramp-plus-selected-gate-plus-readout sequence.
- The ramp strongly changes the state before the gate, so the projected effective-qubit Bell metric becomes much smaller.
- The DQDs are still symmetric.
- The main next extension is still the same one: allow asymmetric DQD1 and DQD2 parameter sets.
