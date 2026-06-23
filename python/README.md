Python layout
=============

## Structure

```
python/
├── helpers/     ← shared support code (importable package)
└── notebooks/   ← experiments, interactive notebooks, and one-off scripts
```

---

## helpers/

The `helpers` package is the shared support layer imported by notebooks and
scripts. It is split into focused modules:

| Module | Contents |
|--------|----------|
| `constants.py` | Physical constants and unit-conversion factors |
| `operators.py` | Single-DQD Pauli operators (numpy and QuTiP variants) |
| `hamiltonian.py` | Single-DQD Hamiltonian, spectrum, eigensystem functions |
| `pulses.py` | Pulse-shape constructors: `square`, `constant`, `piecewise` |
| `system.py` | `DQDsystem` — two DQDs coupled through a microwave resonator |
| `gates.py` | Gate Hamiltonians (`iswap_H`, `xrot_H`, `zrot_H`), operator builders (`drift_H`, `eps1_op`, `eps2_op`, `edsr1_op`, `edsr2_op`), calibration helpers, `build_H` |
| `sequence.py` | `DQDSequenceCompiler` — virtual-Z phase tracking and gate compilation |
| `solver.py` | `run_state`, `run_unitary`, and post-processing projection functions |
| `viz.py` | `plot_pulse_schedule`, `print_summary` |

`helpers/__init__.py` re-exports the full public API. Import everything from
`helpers` directly — do not import private internals or sub-modules:

```python
# Good
from helpers import DQDsystem, DQDSequenceCompiler, run_unitary, drift_H

# Avoid
from helpers.gates import _some_internal
```

QuTiP-dependent names are wrapped in a `try/except ModuleNotFoundError` so
that the pure-numpy path (`constants`, `operators`, `hamiltonian`, `pulses`)
works without QuTiP installed.

---

## notebooks/

Notebooks and scripts are the primary work. Each file is a self-contained
experiment or interactive demonstration. They import from `helpers` and are
not imported by each other (except when a `.ipynb` wraps a companion `.py`
script for interactive widgets).

### Rules

- **New reusable logic goes in `helpers/`** — if the same computation appears
  in two notebooks, extract it to the appropriate `helpers` module.
- **New experiments go in `notebooks/`** — problem-specific code that is not
  shared stays here as a script or notebook.
- **No cross-notebook imports** — notebooks are end-points, not libraries.
  The one exception is a `.ipynb` that wraps its companion `.py` script via
  `from <script_name> import ...`.
- **Reference notebooks** (e.g. `DQDdiagonalized.ipynb`, `DQDlabframe.ipynb`)
  may inline physics code without using `helpers` — they serve as independent
  cross-checks of the helper implementations and should not be refactored away.

### Naming conventions

| Pattern | Meaning |
|---------|---------|
| `<topic>.ipynb` | Interactive notebook (primary deliverable) |
| `<topic>.py` | Runnable script, often the backend for a companion `.ipynb` |

---

## Adding a new experiment

1. Write the experiment as `notebooks/<topic>.py` or directly as `notebooks/<topic>.ipynb`.
2. Import only from `helpers` and standard libraries.
3. If the same physics logic is useful elsewhere, move it into the right `helpers` module and re-import.
