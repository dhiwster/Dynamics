Python layout
=============

Source packages:

- `helpers/`
  Shared Hamiltonian, operator, pulse, solver, and DQD support code.
- `notebooks/`
  Exploratory notebooks and ad hoc experiment scripts.

Current rule of thumb:

- Put reusable simulation code in `helpers/`.
- Keep problem-specific experiments in notebooks or dedicated scripts unless they
  clearly belong in the shared helper layer.
- Keep notebooks and one-off experiment scripts out of the source packages.
