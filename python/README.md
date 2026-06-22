Python layout
=============

Source packages:

- `hamiltonians/`
  Shared Hamiltonian helpers, models, and reusable builders.
- `dqd/`
  Double-quantum-dot system, builders, circuit tooling, and solvers.
- `notebooks/`
  Exploratory notebooks and ad hoc experiment scripts.

Current rule of thumb:

- Put reusable Hamiltonian/model helpers in `hamiltonians/`.
- Put problem-specific implementations in their own top-level package such as `dqd/`.
- Keep notebooks and one-off experiment scripts out of the source packages.
