# Dynamic Simulations Workflow

This repository is the source of truth for code and coordination across machines.
Large simulation outputs are written to a shared OneDrive directory and are intentionally kept out of Git.

## Repo responsibilities

- Version-control Python scripts, Mathematica `.wl` files, and parameter/config files.
- Track simulation assignment and status in `runs.json`.
- Keep machine environments consistent via `requirements.txt`.

## OneDrive responsibilities

- Store heavy outputs such as arrays, plots, logs, and Mathematica binary artifacts.
- Sync those outputs across machines outside of Git.

Suggested shared output path on every machine:

```text
~/OneDrive/simulation_outputs/
```

If the exact filesystem location differs by OS, keep the logical target consistent and expose it through local machine config rather than hardcoding many variants in scripts.

## Recommended daily workflow

1. Start each session with `git pull`.
2. Claim or update a run entry in `runs.json`.
3. Execute simulations and write outputs to the shared OneDrive folder.
4. Commit code, `.wl`, and config changes back to this repository.

## Suggested repository layout

```text
.
|-- local_paths.example.json
|-- README.md
|-- requirements.txt
|-- runs.json
|-- python/
|-- mathematica/
`-- configs/
```

## Local machine paths

Copy `local_paths.example.json` to `local_paths.json` on each machine and adjust paths there.
`local_paths.json` is ignored by Git so every computer can point to its own OneDrive mount location while scripts keep a stable config shape.

## Notes for Mathematica

- Treat notebooks as exploratory frontends when useful.
- Export stable logic to plain-text `.wl` files for reviewable diffs and lower merge friction.

## Notes for dependencies

`requirements.txt` currently contains placeholder pins. Replace them with the exact packages your simulations need and update them deliberately when environments change.
