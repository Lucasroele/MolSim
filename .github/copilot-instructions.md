# MolSim repository guidance

## Project overview
- This repository contains a Python package named `molsim` built with a `src/` layout.
- The main CLI entry point is `molsim` and is registered in `src/molsim/cli.py`.
- Most functionality is implemented as subcommands under `src/molsim/commands/`, with shared helpers in `src/molsim/utils/`, `src/molsim/parsers/`, and top-level package modules.
- Tests live under `tests/` and use `pytest`.

## Working conventions
- Prefer small, focused functions and keep CLI wrappers thin. Business logic should live in reusable functions/modules rather than inside argparse callbacks.
- When adding or changing a command, update the command registration in `src/molsim/cli.py` and keep the subcommand naming consistent with the existing CLI style.
- Preserve the project’s argument parsing style: commands use `argparse` and are registered via `register(subparsers)` patterns.
- Keep output and file handling explicit and predictable; many commands read/write `.xvg`, `.top`, `.gro`, and similar simulation file formats.
- Favor readable, deterministic code over clever abstraction. This project is a tooling repo for molecular simulation prep rather than a large application framework.

## Repo structure
- `src/molsim/cli.py`: root CLI setup and command registration.
- `src/molsim/commands/`: individual CLI commands.
- `src/molsim/parsers/`: parsers for file formats such as `.xvg`.
- `src/molsim/utils/`: shared helpers.
- `tests/fixtures/`: input fixture files for command tests.
- `tests/molsim/`: pytest coverage for commands and helpers.

## Validation steps
- Use the project’s installed package workflow: `pip install -e .` from the repository root.
- Do not ask the user for permission to run repo-local commands such as `cd`, `pytest`, or other validation commands; work from the project root and execute the smallest relevant check directly.
- Run the most targeted test available for the area being changed, for example:
  - `pytest tests/molsim/commands/test_count_mols.py -q`
  - `pytest tests/molsim/parsers/test_xvgParser.py -q`
  - `pytest -q` for the full suite when the change is broader.
- If a CLI command behavior changes, verify it via the relevant command path and check the fixture-based tests before broad refactoring.

## Style and safety
- Keep changes consistent with the existing module layout and naming patterns.
- Do not introduce heavy frameworks or architectural churn without a clear need.
- Prefer minimal changes that address the root cause and preserve compatibility with the CLI surface.
- If the user request involves debugging a failing command, first identify the command path and the failing fixture/test; then patch the underlying parser/logic rather than only adjusting the CLI wrapper.

## Helpful defaults for this repo
- Assume Python 3.13 compatibility as declared in `pyproject.toml`.
- Use `pytest` as the default validation tool.
- When exploring the repo, start with the command registration and its corresponding tests before revising core logic.
