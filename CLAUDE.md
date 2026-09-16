# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

BoolDog (package name `booldog`) is a Python library for building, converting, and simulating
Boolean and semi-quantitative Boolean networks (systems/computational biology). It wraps
`pyboolnet` for the core Boolean algebra (prime implicants) and adds network I/O, ODE-based
continuous relaxations, and simulation-result plotting on top.

## Environment & commands

The project is managed with **uv**. `pyproject.toml` only declares a `[project]` table
(PEP 621) with a setuptools build backend — there is no `[tool.uv]` section needed for
this to work, `uv sync`/`uv run` work directly against it. Requires Python >= 3.12.

```bash
uv sync                        # install core dependencies into the uv-managed venv
uv run python -m unittest discover .   # run tests from within tests/ (see below)
```

Some tests and modules need optional extras that are not part of the base `uv sync`:
- `booldog/io/biomodels.py` (BioModels download API) needs the `biomodels` extra
  (`requests`), needed to run `tests/test_biomodels.py`.
- `tests/test_io.py` exercises the `networks` extra (`igraph`).
- SBML-qual support (`booldog/io/sbml.py`) needs the `sbml` extra (`python-libsbml`, `tabularqual`).

Install everything needed for the full test run with:
```bash
uv sync --all-extras
```

`all` in `[project.optional-dependencies]` is a flat, duplicated list rather than a
self-reference to the other extras (`"BoolDog[sbml,networks,...]"`) — this is a legacy
of a prior Poetry-based setup where that self-referencing PEP 621 pattern broke
`poetry lock`; it's harmless to keep as-is under uv.

### Tests

Tests are plain `unittest` (no pytest config), living in `tests/`. Run all of them from
inside that directory (relative paths to `tests/data` assume this cwd):
```bash
cd tests && uv run python -m unittest discover .
```
Run a single test module/case/method the standard unittest way, e.g.:
```bash
cd tests && uv run python -m unittest test_utils.Test.test_ensure_ndarray
```
`tests/deleteme_test_io.py` is scratch/leftover, not part of the intended suite.

### Lint / format

Per `CONTRIBUTING.md`: formatting uses **yapf** (`.style.yapf`, pep8-based, 79-column limit)
and linting uses **pylint** (`.pylintrc`). Neither is a declared project dependency, so
install them into the venv before use (`uv run pip install yapf pylint`).

### Docs

Sphinx + autoapi + nbsphinx docs live in `docs/source`, built with `uv sync
--extra docs && cd docs && uv run make html` (see `docs/README.md`). Output goes to
`docs/build/html`, which is published by mirroring it into `docs/gh-pages` with
`rsync -a --delete` (served via GitHub Pages from `/docs` on `main`) — don't hand-edit
`docs/gh-pages`; regenerate and re-sync instead.

## Architecture

The central class is `booldog.network.BoolDogModel`, exposed as `booldog.BoolDogModel`.
It carries almost no logic of its own — it's a composition of mixins, each implemented
in its own subpackage:

```
BoolDogModel(
    BooleanNetworkMixin,             # booldog/boolean/boolean.py
    BooleanNetworkModificationMixin, # booldog/boolean/modifications.py
    ContinuousMixin,                 # booldog/continuous/semi_quantitative.py
    BoolDogModelIOFromMixin,         # booldog/io/__init__.py  (from_bnet, from_sbmlqual, ...)
    BoolDogModelIOToMixin,           # booldog/io/__init__.py  (to_bnet, to_networkx, ...)
)
```

When looking for where a `BoolDogModel` method is implemented, check these mixins rather
than `network.py` itself. `network.py` only owns `__init__`, node/index bookkeeping
(`self.nodes`, `self.node_ids`, `self.index`), and the `primes` property.

- **`booldog/classes.py`** — plain dataclasses: `BoolDogNode` (identifier/name/rule) and
  `BoolDogModelInfo` (metadata). A network is fundamentally a dict of `BoolDogNode`s.
- **Prime implicants are the canonical internal representation** of the Boolean rules
  (see `pyboolnet`'s prime-implicant format). `BoolDogModel.primes` is lazily computed
  from node rules via `pyboolnet.file_exchange` and then cached (`_primes_cached`) —
  any mutation of node rules must invalidate this cache (see
  `BooleanNetworkModificationMixin._update_model_object`).
- **`booldog/boolean/`** — `boolean.py` has read-only Boolean-network operations
  (state-space generation, `boolean_simulation`, `steady_states`, parent/child lookups);
  `modifications.py` has structural edits (`add_node`, `remove_node(s)`, `update_node`,
  `modify_network`) plus a `Modification`/`ModificationTypes` audit trail
  (`self.modifications`) recording changes made to a model.
- **`booldog/continuous/`** — semi-quantitative/continuous relaxations of the Boolean
  model. `semi_quantitative.py`'s `ContinuousMixin` is the public entry point
  (`transform_bool_to_continuous`, `continuous_simulation`); `ode_factory.py` builds the
  actual ODE system via an `ode_factory()` function returning one of two `ODE` subclasses:
  `BooleCubeODE` (Boolean cube / Hill-function based) or `SquadODE` (SQUAD-style sigmoidal
  transform). Add new continuous transforms here.
- **`booldog/io/`** — one module per exchange format, each exposing `read_x`/`write_x`
  functions with a consistent contract: readers return a dict of `nodes`/`primes`/
  `modelinfo` consumed by `BoolDogModel(**data)`; writers take the model instance first.
  `BoolDogModelIOFromMixin`/`ToMixin` in `io/__init__.py` are thin dispatchers
  (`_from_reader`/`_to_writer`) around these functions — add a new format by writing the
  reader/writer module and registering a `from_x`/`to_x` classmethod there.
  `interaction_networks.py` + `interaction_logic.py` handle generic activator/inhibitor
  interaction tables/graphs (SIF, GraphML, networkx, igraph) and convert them into
  Boolean rules via a pluggable `LogicBuilder` (e.g. `SquadLogic`).
- **`booldog/simulation_result/`** — result objects returned by simulations, not by
  `BoolDogModel` itself: `BooleanStateSpace`/`BooleanSimulationResult` (from
  `boolean_simulation`, includes state-transition-graph plotting/animation) and
  `ContinuousSimulationResult` (from `continuous_simulation`, includes time-series plotting).
- **`booldog/utils/`** — `boolean_normal_forms.py` (functions ↔ minimal DNF via pyboolnet),
  `decorators.py` (e.g. `@validate_node_argument`, which normalizes a `BoolDogNode`-or-id
  argument and checks membership in `self.node_ids` — used throughout the mixins), `misc.py`
  (array/parameter coercion helpers, re-exported via `from .misc import *`), `logger.py`.
- **`booldog/resources/`** — packaged static assets (Cytoscape style XML, matplotlib
  stylesheet), located via `get_resource_file_path` rather than hardcoded paths.

## Non-package directories

- **`dev/`** — scratch/exploratory scripts and notebooks (SWIG bindings experiments, ad-hoc
  SBML parsing, etc.), not part of the installed package or referenced by it.
- **`squad/`** — a Dockerized build of the third-party SQUAD tool, unrelated to the Python
  package's build/test pipeline; see `squad/README.md` for Docker usage.
- **`tutorials/`** — Jupyter notebooks that double as user-facing docs, embedded into the
  Sphinx build via `nbsphinx`/`nbsphinx_link`.
