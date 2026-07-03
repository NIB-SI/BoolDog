# CI

`tests.yml` runs on push to `main` and on pull requests, but only when a
change touches dependencies (`pyproject.toml`), `booldog/`, `tutorials/`,
`tests/`, or the workflow/scripts themselves — see the `paths:` filter in
the workflow for the exact list. Two jobs, each matrixed across
`{ubuntu, macos, windows}-latest` x Python `{3.12, 3.13}`:

- **test**: installs `.[networks,sbml,biomodels]` and runs
  `python -m unittest discover .` in `tests/`.
- **tutorials**: matrixed additionally over `notebook: [tutorial-basic,
  tutorial-advanced]`. Every notebook is treated the same way: strip any
  `# ci-skip: ...` cells (see below), then execute the result headlessly with
  `jupyter execute`.

The Python/OS matrix is limited to 3.12/3.13 because those are the versions
with full prebuilt-wheel coverage (no compiler toolchain needed) for every
dependency across all three OSes, including `igraph==0.11.9`, which is pinned
exactly due to an upstream GraphML edge-id regression in `igraph>=1.0.0`
(see `booldog/io/interaction_networks.py`). Extend the matrix once newer
Python versions have equivalent wheel coverage.

macOS legs of both jobs are `continue-on-error: true`. GitHub's
`macos-latest` runners are Apple Silicon, and pyboolnet ships an Intel-only
`BNetToPrime` binary that crashes on any `BoolDogModel` instantiation there
([NIB-SI/BoolDog#6](https://github.com/NIB-SI/BoolDog/issues/6)). That's an
upstream pyboolnet issue, not fixable here — macOS stays in the matrix for
visibility (and so it turns green automatically once upstream is fixed),
but doesn't block the pipeline.

## Skipping cells that need external/interactive tools

Some notebook cells can't run in a headless CI runner — e.g.
`tutorial-advanced.ipynb` drives a live Cytoscape Desktop instance via
`py4cytoscape`, which isn't available in CI. Mark such cells by making
`# ci-skip: <reason>` the first line of the cell's source. CI strips any
matching cell (`.github/scripts/strip_ci_skip_cells.py`) before execution;
the tracked notebook itself is untouched, so it still renders normally in
Jupyter, GitHub, and the docs build.

## Adding a new tutorial notebook

Add its basename (without `.ipynb`) to the `notebook` matrix list in
`tests.yml`. It'll automatically go through the same strip-then-execute
pipeline as the existing tutorials.

## Reproducing locally

```bash
pip install .[networks,sbml,biomodels] nbclient ipykernel

# tests
cd tests && python -m unittest discover .

# tutorials
python .github/scripts/strip_ci_skip_cells.py \
    tutorials/tutorial-basic.ipynb tutorials/.ci-tutorial-basic.ipynb
cd tutorials && jupyter execute --kernel_name=python3 --timeout=300 .ci-tutorial-basic.ipynb
```
