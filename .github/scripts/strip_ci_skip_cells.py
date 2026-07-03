"""Write a copy of a notebook with '# ci-skip: ...' cells removed.

Used to execute tutorial notebooks headlessly in CI, skipping cells that
depend on external/interactive tools (e.g. a running Cytoscape Desktop
instance) that CI runners can't provide.

Notebooks that call the live BioModels API also get a hidden setup cell
prepended, so CI replays a recorded cassette instead of depending on
BioModels' availability (same pattern as tests/test_biomodels.py). This
is injected only into the throwaway CI copy, never the tracked notebook -
a reader/learner running the real notebook never sees it or needs vcrpy
installed.
"""
import sys

import nbformat

MARKER = "# ci-skip:"

BIOMODELS_CASSETTE_SETUP = '''\
import os
import vcr
import booldog.io.biomodels as _biomodels

_my_vcr = vcr.VCR(
    cassette_library_dir="../tests/cassettes",
    record_mode="none" if os.environ.get("CI") == "true" else "once",
)
_biomodels.fetch_model_info = _my_vcr.use_cassette(
    "notebook_fetch_model_info.yaml")(_biomodels.fetch_model_info)
_biomodels.fetch_model = _my_vcr.use_cassette(
    "notebook_fetch_model.yaml")(_biomodels.fetch_model)
'''


def main(src_path, dst_path):
    nb = nbformat.read(src_path, as_version=4)
    nb.cells = [
        cell for cell in nb.cells
        if not (cell.cell_type == "code"
                and cell.source.lstrip().startswith(MARKER))
    ]

    uses_biomodels = any(
        cell.cell_type == "code" and "booldog.io.biomodels" in cell.source
        for cell in nb.cells)
    if uses_biomodels:
        nb.cells.insert(0, nbformat.v4.new_code_cell(BIOMODELS_CASSETTE_SETUP))

    nbformat.write(nb, dst_path)


if __name__ == "__main__":
    main(*sys.argv[1:3])
