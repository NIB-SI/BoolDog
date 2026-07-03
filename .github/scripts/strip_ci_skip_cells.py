"""Write a copy of a notebook with '# ci-skip: ...' cells removed.

Used to execute tutorial notebooks headlessly in CI, skipping cells that
depend on external/interactive tools (e.g. a running Cytoscape Desktop
instance) that CI runners can't provide.
"""
import sys

import nbformat

MARKER = "# ci-skip:"


def main(src_path, dst_path):
    nb = nbformat.read(src_path, as_version=4)
    nb.cells = [
        cell for cell in nb.cells
        if not (cell.cell_type == "code"
                and cell.source.lstrip().startswith(MARKER))
    ]
    nbformat.write(nb, dst_path)


if __name__ == "__main__":
    main(*sys.argv[1:3])
