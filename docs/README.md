# BoolDoG documentation

## Visit documentation

[nib-si.github.io/BoolDoG](https://nib-si.github.io/BoolDoG)

## To build

Docs are built with Sphinx (sphinx-rtd-theme, sphinx-autoapi, nbsphinx). These are
declared as the `docs` optional dependency group in the project's `pyproject.toml`,
alongside the other optional extras (`sbml`, `networks`, `graphviz`, `biomodels`).
Install them into a virtual environment of your choice, e.g. with Poetry:

    poetry install --extras docs

or with pip:

    pip install .[docs]

Then build the HTML docs (prefix commands with `poetry run` if using Poetry):

    cd docs
    make html

This writes the built pages to `docs/build/html`. To rebuild from scratch (clearing
Sphinx's cached doctrees first):

    make clean
    make html

## Publish

To push rebuilt docs to GitHub Pages, after reviewing the HTML output, mirror
`build/html` into the tracked `gh-pages` folder (served via GitHub Pages from `/docs`
on `main`). Use `rsync --delete` rather than `cp` so pages removed or renamed since the
last publish don't linger:

    rsync -a --delete build/html/ gh-pages/

Preview what would change first with `rsync -an --delete build/html/ gh-pages/`.

Then commit and push the changes to `docs/gh-pages`.

## Style guide

Follow numpy documentation guidelines

https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html#example-numpy
