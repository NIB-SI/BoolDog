# BoolDoG documentation

## Visit documentation

[nib-si.github.io/BoolDoG](https://nib-si.github.io/BoolDoG)

## Building Locally

Docs are built with Sphinx (sphinx-rtd-theme, sphinx-autoapi, nbsphinx). These are
declared as the `docs` optional dependency group in the project's `pyproject.toml`,
alongside the other optional extras (`sbml`, `networks`, `graphviz`, `biomodels`).
Install them into a virtual environment of your choice, e.g. with Poetry:

    poetry install --extras docs

or with pip:

    pip install .[docs]

Then build the HTML docs (prefix commands with `poetry run` if using Poetry):

    cd docs
    make clean
    make html

This writes the built pages to `docs/build/html`. To view them locally, open `docs/build/html/index.html` in your browser.

## Publishing

Documentation is automatically built and deployed on every push to `main` that modifies files in the `docs/`, 'tutorials', or `booldog/` directories.

The GitHub Actions workflow (.github/workflows/deploy-docs.yml) handles:

1. Building the Sphinx documentation
2. Deploying to the gh-pages branch
3. Publishing to GitHub Pages at https://nib-si.github.io/BoolDog

## Style guide

Follow numpy documentation guidelines:

https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html#example-numpy
