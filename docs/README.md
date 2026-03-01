# BoolDoG documentation

## Visit documentation

[nib-si.github.io/BoolDoG](https://nib-si.github.io/BoolDoG)

## To build

Uses sphinx and sphinx-rtd-theme, and numpydocs for styling.

    cd docs
    conda create -n docs python pip
    conda activate docs
    pip install -r requirements.txt
    pip install -r requirements-no-deps.txt --no-deps

Then to make/remake docs:

    cd docs
    conda activate docs
    sphinx-build -b html -n -T source/ build/

To push remade docs to GitHub pages, after reviewing the HTML pages, copy them to the gh-pages folder:

    cp -r build/* gh-pages

(Then commit them and push.)

## Style guide

Follow numpy documentation guidelines

https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html#example-numpy
