# BoolDoG documentation

## Visit documentation

[nib-si.github.io/BoolDoG](https://nib-si.github.io/BoolDoG)

## To build

Uses sphinx and sphinx-rtd-theme, and numpydocs for styling.

    conda install sphinx numpydocs
    pip install sphinx-rtd-theme

Then to make/remake docs:

    cd docs
    sphinx-apidoc -fMeT  --module-first -o source/api ../booldog ../booldog/utils
    make html
    cp -r build/html/* gh-pages

## Style guide

Follow numpy documentation guidelines

https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html#example-numpy
