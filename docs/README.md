# BoolDoG documentation

## Visit documentation

[nib-si.github.io/BoolDoG](https://nib-si.github.io/BoolDoG)

## To build

Uses sphinx and sphinx-rtd-theme, and numpydocs for styling.

     create -n docs python pip
     conda activate docs
     pip install -r requirements.txt 
     pip install -r requirements-no-deps.txt --no-deps

Then to make/remake docs:

    cd docs
    sphinx-apidoc -fMeT  --module-first -o source/api ../booldog ../booldog/utils
    make html

To push remade docs to github pages site, copy the new html to the gh-pages folder

    cp -r build/html/* gh-pages

(Then commit them and push.)

## Style guide

Follow numpy documentation guidelines

https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html#example-numpy
