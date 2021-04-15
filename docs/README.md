# squad-reboot documentation

## Visit documentation

[nib-si.github.io/squad_reboot](https://nib-si.github.io/squad_reboot)

## To build

Uses sphinx and insegel theme 

    conda install sphinx
    pip install insegel

Then to make/remake docs: 

    cd docs
    sphinx-apidoc -o source ../squad_reboot 
    make html
    cp -r build/html/* gh-pages
