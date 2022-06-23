# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import re
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'BoolDoG'
copyright = '2020-2022 National Institute of Biology, Slovenia'
author = 'Carissa Bleker'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_rtd_theme',
    'sphinx.ext.todo',
    'sphinx.ext.intersphinx',
]

# Autodoc settings
autodoc_default_options = {
    'inherited-members':True,
    'undoc-members':True,
    'private-members':False,
    'special-members': '__init__',
    'exclude-members': '__weakref__'
}
autoclass_content = 'both'
autodoc_member_order = 'bysource'


# autosummary_mock_imports = [
#     'booldog.io.read',
# ]


# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True


# Intersphinx settings
intersphinx_mapping = {
    'pyboolnet': ('https://pyboolnet.readthedocs.io/en/master/', None),
}

# ----------------------------------------------------------------------------

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'style_nav_header_background':'#009739'
}

html_logo = '../figures/logo.png'
html_favicon = '../figures/icon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

def autodoc_skip_member(app, what, name, obj, skip, options):

    # print(app, what, name, obj, skip, options)
    print("HELLLO", name, obj)

    excludes = ['booldog.io.read']

    exclude = None
    for ex in excludes:
        this_exclude = re.findall(f'.*{ex}.*', str(obj))
        if this_exclude:
            exclude = True
            print("SKIP ME!!!", str(obj))
            break

    return exclude

def setup(app):

    app.connect('autodoc-skip-member', autodoc_skip_member)

    # i.o.t add ODE class (inside a function) to the documentation
    import booldog
    example = booldog.ode.ODE_factory(
                '',
                transform='placeholder')
    booldog.ode.ODE = booldog.ode.ODE_factory.ex_class
    booldog.ode.ODE.__name__ = 'ODE'
    booldog.ode.ODE.__module__ = 'booldog.ode'
