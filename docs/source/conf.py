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

# -- Project information -----------------------------------------------------

project = 'BoolDog'
copyright = '2020-2025 National Institute of Biology, Slovenia'
author = 'Carissa Bleker'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'autoapi.extension',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_rtd_theme',
    'sphinx.ext.todo',
    'sphinx.ext.intersphinx',

    # embed jupyter notebooks in the documentation
    "nbsphinx",
    'nbsphinx_link',
]

# Autoapi settings
autoapi_dirs = ['../../booldog']


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

# def autodoc_skip_member(app, what, name, obj, skip, options):

#     # print(app, what, name, obj, skip, options)
#     print("HELLLO", name, obj)

#     excludes = ['booldog.io.read']

#     exclude = None
#     for ex in excludes:
#         this_exclude = re.findall(f'.*{ex}.*', str(obj))
#         if this_exclude:
#             exclude = True
#             print("SKIP ME!!!", str(obj))
#             break

#     return exclude

# def setup(app):

#     app.connect('autodoc-skip-member', autodoc_skip_member)

#     # i.o.t add ODE class (inside a function) to the documentation
#     # import booldog
#     # example = booldog.ode_factory.ode_factory(
#     #             '',
#     #             transform='placeholder')
#     # booldog.ode_factory.ODE = booldog.ode_factory.ode_factory.ex_class
#     # booldog.ode_factory.ODE.__name__ = 'ODE'
#     # booldog.ode_factory.ODE.__module__ = 'booldog.ode'
