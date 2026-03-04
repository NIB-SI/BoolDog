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


def ensure_pandoc_installed(_):
    try:
        import pypandoc
        pandoc_dir = os.path.join(DOCS_DIRECTORY, "bin")
        if pandoc_dir not in os.environ["PATH"].split(os.pathsep):
            os.environ["PATH"] += os.pathsep + pandoc_dir
        pypandoc.ensure_pandoc_installed(
            targetfolder=pandoc_dir,
            delete_installer=True,
        )
    except Exception:
        pass  # pandoc already on PATH (CI) or download failed (local SSL) — continue


def setup(app):
    app.connect("builder-inited", ensure_pandoc_installed)
    # nbsphinx omits image/gif from its MIME-type pipeline, so animated GIFs
    # produced by display(IPython.display.Image(...)) fall back to the
    # text/plain repr.  Three patches are needed:
    #
    #  1. DISPLAY_DATA_PRIORITY_HTML — tell nbsphinx to *select* image/gif output.
    #  2. RST_TEMPLATE — tell the Jinja template to render image/gif via the
    #     standard ``.. image::`` directive (same as image/png / image/jpeg).
    #  3. ExtractOutputPreprocessor.extract_output_types — tell nbconvert to
    #     actually extract the GIF bytes to a file so the directive has a path.
    import logging
    _logger = logging.getLogger(__name__)
    try:
        import nbsphinx as _nbsphinx

        # 1. Priority list
        if "image/gif" not in _nbsphinx.DISPLAY_DATA_PRIORITY_HTML:
            _priority = list(_nbsphinx.DISPLAY_DATA_PRIORITY_HTML)
            try:
                idx = _priority.index("image/jpeg")
                _priority.insert(idx + 1, "image/gif")
            except ValueError:
                import warnings
                warnings.warn(
                    "whippersnappy/conf.py: 'image/jpeg' not found in "
                    "nbsphinx.DISPLAY_DATA_PRIORITY_HTML; appending "
                    "'image/gif' at the end instead. The GIF may still "
                    "render correctly, but priority ordering is unknown.",
                    stacklevel=2,
                )
                _priority.append("image/gif")
            _nbsphinx.DISPLAY_DATA_PRIORITY_HTML = tuple(_priority)

        # 2. RST template — add image/gif alongside the other raster types.
        #    We verify the substitution actually changed something; if the
        #    upstream template text has changed, warn so the breakage is visible
        #    rather than silently producing a broken GIF rendering.
        import warnings as _warnings
        _RST_OLD = "datatype in ['image/svg+xml', 'image/png', 'image/jpeg', 'application/pdf']"
        _RST_NEW = "datatype in ['image/svg+xml', 'image/png', 'image/jpeg', 'image/gif', 'application/pdf']"
        _patched_template = _nbsphinx.RST_TEMPLATE.replace(_RST_OLD, _RST_NEW)
        if _patched_template == _nbsphinx.RST_TEMPLATE:
            _warnings.warn(
                "whippersnappy/conf.py: could not patch nbsphinx.RST_TEMPLATE "
                "to add 'image/gif' support — the expected substring was not "
                "found. Animated GIFs may not render in the documentation. "
                "The nbsphinx template may have changed upstream; please update "
                "the patch in doc/conf.py.",
                stacklevel=2,
            )
        else:
            _nbsphinx.RST_TEMPLATE = _patched_template

        # 3. nbconvert extractor — ExtractOutputPreprocessor hard-codes
        #    {"image/png", "image/jpeg", "application/pdf"} as the types that
        #    get base64-decoded to binary.  image/gif falls into the "else: text"
        #    branch and is written as a raw base64 string, producing a corrupt
        #    file.  Two sub-patches fix this:
        #      3a  add "image/gif" to extract_output_types so the extractor
        #          visits it at all.
        #      3b  wrap preprocess_cell to strip gif data before the parent runs,
        #          then decode it to bytes and inject the result into resources.
        from nbconvert.preprocessors import ExtractOutputPreprocessor as _EOP
        from binascii import a2b_base64 as _a2b

        # 3a — register image/gif in extract_output_types via __init__ patch
        _eop_orig_init = _EOP.__init__
        def _eop_patched_init(self, *args, **kwargs):
            _eop_orig_init(self, *args, **kwargs)
            self.extract_output_types = self.extract_output_types | {"image/gif"}
        _EOP.__init__ = _eop_patched_init

        # 3b — patch preprocess_cell to handle image/gif as binary (base64 → bytes).
        #      The parent hard-codes only png/jpeg/pdf for binary decode; gif falls
        #      into the "else: text" branch and would be written as a raw base64
        #      string (corrupt file).  We strip gif data before calling the parent,
        #      decode it ourselves, and store the binary bytes in resources.
        _eop_orig_preprocess_cell = _EOP.preprocess_cell

        def _eop_patched_preprocess_cell(self, cell, resources, cell_index):
            # Before calling the original, convert any image/gif from base64
            # string to bytes — but the original then hits the
            # `not isinstance(data, str)` → json branch for bytes, so we must
            # pre-decode AND bypass the parent entirely for image/gif outputs.
            #
            # Strategy: strip image/gif from outputs before calling parent,
            # then handle extraction ourselves, then restore.
            import os as _os
            gif_extractions = []  # list of (out, raw_b64) to process after parent

            for out in cell.get("outputs", []):
                if out.get("output_type") not in ("display_data", "execute_result"):
                    continue
                data = out.get("data", {})
                if "image/gif" in data and isinstance(data["image/gif"], str):
                    gif_extractions.append((out, data.pop("image/gif")))

            # Run original preprocessor (without image/gif in data)
            cell, resources = _eop_orig_preprocess_cell(self, cell, resources, cell_index)

            if not gif_extractions:
                return cell, resources

            # Now handle image/gif extractions ourselves
            unique_key = resources.get("unique_key", "output")
            output_files_dir = resources.get("output_files_dir", None)
            if not isinstance(resources.get("outputs"), dict):
                resources["outputs"] = {}

            outputs_list = cell.get("outputs", [])
            for out, raw_b64 in gif_extractions:
                # Restore the b64 string in the cell data for the RST template
                out["data"]["image/gif"] = raw_b64
                # Find the index of this output in the cell
                try:
                    index = outputs_list.index(out)
                except ValueError:
                    index = 0
                # Build filename
                filename = self.output_filename_template.format(
                    unique_key=unique_key,
                    cell_index=cell_index,
                    index=index,
                    extension=".gif",
                )
                if output_files_dir is not None:
                    filename = _os.path.join(output_files_dir, filename)
                # Store binary GIF bytes in resources
                resources["outputs"][filename] = _a2b(raw_b64)
                # Store filename in output metadata so the Jinja template uses it
                if "metadata" not in out:
                    out["metadata"] = {}
                if "filenames" not in out["metadata"]:
                    out["metadata"]["filenames"] = {}
                out["metadata"]["filenames"]["image/gif"] = filename

            return cell, resources

        _EOP.preprocess_cell = _eop_patched_preprocess_cell
    except ImportError as exc:
        _logger.warning(
            "conf.py: could not patch nbsphinx/nbconvert for GIF support "
            "(package not installed): %s. Animated GIFs will not render.",
            exc,
        )
    except AttributeError as exc:
        _logger.warning(
            "conf.py: nbsphinx or nbconvert API has changed and the GIF patch "
            "could not be applied: %s. Animated GIFs will not render. "
            "Please update the patch in doc/conf.py.",
            exc,
        )