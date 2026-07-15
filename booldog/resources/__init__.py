'''Packaged static resource files used by `booldog` (e.g. a Cytoscape visual
style and a matplotlib stylesheet), located on disk relative to this module
rather than via hardcoded paths, so they are found correctly regardless of
where the package is installed.
'''

import os

# visualisation
_CYTOSCAPE_STYLE_XML = "cytoscape_style.xml"
_MPL_STYLE_SHEET = "stylesheet.mplstyle"


def get_resource_file_path(file):
    '''Get the absolute path to a resource file packaged alongside this module.

    Parameters
    ----------
    file : str
        File name of the resource, relative to the `booldog/resources`
        directory (e.g. "cytoscape_style.xml").

    Returns
    -------
    str
        Absolute path to `file` on disk.
    '''
    module_path = os.path.abspath(__file__)
    resource_path = os.path.join(os.path.dirname(module_path), file)
    return resource_path


cytoscape_style_xml = get_resource_file_path(_CYTOSCAPE_STYLE_XML)
'''str : absolute path to the packaged Cytoscape visual style XML file, used
by `booldog.io.cytoscape` and `booldog.simulation_result.boolean_result` when
exporting/plotting networks in Cytoscape.'''
mpl_style_sheet = get_resource_file_path(_MPL_STYLE_SHEET)
'''str : absolute path to the packaged matplotlib stylesheet, used by
`booldog.simulation_result.continuous_result` when plotting simulation
results.'''
