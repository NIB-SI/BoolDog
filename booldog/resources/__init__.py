'''skm_tools resources '''

import os

# visualisation
_CYTOSCAPE_STYLE_XML = "cytoscape_style.xml"
_MPL_STYLE_SHEET = "stylesheet.mplstyle"


def get_resource_file_path(file):
    '''file path on disk '''
    module_path = os.path.abspath(__file__)
    resource_path = os.path.join(os.path.dirname(module_path), file)
    return resource_path


cytoscape_style_xml = get_resource_file_path(_CYTOSCAPE_STYLE_XML)
mpl_style_sheet = get_resource_file_path(_MPL_STYLE_SHEET)
