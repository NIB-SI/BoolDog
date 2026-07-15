'''Helpers for exporting Cytoscape network views to image files via
`py4cytoscape`. Requires a running Cytoscape instance (with the CyREST API
reachable) to talk to.
'''
import py4cytoscape as p4c


def export_network(suid,
                   filename,
                   format="png",
                   overwrite_file=True,
                   zoom=200):
    '''Export a Cytoscape network view to an image file.

    Fits the view to its content, clears any node/edge selection (so
    selection highlighting does not show up in the exported image), and
    exports the current view of the network identified by *suid* via
    `py4cytoscape.network_views.export_image`.

    Parameters
    ----------
    suid : int or str
        The Cytoscape network SUID (or, per `py4cytoscape` convention, a
        network name) identifying the network/view to export.
    filename : str or Path
        Path of the image file to write. Coerced to `str` before being
        passed on to `py4cytoscape`.
    format : str, default "png"
        Image format to export, passed as `py4cytoscape`'s ``type``
        argument (e.g. "png", "pdf", "svg", ...).
    overwrite_file : bool, default True
        Whether to overwrite *filename* if it already exists.
    zoom : int, default 200
        Zoom percentage applied to the exported image.

    Returns
    -------
    None

    Notes
    -----
    `py4cytoscape.network_views.fit_content` fits the view to node
    positions only: it ignores edges/edge labels that extend beyond node
    boundaries, so such elements can be clipped in the exported image.
    This is a known Cytoscape limitation (reported as CSD-979); the
    Cytoscape developers' response was "won't fix".
    '''

    filename = str(filename)

    # fit_content ignores edges that may extend beyond node boundaries
    # Reported bug: CSD-979 - response from Cytoscape dev team is WONT FIX
    p4c.network_views.fit_content(network=suid)
    p4c.network_selection.clear_selection(type='both', network=suid)

    p4c.network_views.export_image(filename=filename,
                                   type=format,
                                   network=suid,
                                   overwrite_file=overwrite_file,
                                   zoom=zoom)
