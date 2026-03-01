import py4cytoscape as p4c


def export_network(suid,
                   filename,
                   format="png",
                   overwrite_file=True,
                   zoom=200):

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
