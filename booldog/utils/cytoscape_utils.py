'''Helpers for exporting Cytoscape network views to image files, and for
extracting node-induced subgraphs, via `py4cytoscape`. Requires a running
Cytoscape instance (with the CyREST API reachable) to talk to.
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


def extract_subnetwork(suid,
                       seed_nodes,
                       boundary_node_types=("species",),
                       direction="predecessors",
                       node_type_column="type",
                       subnetwork_name=None,
                       verify=True):
    '''Extract the node set (and, optionally, a Cytoscape subnetwork) of
    the *hyperedge-closed* neighbourhood of `seed_nodes` in a BoolDog
    logic-circuit network.

    A `logic_circuit`-style network (see
    `BoolDogModelIOToMixin.to_networkx` with ``as_logic_circuit=True``)
    represents each Boolean rule as species nodes feeding into
    ``logical_and``/``logical_or``/``logical_not`` gate nodes, which in
    turn feed the regulated species. A gate's set of incoming edges is a
    single hyperedge: all of its inputs are one unit of meaning, so a
    naive "N-hop neighbours" traversal can silently include only some of
    a gate's inputs and change what the extracted subgraph actually
    means. This function instead walks node-by-node, and whenever it
    reaches a gate node it always includes *all* of that gate's inputs
    (its whole incoming hyperedge, also when walking downstream) as well
    as following its edges in `direction` - so no hyperedge is ever
    partially included - and only stops expanding at nodes whose type is in
    `boundary_node_types` (typically ``"species"``, i.e. regulator
    species one hop upstream of a seed, in their entirety as single
    nodes rather than expanding into their own regulators).

    Parameters
    ----------
    suid : int or str
        The Cytoscape network SUID (or network name) to read from.
    seed_nodes : iterable of str
        Names of the node(s) to grow the subgraph from (e.g. the species
        whose regulators/targets should be included). Always fully
        expanded regardless of their own type, even if that type is
        also in `boundary_node_types`.
    boundary_node_types : iterable of str, default ("species",)
        Node types (per `node_type_column`) at which traversal stops:
        such a node is included in the result but its own
        predecessors/successors are not walked into. Compare against
        the network's ``"type"`` node column, e.g.
        ``{"species", "logical_and", "logical_or", "logical_not"}``.
    direction : {"predecessors", "successors"}, default "predecessors"
        Whether to walk upstream (regulators of `seed_nodes`, following
        edges backward) or downstream (what `seed_nodes` regulate,
        following edges forward).
    node_type_column : str, default "type"
        Name of the node table column identifying each node's type,
        checked against `boundary_node_types`.
    subnetwork_name : str or None, optional
        If given, select the extracted node set and create a Cytoscape
        subnetwork from it (via `py4cytoscape.networks.create_subnetwork`)
        named `subnetwork_name`, clearing any prior selection first.
        Default `None` (no subnetwork is created; `suid`'s own selection
        is left untouched).
    verify : bool, default True
        If `True`, after the walk, independently check for every
        non-boundary (gate-type) node in the result that the sources of
        *all* of its incoming edges in the full, unfiltered edge table
        are also present in the result - i.e. that no hyperedge was
        partially included, whatever the `direction`. Raises
        `ValueError` naming the first violation found if this check
        fails; such a violation would indicate a bug in this function
        rather than anything about the input network.

    Returns
    -------
    node_names : set of str
        Names of all nodes in the extracted subgraph (seeds, gate nodes,
        and boundary-type nodes reached along the way).
    subnetwork_suid : int, optional
        SUID of the created subnetwork; only returned if
        `subnetwork_name` was given.

    Raises
    ------
    ValueError
        If `direction` is not one of ``"predecessors"``/``"successors"``,
        or (when `verify` is `True`) if a hyperedge was found to be only
        partially included in the result.
    '''
    if direction not in ("predecessors", "successors"):
        raise ValueError(
            f"direction must be 'predecessors' or 'successors', got {direction!r}."
        )

    node_table = p4c.tables.get_table_columns(network=suid)
    edge_table = p4c.tables.get_table_columns(table="edge", network=suid)

    name_to_type = dict(zip(node_table["name"], node_table[node_type_column]))

    incoming = {}
    outgoing = {}
    for source, target in zip(edge_table["source"], edge_table["target"]):
        incoming.setdefault(target, []).append(source)
        outgoing.setdefault(source, []).append(target)
    neighbours = incoming if direction == "predecessors" else outgoing

    seed_nodes = set(seed_nodes)
    visited = set()

    def _walk(start):
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            is_boundary = name_to_type.get(node) in boundary_node_types
            if is_boundary and node not in seed_nodes:
                continue
            next_nodes = list(neighbours.get(node, []))
            if not is_boundary:
                # a gate's inputs form one hyperedge: include all of them,
                # also when walking downstream
                next_nodes += incoming.get(node, [])
            stack.extend(n for n in next_nodes if n not in visited)

    for seed in seed_nodes:
        _walk(seed)

    if verify:
        incomplete = _find_incomplete_hyperedge(visited, edge_table,
                                                name_to_type,
                                                boundary_node_types)
        if incomplete is not None:
            node, missing = incomplete
            raise ValueError(
                f"Incomplete hyperedge: {node!r} has an incoming "
                f"edge from {missing!r}, which was not included "
                "in the extracted subgraph. This indicates a bug "
                "in extract_subnetwork rather than the input "
                "network."
            )

    if subnetwork_name is not None:
        p4c.network_selection.clear_selection(type="both", network=suid)
        p4c.network_selection.select_nodes(list(visited),
                                           by_col="name",
                                           network=suid)
        subnetwork_suid = p4c.networks.create_subnetwork(
            subnetwork_name=subnetwork_name, network=suid)
        return visited, subnetwork_suid

    return visited


def _find_incomplete_hyperedge(node_names, edge_table, name_to_type,
                               boundary_node_types):
    '''Find a gate in `node_names` whose inputs are only partly included.

    Built directly from the full edge table (independent of the walk in
    `extract_subnetwork`), this checks every node in `node_names` whose
    type is not in `boundary_node_types` for incoming edges from nodes
    outside `node_names`.

    Parameters
    ----------
    node_names : set of str
        Names of the nodes in the extracted subgraph.
    edge_table : pandas.DataFrame
        The full edge table, with ``"source"`` and ``"target"`` columns of
        node names.
    name_to_type : dict
        Mapping of node name to node type.
    boundary_node_types : iterable of str
        Node types that are not checked (see `extract_subnetwork`).

    Returns
    -------
    tuple of (str, str) or None
        ``(node, missing_source)`` for the first incomplete hyperedge
        found, or `None` if every checked node has all of its inputs in
        `node_names`.
    '''
    for source, target in zip(edge_table["source"], edge_table["target"]):
        if (target in node_names and source not in node_names
                and name_to_type.get(target) not in boundary_node_types):
            return target, source
    return None
