def get_edges(position_table: dict, use_orientation: bool = False) -> list:
    """
    Get the edges present in an assembly based on the position table of the 
    adjacency graph. 

    Args: 
        position_table: The position table from the create_table() step. 
        use_orientation: Boolean, we can reuse this function between ordering 
            and orientation. 
    
    Returns:
        A list of tuples which describe all contig adjacencies in an assembly or
        reference. 
    """
    edges = []
    for index in position_table:
        edge = set()
        for contig, orientation in position_table[index]:
            if use_orientation:
                contig = (contig, orientation)
            edge.add(contig)

        # Skip edges with telomeres.
        if len(edge) < 2:
            continue

        edges.append(edge)

    return edges


def count(first: tuple, second: tuple, orientation: bool = False) -> float:
    """
    Compute the length weighted intersection of two edge lists. 

    Args: 
        first: The reference genome as (contig_table, position_table). 
        second: The assembly as (contig_table, position, table).
    
    Returns:
        A fraction which represents how many contigs were correctly ordered or
        oriented. 
    """
    first_contigs, first_positions = first
    second_contigs, second_positions = second

    first_connections = get_edges(first_positions, orientation)
    second_connections = get_edges(second_positions, orientation)

    weighted_first = 0
    weighted_second = 0
    for edge in first_connections:

        # Add the edge weight.
        current_weight = 0
        for contig in edge:
            if orientation:
                contig = contig[0]
            current_weight += first_contigs[contig]["length"]
        weighted_first += current_weight

        # If the edge is in the other assembly, count it and weight it.
        if edge in second_connections:
            weighted_second += current_weight

    return weighted_second / weighted_first
