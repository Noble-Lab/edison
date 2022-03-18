def count(first: tuple, second: tuple) -> float:
    """
    Since scaffolds are not directly stored and instead are associated with each
    contig, we must extract the total length of scaffolds in each assembly, as 
    well as the length of any intersections.

    Then, we get the maximum weighted Jaccard index for each reference scaffold, 
    which is defined as the length of the intersection divided by the length of 
    the union. To weight this in the averaging step, we then multiply this by 
    the length of the reference scaffold. 

    The sum of these maximum weighted indices are then divided by the total 
    length of the reference scaffolds. 

    Args: 
        first: The reference genome as (contig_table, position_table). 
        second: The assembly as (contig_table, position, table).
    
    Returns:
        A fraction which represents how many contigs were correctly grouped.
    """

    first_contigs, first_positions = first
    second_contigs, second_positions = second

    intersection_scaffolds = {}
    first_scaffolds = {}
    second_scaffolds = {}
    for contig in first_contigs:

        # Get contig length and scaffold information.
        contig_length = first_contigs[contig]["length"]
        first_scaffold_name = first_contigs[contig]["scaffold"]
        second_scaffold_name = second_contigs[contig]["scaffold"]

        # Iterate the count and length on the scaffold of the first assembly.
        if first_scaffold_name not in first_scaffolds:
            first_scaffolds[first_scaffold_name] = contig_length
        else:
            first_scaffolds[first_scaffold_name] += contig_length

        # Iterate the count and length on the scaffold of the second assembly.
        if second_scaffold_name not in second_scaffolds:
            second_scaffolds[second_scaffold_name] = contig_length
        else:
            second_scaffolds[second_scaffold_name] += contig_length

        # Iterate the count of the intersection.
        intersection = (first_scaffold_name, second_scaffold_name)
        if intersection not in intersection_scaffolds:
            intersection_scaffolds[intersection] = contig_length
        else:
            intersection_scaffolds[intersection] += contig_length

    weighted_jaccard = []
    for i in first_scaffolds:

        maximum = 0
        for j in second_scaffolds:

            # Get a value for the intersection.
            if (i, j) not in intersection_scaffolds:
                continue

            # Get a value for the union.
            weighted_union = (
                first_scaffolds[i]
                + second_scaffolds[j]
                - intersection_scaffolds[(i, j)]
            )

            # Append the Jaccard index.
            weighted_index = (
                intersection_scaffolds[(i, j)] / weighted_union
            ) * first_scaffolds[i]

            if weighted_index > maximum:
                maximum = weighted_index

        weighted_jaccard.append(maximum)

    return sum(weighted_jaccard) / sum(first_scaffolds.values())
