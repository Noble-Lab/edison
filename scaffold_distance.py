import argparse
import pathlib

# Local imports.
import grouping
import ordering
import agp


def setup():
    parser = argparse.ArgumentParser(
        description="Find the scaffolding accuracy based on two AGP files."
    )

    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        type=pathlib.Path,
        help="Reference AGP file.",
    )

    parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        type=pathlib.Path,
        help="Assembly AGP file.",
    )

    return parser.parse_args()


def get_orientation(orientation: str, left: bool) -> str:
    """
    Determine which side of the contig is head and tails. If the orientation is
    + and it is on the left, then we know the tail is involved in the junction.
    We compute all four possibilities. 

    Args:
        orientation: + or - 
        left: If a contig is on the left or right side of a pair. 

    Returns: 
        heads or tails 
    """
    if (orientation == "+") & (left):
        return "tail"
    elif (orientation == "+") & (not left):
        return "head"
    elif (orientation == "-") & (left):
        return "head"
    elif (orientation == "-") & (not left):
        return "tail"



def update_table(contig: str):
    """
    Update the adjacency graph once we reach a contig during the DCJ algorithm.

    Args:
        contig: The contig name. 
    """

    contig_table = contig["contig_table"]
    position_table = contig["position_table"]
    contig_name = contig["contig_name"]
    scaffold_name = contig["scaffold_name"]
    length = contig["length"]
    index = contig["index"]
    orientation = contig["orientation"]

    # Update contig to index table.
    if contig_name not in contig_table:
        contig_table[contig_name] = {}
    contig_table[contig_name][orientation] = index
    contig_table[contig_name]["scaffold"] = scaffold_name
    contig_table[contig_name]["length"] = length

    # Update index to contig table.
    if index not in position_table:
        position_table[contig["index"]] = [(contig_name, orientation)]
    else:
        position_table[index].append((contig_name, orientation))


def create_table(assembly: dict) -> tuple:
    """
    Create the DCJ adjacency graph by creating two dictionaries. 

    Args:
        assembly: The dictionary created by reading the assembly AGP.

    Returns:
        Two dictionaries which together represent the adjacency matrix. One 
        dictionary indicates how contigs are connected, and the other records
        where they are in the assembly.  
    """
    
    index = 0
    contig_table = {}
    position_table = {}

    # For each scaffold in the assembly -
    for scaffold_name in assembly:
        contigs = list(assembly[scaffold_name].keys())

        # Put in first telomere.
        length = (
            assembly[scaffold_name][contigs[0]]["end"]
            - assembly[scaffold_name][contigs[0]]["start"]
            + 1
        )
        update_table(
            {
                "contig_table": contig_table,
                "position_table": position_table,
                "contig_name": contigs[0],
                "scaffold_name": scaffold_name,
                "length": length,
                "index": index,
                "orientation": get_orientation(
                    orientation=assembly[scaffold_name][contigs[0]]["orientation"],
                    left=False,
                ),
            }
        )

        index += 1

        # Put in the remaining edges.
        for i in range(len(contigs) - 1):
            length = (
                assembly[scaffold_name][contigs[i]]["end"]
                - assembly[scaffold_name][contigs[i]]["start"]
                + 1
            )
            update_table(
                {
                    "contig_table": contig_table,
                    "position_table": position_table,
                    "contig_name": contigs[i],
                    "scaffold_name": scaffold_name,
                    "length": length,
                    "index": index,
                    "orientation": get_orientation(
                        orientation=assembly[scaffold_name][contigs[i]]["orientation"],
                        left=True,
                    ),
                }
            )
            length = (
                assembly[scaffold_name][contigs[i + 1]]["end"]
                - assembly[scaffold_name][contigs[i + 1]]["start"]
                + 1
            )
            update_table(
                {
                    "contig_table": contig_table,
                    "position_table": position_table,
                    "contig_name": contigs[i + 1],
                    "scaffold_name": scaffold_name,
                    "length": length,
                    "index": index,
                    "orientation": get_orientation(
                        orientation=assembly[scaffold_name][contigs[i + 1]][
                            "orientation"
                        ],
                        left=False,
                    ),
                }
            )

            index += 1

        # Put in the second telomere.
        length = (
            assembly[scaffold_name][contigs[-1]]["end"]
            - assembly[scaffold_name][contigs[-1]]["start"]
            + 1
        )
        update_table(
            {
                "contig_table": contig_table,
                "position_table": position_table,
                "contig_name": contigs[-1],
                "scaffold_name": scaffold_name,
                "length": length,
                "index": index,
                "orientation": get_orientation(
                    orientation=assembly[scaffold_name][contigs[-1]]["orientation"],
                    left=True,
                ),
            }
        )

        index += 1

    return contig_table, position_table


def follow_path(contigs: dict, positions: dict, assembly: str, index: int, visited: list) -> tuple:
    """
    Create a single path in the adjacency graph by following a contig from the 
    assembly to the reference back and forth. We also record the total contig 
    length of a path based on the two longest contigs in the path. Lengths of 
    even paths are ignored since they do not contribute to the edit distance. 

    Args:
        contigs: A dictionary containing contig adjacency information for 
            both the assembly and the reference. 
        positions: A dictionary containing contig position information for both
            the assembly and the reference.
        assembly: First or second, this lets us keep track of which side of the 
            adjacency graph we are on. 
        index: The index of the starting contig in the positions dictionary. 
        visited: List of all visited contigs.

    Returns:
        The type of path - cycle, even, or odd is returned, as well as the
        number of contigs in the path. In addition, the list of visited nodes is
        updated. 
    """

    # Get ready to switch between first and second assemblies.
    switch = {"first": "second", "second": "first"}

    # Start with the telomere or cycle component.
    previous_edge = positions[assembly][index]
    contig_name, orientation = previous_edge[0]
    length = 0
    path = [(assembly, index)]
    follow = True
    while follow:
        if assembly == "first":
            edge_length = 0
            for contig_name, orientation in previous_edge:
                edge_length += contigs[assembly][contig_name]["length"]
            if edge_length > length:
                length = edge_length

        # Find where that edge leads.
        assembly = switch[assembly]
        index = contigs[assembly][contig_name][orientation]

        # If we have already visited this edge then return.
        if (assembly, index) in path:
            visited += path
            return "cycle", length

        # Add the edge to the path.
        path.append((assembly, index))

        # Find the other component of the edge.
        current_edge = positions[assembly][index]
        other = list(set(current_edge) - set([(contig_name, orientation)]))

        # If there is no other component, then it was a telomere and we should
        # stop the path. Otherwise, continue following it.
        if len(other) > 0:
            previous_edge = other
            contig_name, orientation = previous_edge[0]
        else:
            follow = False

    # The number of edges is one less than the number of edges in paths.
    visited += path
    if ((len(path) - 1) % 2) == 0:
        return "even", 0
    return "odd", length


def count_paths(first: tuple, second: tuple) -> tuple:
    """
    Count all types of paths following the DCJ algorithm. 

    Args:
        first: The reference genome as (contig_table, position_table). 
        second: The assembly as (contig_table, position, table).

    Returns:
        A dictionary showing the counts for each type of path, the total number
        of contigs, and the length of contigs placed correctly. 
    """

    first_contigs, first_positions = first
    second_contigs, second_positions = second
    contigs = {"first": first_contigs, "second": second_contigs}
    positions = {"first": first_positions, "second": second_positions}

    # List of unvisited nodes.
    visited = []

    # Keep track of path types.
    paths = {"odd": 0, "even": 0, "cycle": 0}

    # Keep track of total correct length.
    correct_length = 0

    # Follow the paths through each telomere on the first assembly.
    for index in first_positions:
        node = first_positions[index]
        if (len(node) == 1) and (("first", index) not in visited):
            path_type, length = follow_path(contigs, positions, "first", index, visited)
            paths[path_type] += 1
            correct_length += length

    # Follow the paths through each telomere on the second assembly.
    for index in second_positions:
        node = second_positions[index]
        if (len(node) == 1) and (("second", index) not in visited):
            path_type, length = follow_path(contigs, positions, "second", index, visited)
            paths[path_type] += 1
            correct_length += length

    # Follow the paths through unvisited nodes, which are cycles.
    for index in first_positions:
        node = first_positions[index]
        if ("first", index) not in visited:
            path_type, length = follow_path(contigs, positions, "first", index, visited)
            paths[path_type] += 1
            correct_length += length

    return paths, len(first_contigs), correct_length


def perfect_weighted_distance(reference: dict) -> int:
    """
    The nature of the algorithm causes some double counting to occur when we are
    recording contig lengths. Here we find out what the "length" of a completely
    correct assembly would be, and use that as the max length. 

    Args:
        reference: The reference AGP dictionary. 

    Returns:
        The maximum possible length from a completely correct assembly. 
    """

    contigs, positions = reference
    length = 0

    for index in positions:

        # Get the total length for each edge.
        for contig_name, orientation in positions[index]:
            length += contigs[contig_name]["length"]

    return length

def distance(reference, assembly):
    """
    From the Bergeron 2006 paper, the edit distance using the Double Cut or Join
    operation is d(a,b) = N - C - I/2 where N is the number of contigs, c is the
    number of cycles, and I is the number of odd paths.

    Args:
        reference: The AGP dictionary for the reference genome. 
        assembly: The AGP dictionary for the assembly. 

    Returns:
        A dictionary of scores. 
    """

    # Create the tables for the adjacency graph.
    reference = create_table(reference)
    assembly = create_table(assembly)

    # Make sure the reference and assembly contain the same contigs.
    if set(reference[0].keys()) != set(assembly[0].keys()):
        print("Reference and Assembly have unequal sets of contigs!")
        return

    print()
    print("   Scaffolding Statistics    ")
    print("-----------------------------")

    # Grouping distance.
    group = grouping.count(first=reference, second=assembly)
    print(f"            Grouping: {group*100:>5.2f}%")

    order = ordering.count(first=reference, second=assembly)
    print(f"            Ordering: {order*100:>5.2f}%")

    orientation = ordering.count(first=reference, second=assembly, orientation=True)
    print(f"         Orientation: {orientation*100:>5.2f}%")

    # Count the number of paths and calculate edit distance.
    paths, total, correct_length = count_paths(reference, assembly)
    edit_distance = int(total - paths["cycle"] - paths["odd"] / 2)
    accuracy = correct_length / perfect_weighted_distance(reference)

    print(f"            Accuracy: {accuracy*100:>5.2f}%")
    print(f"       Edit Distance: {edit_distance:>6}")
    print()

    distances = {
        "grouping": group,
        "order": order,
        "orientation": orientation,
        "accuracy": accuracy,
        "edit_distance": edit_distance,
    }

    return distances

def main(reference, assembly):
    reference_contigs = agp.read(reference)
    assembly_contigs = agp.read(assembly)

    distances = distance(reference_contigs, assembly_contigs)

    return distances


if __name__ == "__main__":
    arguments = setup()
    main(reference=arguments.reference, assembly=arguments.assembly)
