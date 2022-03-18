# Parse an AGP file and return a dictionary of scaffolds and contigs.
def read(filename):
    assembly = {}
    with open(filename) as infile:
        for line in infile:

            # Skip commented lines.
            if line[0] == "#":
                continue

            # Parse each line.
            split = line.strip().split("\t")
            scaffold_name = split[0]
            segement_type = split[4]
            contig_name = split[5]
            orientation = split[8]

            # Skip gaps.
            if segement_type == "N":
                continue

            # Account for 1-based indexing.
            start = int(split[6])
            end = int(split[7])

            # Add the scaffold if it doesn't already exist.
            if scaffold_name not in assembly:
                assembly[scaffold_name] = {}

            # Add the contig.
            assembly[scaffold_name][contig_name] = {
                "start": start,
                "end": end,
                "orientation": orientation,
            }

    return assembly


def write(assembly, filename):
    with open(filename, "w") as outfile:
        for scaffold_name in assembly:
            contigs = assembly[scaffold_name]
            part_number = 1
            object_beg = 1
            object_end = 0
            total_gap = 0

            i = 0
            for contig_name in contigs:
                length = contigs[contig_name]["end"] - contigs[contig_name]["start"]
                object_end += length

                line = (
                    f"{scaffold_name}\t"
                    f"{object_beg}\t"
                    f"{object_end}\t"
                    f"{part_number}\t"
                    "W\t"
                    f"{contig_name}\t"
                    f"{contigs[contig_name]['start']}\t"
                    f"{contigs[contig_name]['end']}\t"
                    f"{contigs[contig_name]['orientation']}\n"
                )
                outfile.write(line)

                # Iterate object_beg.

                object_beg += length

                # As long as we are before the last contig, add gaps.
                if i < (len(contigs.keys()) - 1):
                    gap_length = 1000
                    object_end += gap_length
                    total_gap += gap_length
                    part_number += 1

                    line = (
                        f"{scaffold_name}\t"
                        f"{object_beg}\t"
                        f"{object_end}\t"
                        f"{part_number}\t"
                        "N\t"
                        f"{gap_length}\t"
                        "scaffold\t"
                        "yes\t"
                        "proximity_ligation\n"
                    )
                    outfile.write(line)

                    # Iterate object_beg.
                    object_beg += gap_length

                i += 1
                part_number += 1
