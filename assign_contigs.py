from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from matplotlib import patches
import pandas as pd
import numpy as np
import subprocess
import argparse
import pathlib

# Local imports.
import fasta
import agp


def setup():
    parser = argparse.ArgumentParser(
        description="Run nucmer to find contig chromosome assignments."
    )

    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        type=pathlib.Path,
        help="Reference genome sequence.",
    )

    parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        type=pathlib.Path,
        help="Assembly file.",
    )

    parser.add_argument(
        "-s",
        "--scaffolds",
        required=True,
        type=pathlib.Path,
        help="Scaffolder AGP file.",
    )

    parser.add_argument(
        "-t", 
        "--threads", 
        default="1", 
        help="Number of threads to use."
    )

    parser.add_argument(
        "-l",
        "--min-length",
        default=1000,
        help="Minimum alignment length to be considered. (bp)",
    )

    parser.add_argument(
        "-c",
        "--min-coverage",
        default=0.2,
        help="Minimum alignment percent to be considered. (0 - 1)",
    )

    return parser.parse_args()


def run_mummer(reference: pathlib.Path, assembly: pathlib.Path, threads: int):
    """
    MUMmer produces a coordinates file which describes the alignment of 
    sequences. We parse this file and produce a pandas dataframe. 

    Args:
        filename: The file path to the MUMmer show-coords coordinate file. 

    Returns: 
        A pandas dataframe containing all the columns of the MUMmer output. 
    """

    # Run nucmer.
    print("Running MUMmer.")
    subprocess.call(["nucmer", "-t", threads, reference, assembly])

    # Run delta-filter.
    with open("filter.delta", "w") as outfile:
        subprocess.call(["delta-filter", "-1", "out.delta"], stdout=outfile)

    # Run show-coords.
    with open("coordinates.txt", "w") as outfile:
        subprocess.call(
            ["show-coords", "-c", "-l", "-r", "-d", "filter.delta"],
            stdout=outfile,
        )

    # Read the MUMmer output. 
    alignment = read_mummer("coordinates.txt")

    # Clean up.
    pathlib.Path.unlink(pathlib.Path("out.delta"))
    pathlib.Path.unlink(pathlib.Path("filter.delta"))
    pathlib.Path.unlink(pathlib.Path("coordinates.txt"))

    return alignment


def read_mummer(filename: pathlib.Path) -> pd.DataFrame:
    """
    MUMmer produces a coordinates file which describes the alignment of 
    sequences. We parse this file and produce a pandas dataframe. 

    Args:
        filename: The file path to the MUMmer show-coords coordinate file. 

    Returns: 
        A pandas dataframe containing all the columns of the MUMmer output. 
    """

    with open(filename) as infile:
        lines = infile.readlines()

        # All the fields.
        alignment = {
            "reference_name": [],
            "query_name": [],
            "reference_start": [],
            "reference_end": [],
            "query_start": [],
            "query_end": [],
            "reference_alignment_length": [],
            "query_alignment_length": [],
            "percent_identity": [],
            "reference_length": [],
            "query_length": [],
            "reference_coverage": [],
            "query_coverage": [],
            "reference_orientation": [],
            "query_orientation": [],
        }

        for line in lines[5:]:
            split = line.split()

            # Where the fields show up in the file.
            alignment["reference_start"].append(int(split[0]))
            alignment["reference_end"].append(int(split[1]))
            alignment["query_start"].append(int(split[3]))
            alignment["query_end"].append(int(split[4]))
            alignment["reference_alignment_length"].append(int(split[6]))
            alignment["query_alignment_length"].append(int(split[7]))
            alignment["percent_identity"].append(float(split[9]))
            alignment["reference_length"].append(int(split[11]))
            alignment["query_length"].append(int(split[12]))
            alignment["reference_coverage"].append(float(split[14]))
            alignment["query_coverage"].append(float(split[15]))
            alignment["reference_orientation"].append(int(split[17]))
            alignment["query_orientation"].append(int(split[18]))
            alignment["reference_name"].append(split[19])
            alignment["query_name"].append(split[20])

    alignment = pd.DataFrame(alignment)
    alignment.set_index(["query_name", "reference_name"], inplace=True)
    alignment.sort_index(level=["query_name", "reference_name"], inplace=True)

    return alignment


def remove_alignments(group, min_length: int=1000, min_coverage: float=0.2):
    """
    Since MUMmer aligns sections of a sequence, it can produce multiple regions
    of alignment for a given sequence. We filter out low alignments using this 
    function. 

    Args:
        group: A pandas group_by group. 
        min_length: The minimum length a sequence needs to align to a reference,
            by summing all of its constituent alignments. 
        min_coverage: The minimum percent of the query sequence that needs to
            have aligned to the reference sequence. This is to prevent long 
            sequences which meet the min_length requirement from passing. 

    Returns: 
         True or false, which helps with the filtering function used by pandas. 
    """
    alignment_length = group["query_alignment_length"].sum()
    alignment_coverage = alignment_length / group["query_length"].iloc[0]

    if (alignment_coverage > min_coverage) and (alignment_length > min_length):
        return True
    else:
        return False


def get_assignments(df: pd.DataFrame, reference_lengths: dict, assembly_lengths: dict) -> tuple:
    """
    Creates an AGP compatible dictionary for how contigs ought to be scaffolded
    based on their alignments to a reference genome. 

    Args:
        df: The pandas dataframe with the alignments between the assembly and 
            reference genome. 
        reference_lengths: A dictionary containing sequence lengths for the 
            reference genome. 
        assembly_lengths: A dictionary containing sequence lengths for the 
            assembly. 

    Returns: 
         A dictionary which describes the optimal scaffolding of contigs based
         on their alignments. It is used to create an AGP file. Also returned is
         the number of reference chromosomes which were skipped due to low 
         alignments. 
    """

    # Keep a list of contigs in the assembly.
    remaining_contigs = list(assembly_lengths.keys())

    skipped = 0
    assembly = {}
    for scaffold_name in reference_lengths:

        # If reference scaffold has no alignments.
        if scaffold_name not in df.index.get_level_values("reference_name").unique():
            skipped += 1
            continue

        # Get the dataframe pertaining to a particular scaffold.
        segment = df.xs(scaffold_name, level="reference_name").sort_values(
            "mean_weighted_center"
        )[["orientation", "query_length"]]

        # Get rid of duplicate alignments.
        segment = segment.loc[~segment.index.duplicated()]

        # Add each contig.
        assembly[scaffold_name] = {}
        for contig in segment.itertuples():
            assembly[scaffold_name][contig.Index] = {
                "start": 1,
                "end": contig.query_length,
                "orientation": contig.orientation,
            }

            # Remove from remaining contigs.
            remaining_contigs.remove(contig.Index)

    # Add unscaffolded contigs.
    i = 1
    for contig in remaining_contigs:
        scaffold_name = f"unscaffolded_{i}"
        assembly[scaffold_name] = {
            contig: {
                "start": 1,
                "end": assembly_lengths[contig],
                "orientation": "+",
            }
        }
        i += 1

    return assembly, skipped

def plot(df: pd.DataFrame, reference: dict, assembly: dict, filename: pathlib.Path):
    """
    Plot how contigs are scaffolded on the reference genome compared to how they 
    are scaffolded in the assembly. 

    Args:
        df: The pandas dataframe with the alignments between the assembly and 
            reference genome. 
        reference: An AGP compatible dictionary showing ideal scaffolding. 
        assembly: An AGP compatible dictionary showing actual scaffolding.
    """

    # Assign scaffold names to the dataframe since they get lost during the 
    # MUMmer alignment. 
    for scaffold_name in assembly:
        for query_name in assembly[scaffold_name]:
            if query_name in df.index.get_level_values("query_name").unique():
                df.loc[query_name, "scaffold_name"] = scaffold_name

    # Compute the center of each scaffold so that they can be ordered by how 
    # they appear on a chromosome.  
    df["scaffold_weighted_center"] = df.groupby(by=["reference_name", "scaffold_name"])[
        "mean_weighted_center"
    ].transform(np.mean)
    df.sort_values(["reference_name", "scaffold_weighted_center"], inplace=True)

    # Create a pdf. 
    with PdfPages(filename) as pdf:

        # For each reference chromosome - 
        for reference_name in reference:

            i = 0
            ci = 0
            minimum = 0
            reference_maximum = reference[reference_name]
            maximum = reference_maximum

            # Rotate five colors. 
            light_colors = [
                "#64b0e4",
                "#ff9e49",
                "#6bd66b",
                "#e36869",
                "#a782c9",
            ]
            dark_colors = [
                "#1f77b4",
                "#ff7f0e",
                "#2ca02c",
                "#d62728",
                "#9467bd",
            ]

            # If reference scaffold has alignments.
            if reference_name in df.index.get_level_values("reference_name").unique():

                # Extract all the query sequences associated with it.
                section = df.xs(reference_name, level="reference_name")

                # Set figure size.
                num_contigs = len(section.index.get_level_values("query_name").unique())
                height = num_contigs * 0.2 + 1
                if height > 100:
                    height = 100
                plt.figure(figsize=(8, height), dpi=100)

                # For each scaffold.
                for scaffold_name, scaffold in section.groupby(
                    "scaffold_name", sort=False
                ):
                    j = 0

                    # For each contig in the scaffold.
                    for query_name in assembly[scaffold_name]:

                        # If the contig has been aligned.
                        if (
                            query_name
                            in section.index.get_level_values("query_name").unique()
                        ):
                            contig = df.xs(query_name, level="query_name")

                            # Get the start and end of the query coordinates.
                            min_query = contig[["query_start", "query_end"]].min().min()
                            max_query = contig[["query_start", "query_end"]].max().max()
                            query_length = contig["query_length"][0]

                            # Get the start and end of the reference coordinates.
                            reference_start = contig["reference_start"].min()
                            reference_end = contig["reference_end"].max()
                            reference_length = reference_end - reference_start

                            # Calculate how much of the query was unmapped.
                            preceeding_start = reference_start - min_query
                            proceeding_length = query_length - max_query
                            proceeding_end = reference_end + proceeding_length

                            # Keep track of the global minumum coordinate.
                            if preceeding_start < minimum:
                                minimum = preceeding_start

                            # Keep track of the global maximum coordinate.
                            if proceeding_end > maximum:
                                maximum = proceeding_end

                            # Unaligned segment of contig before main alignment.
                            preceeding_segment = patches.Rectangle(
                                xy=(preceeding_start, i),
                                width=min_query,
                                height=0.75,
                                facecolor=light_colors[ci % 5],
                                edgecolor="black",
                            )

                            # Aligned segment of the contig.
                            main_segment = patches.Rectangle(
                                xy=(reference_start, i),
                                width=reference_length,
                                height=0.75,
                                facecolor=dark_colors[ci % 5],
                                edgecolor="black",
                            )

                            # Unaligned segment of the contig after main alignment.
                            proceeding_segment = patches.Rectangle(
                                xy=(reference_end, i),
                                width=proceeding_length,
                                height=0.75,
                                facecolor=light_colors[ci % 5],
                                edgecolor="black",
                            )

                            # Add the rectangle patches.
                            plt.gca().add_patch(preceeding_segment)
                            plt.gca().add_patch(main_segment)
                            plt.gca().add_patch(proceeding_segment)

                            # Draw arrow for orientation.
                            arrowprops = {
                                "width": 1,
                                "shrink": 0.1,
                                "facecolor": "black",
                                "headwidth": 7,
                                "headlength": 7,
                            }
                            if contig["orientation"].iloc[0] == "+":
                                xy = (reference_end, i + 0.375)
                                xytext = (reference_start, i + 0.375)
                            else:
                                xy = (reference_start, i + 0.375)
                                xytext = (reference_end, i + 0.375)

                            plt.annotate(
                                text="",
                                xy=xy,
                                xytext=xytext,
                                arrowprops=arrowprops,
                            )

                            # Add scaffold label.
                            if j == 0:
                                x = (
                                    reference_end
                                    + proceeding_length
                                    + (reference_maximum * 0.01)
                                )
                                text = plt.text(
                                    x=x,
                                    y=i + 0.375,
                                    s=scaffold_name,
                                    va="center",
                                )

                                plt.gcf().canvas.draw()

                                text_end = (
                                    text.get_window_extent()
                                    .transformed(plt.gca().transData.inverted())
                                    .x1
                                )

                                if text_end > maximum:
                                    maximum = text_end

                            # Place the next contig below.
                            i -= 1

                            # Count contigs in scaffold.
                            j += 1

                    # Count color changes.
                    ci += 1
            else:
                plt.figure(figsize=(8, 2), dpi=100)

            # Add the reference sequence.
            main_segment = patches.Rectangle(
                xy=(0, i),
                width=reference_maximum,
                height=0.75,
                facecolor="gray",
                edgecolor="black",
            )

            plt.gca().add_patch(main_segment)
            plt.axis("off")

            # Add a 5% margin on the figure.
            plt.xlim([minimum - (maximum * 0.05), maximum + (maximum * 0.05)])
            plt.ylim([i, 2])
            plt.title(reference_name)
            plt.tight_layout()
            pdf.savefig()
            plt.close()


# Be sure to have the nucmer, delta-filter, and show-coords, in your path.
def main(
    reference,
    assembly,
    scaffold_file,
    outdir,
    threads,
    min_length,
    min_coverage,
):

    alignment = run_mummer(reference=reference, assembly=assembly, threads=threads)

    # If something goes wrong in the MUMmer process, stop the program. 
    if alignment is None:
        return

    contig_lengths = fasta.lengths(assembly)

    # Count the number of aligned contigs. 
    aligned_contigs = len(alignment.index.get_level_values("query_name").unique())

    # Filter out contigs with low alignment coverage or low alignment lengths.
    filtered = alignment.groupby(by=["query_name", "reference_name"]).filter(
        remove_alignments, min_length=min_length, min_coverage=min_coverage
    )
    filtered_contigs = len(filtered.index.get_level_values("query_name").unique())

    # Add orientation.
    filtered["signed_query_alignment_length"] = (
        filtered["query_alignment_length"] * filtered["query_orientation"]
    )
    filtered["orientation"] = np.where(
        filtered.groupby(by=["query_name", "reference_name"])[
            "signed_query_alignment_length"
        ].transform(np.mean)
        > 0,
        "+",
        "-",
    )

    # Add total alignment length for each contig and reference pair.
    filtered["total_reference_alignment_length"] = filtered.groupby(
        by=["query_name", "reference_name"]
    )["reference_alignment_length"].transform(sum)

    # Find the center of each group of alignments weighted by the length of the
    # alignment.
    filtered["weighted_center"] = filtered[["reference_start", "reference_end"]].mean(
        axis=1
    ) * (
        filtered["reference_alignment_length"]
        / filtered["total_reference_alignment_length"]
    )
    filtered["mean_weighted_center"] = filtered.groupby(
        by=["query_name", "reference_name"]
    )["weighted_center"].transform(sum)

    # Only keep the query reference pair with the largest total alignment.
    filtered = filtered.loc[
        filtered.groupby(by="query_name")["total_reference_alignment_length"].idxmax()
    ]

    # Sort by center of contig alignments.
    filtered.sort_values(["reference_name", "mean_weighted_center"], inplace=True)

    # Write the AGP file.
    print("Writing AGP.")
    reference_lengths = fasta.get_lengths(reference)
    scaffolds, skipped = get_assignments(
        df=filtered,
        reference_lengths=reference_lengths,
        assembly_lengths=contig_lengths,
    )
    stem = "_".join(assembly.stem.split("_")[:-1])
    agp_file = outdir.joinpath(f"{stem}_reference.agp")
    agp.write(scaffolds, agp_file)

    # Plot alignments.
    print("Plotting alignments.")
    scaffold_file = agp.read(scaffold_file)
    filename = outdir.joinpath(f"{stem}_alignment.pdf")
    plot(filtered, reference_lengths, scaffold_file, filename)

    # Alignment statistics.
    if len(filtered) > 0:
        minimum_alignment = (
            filtered["total_reference_alignment_length"] / filtered["query_length"]
        ).min() * 100

        average_alignment = (
            filtered.groupby(by="query_name")
            .first()["total_reference_alignment_length"]
            .sum()
            / filtered.groupby(by="query_name").first()["query_length"].sum()
            * 100
        )

    print()
    print("     Contig Statistics     ")
    print("-----------------------------")
    print(f"            Starting: {len(contig_lengths):>6}")
    print(f"  Missing Alignments: {len(contig_lengths) - aligned_contigs:>6}")
    print(f"      Low Alignments: {aligned_contigs - filtered_contigs:>6}")
    print(f"   Missing Reference: {skipped:>6}")
    print(f"   Minimum Alignment: {minimum_alignment:>5.2f}%")
    print(f"   Average Alignment: {average_alignment:>5.2f}%")

    if len(filtered) == 0:
        print("\nWarning, there are no alignments.")
        return

    if average_alignment < 50:
        print("\nWarning, the average query alignment coverage is low!")
        print("The assembly does not match the reference well.")

    return agp_file


if __name__ == "__main__":
    arguments = setup()
    main(
        reference=arguments.reference,
        assembly=arguments.assembly,
        scaffold_file=arguments.scaffolds,
        threads=arguments.threads,
        min_length=arguments.min_length,
        min_coverage=arguments.min_coverage,
        outdir=pathlib.Path.cwd(),
    )
