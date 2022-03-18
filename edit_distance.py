from ast import boolop
import subprocess
import argparse
import pathlib
import json

# Local imports.
import scaffold_distance
import break_scaffolds
import assign_contigs


def setup():
    parser = argparse.ArgumentParser(
        description="Determine the accuracy of a scaffolding given a reference."
    )

    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        type=pathlib.Path,
        help="Reference AGP or FASTA file.",
    )

    parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        type=pathlib.Path,
        help="Assembly AGP or FASTA file.",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help="Output directory.",
    )

    parser.add_argument(
        "-t", "--threads", default="1", help="Number of threads to use."
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


def check_files(reference: pathlib.Path, assembly: pathlib.Path):
    """
    If two fasta files are given, preprocess.
    If two AGP files are given, don't preprocess.
    If some other files are given, end.

    Args: 
        reference: File path to reference genome. 
        assembly: File path to assembly genome. 
    """

    # Check file suffixes.
    if reference.suffix not in [".agp", ".fasta", ".fa"]:
        print("Reference in incorrect file format.")
        return
    if assembly.suffix not in [".agp", ".fasta", ".fa"]:
        print("Assembly in incorrect file format.")

    # Check reference suffix.
    if reference.suffix == ".agp":
        reference_agp = True
    else:
        reference_agp = False

    # Check assembly suffix.
    if assembly.suffix == ".agp":
        assembly_agp = True
    else:
        assembly_agp = False

    # Check they are both the same suffix.
    if reference_agp and assembly_agp:
        return False
    elif (not reference_agp) and (not assembly_agp):
        return True
    else:
        print("Reference and assembly need to be in the same file format.")
        return

def check_agp(assembly: pathlib.Path):
    """
    Check if the AGP files already exist. 

    Args: 
        assembly: File path to assembly genome. 
    """

    reference = pathlib.Path(f"{assembly.stem}_reference.agp")
    assembly = pathlib.Path(f"{assembly.stem}_assembly.agp")

    if reference.is_file() and assembly.is_file():
        return reference, assembly


def check_mummer() :
    """
    Check for mummer in the environment.
    """
    code = subprocess.call(
        ["which", "mummer"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if code:
        print("Mummer is not in the environment.")
    else:
        return True


def main():
    arguments = setup()

    # Check to make sure we are given the right files.
    preprocess = check_files(reference=arguments.reference, assembly=arguments.assembly)
    if preprocess is None:
        return
    elif preprocess is False:
        reference_agp = arguments.reference
        assembly_agp = arguments.assembly

    # Check if the AGP files already exist if given FASTA files.
    if check_agp(arguments.assembly):
        reference_agp, assembly_agp = check_agp(assembly=arguments.assembly)
        preprocess = False

    if preprocess:
        if not check_mummer():
            return

        # Break the assembly scaffolds into contigs.
        assembly_fasta, assembly_agp = break_scaffolds.main(
            infile=arguments.assembly, outdir=arguments.outdir
        )

        # Find the true order of the assembly contigs.
        reference_agp = assign_contigs.main(
            reference=arguments.reference,
            assembly=assembly_fasta,
            scaffold_file=assembly_agp,
            outdir=arguments.outdir,
            threads=arguments.threads,
            min_length=arguments.min_length,
            min_coverage=arguments.min_coverage,
        )

        # If assign_contigs fails.
        if reference_agp is None:
            return

        # Remove the contigs FASTA file.
        pathlib.Path.unlink(assembly_fasta)

    # Compare the two AGP files.
    distances = scaffold_distance.main(reference=reference_agp, assembly=assembly_agp)

    filename = arguments.outdir.joinpath(f"{arguments.assembly.stem}.json")
    with open(filename, "w") as outfile:
        json.dump(distances, outfile, indent=4)


if __name__ == "__main__":
    main()
