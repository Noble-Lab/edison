import argparse
import pathlib

# Local imports.
import fasta
import agp

def setup():
    parser = argparse.ArgumentParser(
        description="Split a scaffold into its component contigs."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=pathlib.Path,
        help="Input FASTA file with scaffolded sequences.",
    )

    return parser.parse_args()

def break_scaffold(sequence: str, n: int = 10) -> list:
    """
    We want to split scaffolds whenever there is a run of Ns. To achieve this,
    we scan a sequence from left to right and keep track of the number of
    consecutive Ns. If the number ever exceeds the threshold, we save the
    sequence, and wait until the next non N base pair is seen to start saving
    the next sequence.

    Args:
        sequence: A sequence of nucleotides.
        n: The number of N's in a row we need to observe before breaking a 
            scaffold. 

    Returns: 
        A list of strings, where each string is set of nucleotides found before 
        or after a series of N's. 
    """

    contigs = []
    n_count = 0
    current = ""
    for base in sequence:

        current += base

        # If the base is an N, add to the N counter, otherwise reset it.
        if base != "N":
            n_count = 0
        else:
            n_count += 1

        # If the N count is greater than the threshold add the contig.
        if n_count >= n:
            if len(current[:-10]) > 0:
                contigs.append(current[:-10])
            current = ""

    # Add the last contig if it was not a run of Ns.
    if len(current) > 0:
        contigs.append(current)

    return contigs


def get_contigs(sequences: dict, n: int = 10) -> tuple:
    """
    Given a set of scaffolded sequences, break the sequences at runs of Ns and 
    then return a list of subsequences and the scaffolds they belonged to. 

    Args:
        sequences: A sequence of nucleotides.
        n: The number of N's in a row we need to observe before breaking a 
            scaffold. 

    Returns: 
        Two dictionaries. First, a dictionary representation of the scaffolded 
        assembly which can be later read by agp.write() to create an AGP file. 
        Second, a dictionary with contig names and sequences so that we can use
        them later to align against the reference genome. 
    """
    contigs = {}
    assembly = {}

    for scaffold in sequences:

        # If there are spaces in the seqeunce header, drop everything after the 
        # first one. 
        scaffold_name = scaffold.split(" ")[0]
        assembly[scaffold_name] = {}

        # Split the scaffold sequence at Ns.
        contig_sequences = break_scaffold(sequences[scaffold])

        # Record each contig found in the scaffold. 
        contig_index = 1
        for sequence in contig_sequences:
            contig_name = f"{scaffold_name}.{contig_index}"
            contigs[contig_name] = sequence
            assembly[scaffold_name][contig_name] = {
                "start": 1,
                "end": len(sequence),
                "orientation": "+",
            }

            # Iterate contig index.
            contig_index += 1

    return assembly, contigs

def main(infile, outdir):
    # Split scaffolds into contigs.
    print("Generating contigs from scaffolds.")
    scaffolds = fasta.read(infile)
    assembly, contigs = get_contigs(scaffolds)

    # Write contigs to FASTA.
    fasta_file = pathlib.Path(f"{infile.stem}_contigs.fasta")
    fasta.write(contigs, fasta_file)

    # Write the AGP file.
    agp_file = outdir.joinpath(f"{infile.stem}_assembly.agp")
    agp.write(assembly, agp_file)

    return fasta_file, agp_file


if __name__ == "__main__":
    arguments = setup()
    main(infile=arguments.input, outdir=pathlib.Path.cwd())
