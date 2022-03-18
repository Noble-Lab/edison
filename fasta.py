def read(filename: pathlib.Path) -> dict:
    """
    Read a FASTA file and return a dictionary containing the sequence names and
    the sequences.

    Args:
        filename: The file path to the FASTA file.

    Returns:
        A dictionary containing sequence names to sequence sequences.
    """

    sequences = {}
    with open(filename) as infile:
        sequence = ""

        for line in infile:
            # If it is a header line -
            if ">" in line:

                # Add sequence to dictionary.
                if len(sequence) > 0:
                    sequences[name] = sequence
                sequence = ""

                # Trim header and newline.
                name = line[1:-1]

            else:
                sequence += line.strip()

    # Add last sequence.
    sequences[name] = sequence

    return sequences


def write(sequences: dict, filename: pathlib.Path):
    """
    Write a dictionary contain sequence names and sequences to a FASTA file.  

    Args:
        sequences: The dictionary containing sequence names and sequences. 
        filename: Output filepath.
    """
    with open(filename, "w") as outfile:
        for name in sequences:

            # Write the header.
            header = f">{name}\n"
            outfile.write(header)
            sequence = sequences[name]

            # Write 80 characters per line.
            i = 0
            while i <= (len(sequence) - 80):
                line = f"{sequence[i:i+80]}\n"
                outfile.write(line)
                i += 80

            # Write the last line.
            line = f"{sequence[i:]}\n"
            outfile.write(line)


def lengths(filename: pathlib.Path) -> dict:
    """
    Read a FASTA file and return a dictionary containing the sequence names and
    the sequence lengths. This function has the benefit of reading line by line 
    and therefore is memory efficient. 

    Args:
        filename: The file path of the FASTA file. 

    Returns:
        A dictionary containing sequence names and sequence lengths.
    """

    lengths = {}
    with open(filename) as infile:
        for line in infile:
            # Add a new sequence with each new header.
            if ">" == line[0]:
                name = line[1:-1]
                lengths[name] = 0

            # Add to the running total.
            else:
                lengths[name] += len(line.strip())
    return lengths