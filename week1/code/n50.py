#!/usr/bin/env codon
import sys
from typing import List


def parse_fasta(file_path: str) -> List[str]:
    """
    Parse a FASTA file and return a list of sequences.
    """
    sequences: List[str] = []
    current_seq: str = ""

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line

    if current_seq:
        sequences.append(current_seq)

    return sequences


def calculate_n50(sequences: List[str]) -> int:
    """
    Calculate N50 from a list of sequences.
    """
    if not sequences:
        return 0

    # Get lengths of all sequences
    lengths: List[int] = [len(seq) for seq in sequences]

    # Sort in descending order
    lengths.sort(reverse=True)

    # Calculate total length
    total: int = sum(lengths)
    half: int = total // 2

    # Find N50
    running: int = 0
    for length in lengths:
        running += length
        if running >= half:
            return length

    return 0


def main():
    if len(sys.argv) != 2:
        print("Usage: codon n50.codon <contig.fasta>")
        sys.exit(1)

    fasta_file: str = sys.argv[1]
    sequences: List[str] = parse_fasta(fasta_file)
    n50: int = calculate_n50(sequences)
    print(f"N50 for {fasta_file}: {n50}")


if __name__ == "__main__":
    main()