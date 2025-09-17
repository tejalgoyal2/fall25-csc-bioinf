#!/usr/bin/env codon
import os
from typing import List, Tuple


def read_fasta(path: str, name: str) -> List[str]:
    """
    Read a FASTA file and return a list of sequences.
    """
    data: List[str] = []
    with open(os.path.join(path, name), 'r') as f:
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    data.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line

        if current_seq:
            data.append(current_seq)

    print(name, len(data), len(data[0]) if data else 0)
    return data


def read_data(path: str) -> Tuple[List[str], List[str], List[str]]:
    """
    Read all data files for the dataset.
    """
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1