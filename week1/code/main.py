#!/usr/bin/env codon
import sys
import os
from dbg import DBG
from utils import read_data


def main():
    if len(sys.argv) < 2:
        print("Usage: codon main.codon <data_dir>")
        sys.exit(1)

    data_dir = sys.argv[1]
    short1, short2, long1 = read_data(os.path.join('./', data_dir))

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    output_path = os.path.join('./', data_dir, 'contig.fasta')
    with open(output_path, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(f'>contig_{i}\n')
            f.write(f'{c}\n')


if __name__ == "__main__":
    main()