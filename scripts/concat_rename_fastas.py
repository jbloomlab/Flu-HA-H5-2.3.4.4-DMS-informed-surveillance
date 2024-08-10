import argparse

import Bio.SeqIO

import pandas as pd


def concat_rename_fastas(fastadir, sample_rename_tsv, fasta):
    sample_rename = pd.read_csv(
        sample_rename_tsv, sep="\t", header=None, names=["orig", "new"]
    )

    with open(fasta, "w") as f:
        for tup in sample_rename.itertuples():
            seq = Bio.SeqIO.read(f"{fastadir}/{tup.new}.fa", format="fasta")
            assert seq.id == tup.new, f"{tup.new=}\n{seq=}"
            f.write(f">{tup.orig}\n{str(seq.seq)}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Concatenate and rename FASTAs",
    )
    parser.add_argument(
        "-i",
        help="directory with input FASTAs",
        dest="fastadir",
        required=True,
    )
    parser.add_argument(
        "-s",
        help="TSV file (no header) with original and renamed names",
        dest="sample_rename",
        required=True,
    )
    parser.add_argument(
        "-o",
        help="output concatenated FASTA",
        dest="fasta",
        required=True,
    )

    args = parser.parse_args()

    concat_rename_fasta(args.fastadir, args.sample_rename, args.fasta)
