import argparse
import re

import pandas as pd


def sequential_names(input_file, output_file):
    records = []
    with open(input_file) as f:
        i = 0
        for line in f.readlines():
            if line and not line.isspace():
                i += 1
                name = line.strip()
                if re.search(r"\s", name):
                    raise ValueError(f"space in {name=}")
                records.append((name, f"seq_{i}"))
    pd.DataFrame(records).to_csv(output_file, sep="\t", header=False, index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Rename sequences to sequential number names"
    )
    parser.add_argument(
        "-i",
        help="input file listing input names one per line",
        dest="input_file",
        required=True,
    )
    parser.add_argument(
        "-o",
        help="output TSV file (no header) with original and renamed names",
        dest="output_file",
        required=True,
    )

    args = parser.parse_args()

    sequential_names(args.input_file, args.output_file)
