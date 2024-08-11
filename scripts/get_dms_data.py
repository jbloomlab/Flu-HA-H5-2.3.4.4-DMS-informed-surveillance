import io
import sys
import urllib

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# get DMS protein, removing trailing stop codon if there is one
with urllib.request.urlopen(snakemake.params.prot) as response:
    fasta_content = response.read().decode("utf-8")
fasta_io = io.StringIO(fasta_content)
prot = Bio.SeqIO.read(fasta_io, "fasta")
if prot[-1] == "*":
    prot = prot[: -1]
Bio.SeqIO.write([prot], snakemake.output.prot, "fasta")

phenotypes = (
    pd.read_csv(snakemake.params.phenotypes)
    .merge(
        (
            pd.read_csv(snakemake.params.escape_by_species)
            .assign(serum=lambda x: x["antibody"] + "_sera_escape")
            .pivot_table(
                index=["site", "wildtype", "mutant"],
                values="escape",
                columns="serum",
            )
            .reset_index()
        ),
        on=["site", "wildtype", "mutant"],
        validate="one_to_one",
        how="outer",
    )
    .rename(
        columns={
            "site": "mature_H3_site",
            "SA26 usage increase": "SA26_usage_increase",
            "entry in 293T cells": "entry_in_293T_cells"
        }
    )
    [
        [
            "sequential_site",
            "mature_H3_site",
            "mature_H5_site",
            "wildtype",
            "mutant",
            "entry_in_293T_cells",
            "stability",
            "SA26_usage_increase",
            "ferret_sera_escape",
            "mouse_sera_escape",
        ]
    ]
    .assign(
        stability_increase=lambda x: x["stability"].clip(lower=0),
        ferret_sera_escape_increase=lambda x: x["ferret_sera_escape"].clip(lower=0),
        mouse_sera_escape_increase=lambda x: x["mouse_sera_escape"].clip(lower=0),
    )
    .sort_values(["sequential_site", "mutant"])
)

assert len(phenotypes) == len(phenotypes.groupby(["sequential_site", "mutant"]))

# make sure wildtypes in protein match those in phenotypes
wts = phenotypes.set_index("sequential_site")["wildtype"].to_dict()
wts_prot = dict(enumerate(str(prot.seq), start=1))
assert all(aa == wts_prot[i] for (i, aa) in wts.items())

phenotypes.to_csv(
    snakemake.output.phenotypes,
    index=False,
    float_format=f"%.{snakemake.params.decimal_scale}f",
)
