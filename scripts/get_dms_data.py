import sys
import urllib

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

urllib.request.urlretrieve(snakemake.params.prot, snakemake.output.prot)

prot = Bio.SeqIO.read(snakemake.output.prot, "fasta")

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

phenotypes.to_csv(snakemake.output.phenotypes, index=False, float_format="%.4g")
