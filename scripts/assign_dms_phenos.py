import datetime
import re
import sys

import Bio.SeqIO

import pandas as pd


def decimal_year(date_string):
    if date_string == "?":
        return pd.NA
    elif re.fullmatch("\d{4}\-\d{2}\-\d{2}", date_string):
        date_format = "%Y-%m-%d"
    elif re.fullmatch("\d{4}\-\d{2}", date_string):
        date_format = "%Y-%m"
    elif re.fullmatch("\d{4}", date_string):
        date_format = "%Y"
    else:
        raise ValueError(f"unrecognized {date_string=}")
    date = datetime.datetime.strptime(date_string, date_format)

    # Get the start of the year and the start of the next year
    start_of_year = datetime.datetime(year=date.year, month=1, day=1)
    start_of_next_year = datetime.datetime(year=date.year + 1, month=1, day=1)

    # Calculate the number of days since the start of the year
    days_since_start_of_year = (date - start_of_year).days

    # Calculate the total number of days in the year
    total_days_in_year = (start_of_next_year - start_of_year).days

    # Calculate the decimal year
    decimal_year = date.year + days_since_start_of_year / total_days_in_year

    return decimal_year


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

decimal_scale = snakemake.params.decimal_scale
site_numbering_schemes = snakemake.params.site_numbering_schemes

aas = "ACDEFGHIKLMNPQRSTVWY*-"

dms_prot = Bio.SeqIO.read(snakemake.input.dms_prot, "fasta")
dms_prot_str = str(dms_prot.seq)
assert all(aa in aas for aa in dms_prot_str), dms_prot_str

# read phenotypes and map sequential mutation to phenotypes and other numbering schemes
phenotypes = pd.read_csv(snakemake.input.phenotypes)
phenos = [
    c for c in phenotypes.columns
    if (c not in ["wildtype", "mutant"]) and not c.endswith("_site")
]
assert "sequential_site" in phenotypes.columns
phenotypes["mutation_sequential"] = (
    phenotypes["wildtype"]
    + phenotypes["sequential_site"].astype(str)
    + phenotypes["mutant"]
)
mut_to_pheno = phenotypes.set_index("mutation_sequential")[phenos].to_dict()

# map sequential mutation to other numbering schemes
sequential_to_mutation = (
    pd.read_csv(snakemake.input.site_numbering_map)
    .assign(_i=lambda x: x["sequential_site"])
    .rename(
        columns={f"{v}_site": f"mutation_{k}" for (k, v) in site_numbering_schemes.items()}
    )
    .merge(pd.DataFrame({"mutant": list(aas)}), how="cross")
    .melt(
        id_vars=["wildtype", "mutant", "_i"],
        value_vars=[f"mutation_{scheme}" for scheme in site_numbering_schemes],
        var_name="numbering_scheme",
        value_name="site",
    )
    .assign(mutation=lambda x: x["wildtype"] + x["site"].astype(str) + x["mutant"])
    .pivot_table(
        index=["wildtype", "mutant", "_i"],
        values="mutation",
        columns="numbering_scheme",
        aggfunc="".join,
    )
    .reset_index()
    .assign(index=lambda x: x["mutation_sequential"])
    .set_index("index")
    [[f"mutation_{scheme}" for scheme in site_numbering_schemes]]
    .to_dict()
)

# build data frame of mutations in each sequence
records = []
for seq in Bio.SeqIO.parse(snakemake.input.alignment, "fasta"):
    assert len(seq) == len(dms_prot)
    seq_str = str(seq.seq)
    assert all(aa in aas for aa in seq_str), seq_str
    assert seq.id == seq.name == seq.description, f"{seq.id=}\n{seq.name=}\n{seq.description=}"

    date = decimal_year(seq.id.split("|")[-1])

    frac_aligned = sum(aa != "-" for aa in seq_str) / len(seq_str)
    aa_mutations = [
        f"{aa_dms}{r}{aa_seq}"
        for (r, (aa_dms, aa_seq)) in enumerate(zip(dms_prot_str, seq_str), start=1)
        if (aa_dms != aa_seq) and (aa_seq != "-")
    ]
    frac_divergence = len(aa_mutations) / len(seq_str)

    record = [seq.id, date, frac_aligned, frac_divergence, aa_mutations]
    for scheme in site_numbering_schemes:
        d = sequential_to_mutation[f"mutation_{scheme}"]
        record.append(", ".join(d[m] for m in aa_mutations))

    records.append(record)

dms_df = pd.DataFrame(
    records,
    columns=[
        "name",
        "date",
        "frac_aligned",
        "frac_divergence",
        "aa_mutations",
        *[f"mutations_{scheme}" for scheme in site_numbering_schemes],
    ],
)

# assign phenotypes to each sequence along with mutations and constituent effects
dms_phenos_df = []
for pheno in phenos:
    pheno_df = dms_df.assign(
        phenotype_name=pheno,
        phenotype_value=lambda x: x["aa_mutations"].map(
            lambda ms: sum(
                mut_to_pheno[pheno][m]
                for m in ms
                if (m in mut_to_pheno[pheno]) and pd.notnull(mut_to_pheno[pheno][m])
            )
        ),
    )
    dms_phenos_df.append(pheno_df)
    
dms_phenos_df = (
    pd.concat(dms_phenos_df, ignore_index=True)
    .drop(columns="aa_mutations")
    .pivot_table(
        index=[c for c in dms_df.columns if c != "aa_mutations"],
        columns="phenotype_name",
        values="phenotype_value",
    )
    .reset_index()
)

dms_phenos_df.to_csv(
    snakemake.output.tsv, index=False, float_format=f"%.{decimal_scale}f", sep="\t"
)
