import re
import sys

import Bio.SeqIO

import pandas as pd


#sys.stderr = sys.stdout = open(snakemake.log[0], "w")

dms_prot = Bio.SeqIO.read(snakemake.input.dms_prot, "fasta")
dms_prot_str = str(dms_prot.seq)
assert re.fullmatch("[A-Z]+", dms_prot_str), dms_prot_str

records = []
for seq in Bio.SeqIO.parse(snakemake.input.alignment, "fasta"):
    assert len(seq) == len(dms_prot)
    seq_str = str(seq.seq)
    assert re.fullmatch(r"[A-Z\-]+", seq_str)
    assert seq.id == seq.name == seq.description, f"{seq.id=}\n{seq.name=}\n{seq.description=}"

    frac_aligned = sum(aa != "-" for aa in seq_str) / len(seq_str)
    aa_mutations = [
        f"{aa_dms}{r}{aa_seq}"
        for (r, (aa_dms, aa_seq)) in enumerate(zip(dms_prot_str, seq_str), start=1)
        if (aa_dms != aa_seq) and (aa_seq != "-")
    ]
    frac_divergence = len(aa_mutations) / len(seq_str)

    records.append((seq.id, aa_mutations, frac_aligned, frac_ident))

dms_df = pd.DataFrame(
    records,
    columns=["name", "aa_mutations", "frac_aligned", "frac_divergence"],
)

print(dms_df)
