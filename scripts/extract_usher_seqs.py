import argparse

import Bio.Seq
import Bio.SeqIO

import bte

import pandas as pd


def _make_seq(seqlist, muts):
    seq = seqlist.copy()
    for mut in muts:
        wt = mut[0]
        assert wt in {"A", "C", "T", "G"}, mut
        i = int(mut[1: -1])
        m = mut[-1]
        assert m in {"A", "C", "T", "G"}, mut
        assert seq[i - 1] == wt, f"{seq[i - 1]=}, {i=}, {wt=}, {mut=}"
        seq[i - 1] = m
    return "".join(seq)


def extract_usher_seqs(
    tree,
    ref,
    ref_cds_coords,
    nt_seqs,
    cds_seqs,
    prot_seqs,
):
    t = bte.MATree(tree)

    records = []
    for nid in t.get_leaves_ids():
        muts = t.get_haplotype(nid)
        records.append((nid, muts))

    muts_df = pd.DataFrame(records, columns=["name", "muts"])
    assert len(muts_df) == t.count_leaves() == muts_df["name"].nunique()
    print(f"Read {len(muts_df)} sequences from {tree}")

    # get reference sequence and make sure there is single entry w no mutations
    refseq = Bio.SeqIO.read(ref, "fasta")
    muts_df_ref = muts_df[muts_df["name"].str.contains(refseq.id)]
    assert len(muts_df_ref) == 1, muts_df_ref
    assert (muts_df_ref["muts"].map(len) == 0).all(), muts_df_ref
    refseqlist = list(str(refseq.seq))
    refprotlen = (ref_cds_coords[1] - ref_cds_coords[0] + 1) / 3
    print(f"Read reference sequence of {len(refseq)} nucleotides from {ref}")
    print(f"The CDS spanning {ref_cds_coords} encodes protein of length {refprotlen}")

    # get all other sequences relative to reference, then get protein and CDS
    muts_df = muts_df.assign(
        seq=lambda x: x["muts"].map(lambda m: _make_seq(refseqlist, m)),
        cds=lambda x: x["seq"].str[ref_cds_coords[0] - 1: ref_cds_coords[1]],
        prot=lambda x: x["cds"].map(lambda s: str(Bio.Seq.Seq(s).translate(to_stop=True))),
        protlen=lambda x: x["prot"].map(len),
    )
    print("Distribution of protein lengths:")
    print(muts_df.groupby("protlen").aggregate(n=pd.NamedAgg("name", "count")))

    print(f"Writing sequences to {nt_seqs=}, {cds_seqs=}, {prot_seqs=}")
    with (
        open(nt_seqs, "w") as f_nt,
        open(cds_seqs, "w") as f_cds,
        open(prot_seqs, "w") as f_prot
    ):
        for tup in muts_df.itertuples():
            f_nt.write(f">{tup.name}\n{tup.seq}\n")
            f_cds.write(f">{tup.name}\n{tup.cds}\n")
            f_prot.write(f">{tup.name}\n{tup.prot}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="extract sequences from UShER MAT",
    )
    parser.add_argument(
        "--tree",
        help="input UShER MAT",
        dest="tree",
        required=True,
    )
    parser.add_argument(
        "--ref",
        help="FASTA file with reference",
        dest="ref",
        required=True,
    )
    parser.add_argument(
        "--ref_cds_coords",
        help="coordinates of CDS in reference",
        dest="ref_cds_coords",
        nargs=2,
        type=int,
        metavar=("START", "END"),
        required=True,
    )
    parser.add_argument(
        "--nt_seqs",
        help="output FASTA with extracted nucleotide sequences",
        dest="nt_seqs",
        required=True,
    )
    parser.add_argument(
        "--cds_seqs",
        help="output FASTA with extracted coding sequences",
        dest="cds_seqs",
        required=True,
    )
    parser.add_argument(
        "--prot_seqs",
        help="output FASTA with extracted protein sequences",
        dest="prot_seqs",
        required=True,
    )

    args = parser.parse_args()

    extract_usher_seqs(**vars(args))
