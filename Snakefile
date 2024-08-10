"""``snakemake`` file that runs analysis."""


from snakemake.utils import min_version


min_version("8.16")


# global conda environment to run `snakemake`
conda:
    "envs/global.yml"


configfile: "config.yaml"


rule all:
    input:
        "results/usher_extracted_seqs/prot_seqs.fa",
        "results/dms_data/prot.fa",
        "results/dms_data/phenotypes.csv",


rule usher_prebuilt:
    """Get pre-built UShER tree files."""
    params:
        tree=config["usher_prebuilt"]["tree"],
        version_info=config["usher_prebuilt"]["version_info"],
        ref=config["usher_prebuilt"]["ref"],
    output:
        tree="results/usher_prebuilt/tree.pb.gz",
        version_info="results/usher_prebuilt/version_info.txt",
        ref="results/usher_prebuilt/ref.fa",
    conda:
        "envs/entrez-direct.yml"
    log:
        "results/logs/usher_prebuilt.txt",
    shell:
        """
        curl -o {output.tree} {params.tree} &> {log}
        curl -o {output.version_info} {params.version_info} &>> {log}
        efetch -db nucleotide -id {params.ref} -format fasta > {output.ref} 2>> {log}
        """


rule extract_usher_seqs:
    """Extract sequences from the UShER tree."""
    input:
        tree=rules.usher_prebuilt.output.tree,
        ref=rules.usher_prebuilt.output.ref,
    output:
        nt_seqs="results/usher_extracted_seqs/nt_seqs.fa",
        cds_seqs="results/usher_extracted_seqs/cds_seqs.fa",
        prot_seqs="results/usher_extracted_seqs/prot_seqs.fa",
    params:
        ref_cds_coords=config["usher_prebuilt"]["ref_cds_coords"],
    conda:
        "envs/bte.yml"
    log:
        "results/logs/extract_usher_seqs.txt",
    shell:
        """
        python scripts/extract_usher_seqs.py \
            --tree {input.tree} \
            --ref {input.ref} \
            --ref_cds_coords {params.ref_cds_coords} \
            --nt_seqs {output.nt_seqs} \
            --cds_seqs {output.cds_seqs} \
            --prot_seqs {output.prot_seqs} \
            &> {log}
        """


rule get_dms_data:
    """Get the DMS data."""
    params:
        prot=config["dms"]["prot"],
        phenotypes=config["dms"]["phenotypes"],
        escape_by_species=config["dms"]["escape_by_species"],
    output:
        prot="results/dms_data/prot.fa",
        phenotypes="results/dms_data/phenotypes.csv",
    conda:
        "envs/python.yml",
    log:
        "results/logs/get_dms_data.txt",
    script:
        "scripts/get_dms_data.py"
