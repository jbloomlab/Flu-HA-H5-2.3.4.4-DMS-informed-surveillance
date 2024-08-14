"""``snakemake`` file that runs analysis."""


import io

import ruamel.yaml as yaml

from snakemake.utils import min_version


min_version("8.16")


def yaml_str(d):
    """Convert a dictionary to a YAML str."""
    stream = io.StringIO()
    yaml.YAML().dump(d, stream)
    return stream.getvalue()


# global conda environment to run `snakemake`
conda:
    "envs/global.yml"


configfile: "config.yaml"


rule all:
    input:
        "results/usher_extracted_seqs/prot_seqs.fa",
        "results/dms_data/prot.fa",
        "results/dms_data/phenotypes.csv",
        "results/usher_prots_w_dms/alignment.fa",
        "results/usher_prots_w_dms/seq_phenotypes.tsv",
        expand(
            "results/dms_analyses/{analysis_filter}/analyze_dms.ipynb",
            analysis_filter=config["analysis_filters"],
        ),


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
        pyscript="scripts/extract_user_seqs.py",
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
        python {input.pyscript} \
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
        prot=config["dms"]["prot"],  # any terminal stop codon is removed
        phenotypes=config["dms"]["phenotypes"],
        escape_by_species=config["dms"]["escape_by_species"],
        site_numbering_map=config["dms"]["site_numbering_map"],
        decimal_scale=config["dms_decimal_scale"],
        site_numbering_schemes=config["dms"]["site_numbering_schemes"],
    output:
        prot="results/dms_data/prot.fa",
        phenotypes="results/dms_data/phenotypes.csv",
        site_numbering_map="results/dms_data/site_numbering_map.csv",
    conda:
        "envs/python.yml"
    log:
        "results/logs/get_dms_data.txt",
    script:
        "scripts/get_dms_data.py"


rule align_to_dms:
    """Align all the proteins to the DMS sequences."""
    input:
        dms_prot=rules.get_dms_data.output.prot,
        prots=rules.extract_usher_seqs.output.prot_seqs,
    output:
        alignment="results/usher_prots_w_dms/alignment.fa",
    conda:
        "envs/mafft.yml"
    threads: 4
    log:
        "results/align_to_dms.txt",
    shell:
        # after aligning, we use awk to remove the first sequence since it is dms_prot
        """
        mafft \
            --thread {threads} \
            --auto \
            --keeplength \
            --add {input.prots} \
            {input.dms_prot} \
            2> {log} \
            | awk 'BEGIN {{RS=">"; ORS=""}} NR>2 {{print ">" $0}}' \
            > {output.alignment} \
            2>> {log}
        """


rule assign_dms_phenos:
    """Assign deep mutational scanning phenotypes to all sequences."""
    input:
        alignment=rules.align_to_dms.output.alignment,
        dms_prot=rules.get_dms_data.output.prot,
        phenotypes=rules.get_dms_data.output.phenotypes,
        site_numbering_map=rules.get_dms_data.output.site_numbering_map,
    output:
        tsv="results/usher_prots_w_dms/seq_phenotypes.tsv",
    params:
        decimal_scale=config["dms_decimal_scale"],
        site_numbering_schemes=config["dms"]["site_numbering_schemes"],
    conda:
        "envs/python.yml"
    log:
        "results/logs/assign_dms_phenos.txt",
    script:
        "scripts/assign_dms_phenos.py"


rule analyze_dms:
    """Analyze the deep mutational scanning phenotypes."""
    input:
        tsv=rules.assign_dms_phenos.output.tsv,
        nb="notebooks/analyze_dms.ipynb",
    output:
        nb="results/dms_analyses/{analysis_filter}/analyze_dms.ipynb",
        filter_chart="results/dms_analyses/{analysis_filter}/filter_chart.html",
        pheno_scatter_chart="results/dms_analyses/{analysis_filter}/pheno_scatter_chart.html",
        samples="results/dms_analyses/{analysis_filter}/samples.txt",
    params:
        params_yaml=lambda wc, input, output: yaml_str(
            {
                "analysis_filters": config["analysis_filters"][wc.analysis_filter],
                "decimal_scale": config["dms_decimal_scale"],
                "pheno_scatter_init_phenos": config["pheno_scatter_init_phenos"],
                "site_numbering_schemes": list(config["dms"]["site_numbering_schemes"]),
                "input_tsv": input.tsv,
                "filter_chart_html": output.filter_chart,
                "pheno_scatter_chart_html": output.pheno_scatter_chart,
                "samples_txt": output.samples,
            },
        ),
    conda:
        "envs/python.yml"
    log:
        "results/logs/analyze_dms_{analysis_filter}.txt",
    shell:
        "papermill {input.nb} {output.nb} -y '{params.params_yaml}' &> {log}"
