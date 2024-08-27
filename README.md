# Deep mutational scanning phenotypes of influenza H5 clade 2.3.4.4b HAs

This repository visualizes the phenotypes of H5 HAs calculated as the sum of the effects of the constituent mutations as measured by the pseudovirus deep mutational scanning (DMS) of [Dadonaite et al (2024)](https://doi.org/10.1101/2024.05.23.595634).
This visualization workflow was written by Jesse Bloom.

Briefly, the H5 HA protein sequences are extracted from Angie Hinrich's [pre-built UShER tree](https://hgdownload.gi.ucsc.edu/hubs/GCF/000/864/105/GCF_000864105.1/UShER_NC_007362.1/).
The protein sequences are then aligned to the protein used in the DMS to identify all mutations relative to that strain, and phenotypes are calculated as the sum of the effects of the constituent mutations in the [DMS data](https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/).

The data are the provided in several visual forms.
Because the extrapolations of the DMS measurements are expected to be most accurate for sequences more closely related to the one used in the DMS, visualizations are made at several different protein identity cutoffs relative to the DMS strain.

## Visualizations and data

### Interactive scatter plots of DMS phenotypes
These plots group strains by their protein haplotype, and then allow you to plot points representing the different haplotypes.
You can mouseover points for details, adjust which DMS phenotypes are shown on the x- and y-axes, filter by divergence from the DMS protein, and scale or hide points based on how many sequences have that haplotypes.

The plots are at:

 - Strains from cattle with no more than 10% divergence from DMS strain: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/cattle_pheno_scatter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/cattle_pheno_scatter_chart.html)
 - Proteins with no more than 4% divergence from DMS strain: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.04_pheno_scatter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.04_pheno_scatter_chart.html)
 - Proteins with no more than 8% divergence from DMS strain: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.08_pheno_scatter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.08_pheno_scatter_chart.html)
 - Proteins with no more than 12% divergence from DMS strain: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.12_pheno_scatter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.12_pheno_scatter_chart.html)

### Nextstrain visualizations of the trees
These trees are the ones built by UShER on the nucleotide sequences, but points are colored by the DMS phenotypes based on the proteins.
You can use the dropdown at left to choose which phenotype to use to color the trees, and mouseover points for details:

  - Strains from cattle with no more than 10% divergence from the DMS strain: [https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/cattle](https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/cattle)
  - Strains with no more than 4% protein divergence from the DMS strain: [https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.04](https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.04)
  - Strains with no more than 8% protein divergence from the DMS strain: [https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.08](https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.08)
  - Strains with no more than 12% protein divergence from the DMS strain: [https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.12](https://nextstrain.org/community/jbloomlab/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.12)

### Plots showing filtering of all strains in UShER tree for each sequence set
These plots show how many sequences were retained/filtered for each sequence set shown above:

  - For the 4% max divergence set: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.04_filter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.04_filter_chart.html)
  - For the 8% max divergence set: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.08_filter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.08_filter_chart.html)
  - For the 12% max divergence set: [https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.12_filter_chart.html](https://jbloomlab.github.io/Flu-HA-H5-2.3.4.4-DMS-informed-surveillance/max-diff-0.12_filter_chart.html)

### Version of UShER pre-built tree used to get data
The file [results/usher_prebuilt/version_info.txt](results/usher_prebuilt/version_info.txt) indicates which version of the [pre-built UShER tree](https://hgdownload.gi.ucsc.edu/hubs/GCF/000/864/105/GCF_000864105.1/UShER_NC_007362.1/) was used.

### Numerical and sequence data
Here are numerical and sequence data tracked in this repo:

 - [results/dms_data/prot.fa](results/dms_data/prot.fa): sequence of protein used in DMS
 - [results/dms_data/phenotypes.csv](results/dms_data/phenotypes.csv): mutation effects on phenotypes as measured in DMS
 - [results/usher_prots_w_dms/seq_phenotypes.tsv](results/usher_prots_w_dms/seq_phenotypes.tsv): phenotypes assigned to sequences under additive models of the effects of constituent mutations
 - [results/haplotypes/](results/haplotypes): phenotypes of each protein haplotype along with number of sequences with that haplotype at each protein divergence cutoff

## DMS phenotypes and sequence numbering
The DMS phenotypes shown are those measured in [Dadonaite et al (2024)](https://doi.org/10.1101/2024.05.23.595634), and are escape from neutralization by ferret or mouse sera, HA stability, and a2-6-linked sialic acid usage.
Phenotypes suffixed by `_increase` are those that only include mutations with positive effects on the phenotype; this is of relevance as some mutations are measured to decrease (rather than increase) neutralization escape or stability but the biological relevance of those measurements are not always clear.

See the [interactive DMS data page](https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/) for more details on per-mutation effects.

Mutations are shown in three numbering schemes:
 - sequential numbering of the protein used in DMS
 - mature H3 (ectodomain) numbering
 - mature H5 (ectodomain) numbering of the proten used in DMS

See [https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/numbering.html](https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/numbering.html) for more details on sequence numbering.

Note also that UShER itself numbers with respect to A/Goose Guangdong/1/1996 which has an insertion in HA2 relative to the DMS strain.

## Caveats

1. The sequences used here are extracted from the UShER tree, which can interpolate some missing amino acids during phylogenetic placement. In addition, indels are not included. **Therefore, be sure to directly check the identities of any sequences of high interest.**

2. The UShER tree sequence set is not as highly QC-ed as for example the Nextstrain tree to remove erroneous or duplicated sequences, so there may be some spurious or low-quality sequences.

3. The DMS phenotype of sequences are just calculated as the additive effects of the mutations as measured in DMS. Particularly for sequences with many mutations, this additive model could break down and the phenotypes for mmore highly mutated sequences (relative to the DMS strain) may be less relevant.

## Workflow and running analysis
The analysis is run by a [snakemake](https://snakemake.readthedocs.io/) pipeline.
The pipeline is in [Snakefile](Snakefile), and the configuration is in [config.yaml](config.yaml).
The `conda` envs used by the pipeline are in [./envs/](envs).
To run the pipeline, build the `Flu-HA-H5-2.3.4.4-DMS-informed-surveillance` environment in [envs/global.yml](envs/global.yml), activate it, and run [Snakefile](Snakefile) with:

    snakemake -j 4 --software-deployment-method conda

On the Hutch cluster, you can just run the pipeline with the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash).

Code used by the pipeline is in [./notebooks/](notebooks) and [./scripts/](scripts).

Numerical and sequence results are placed in [./results/](results), with only some files tracked in the GitHub repo.

The interactive Altair plots are placed in [./docs/](docs), where they can be visualized using GitHub Pages.

The Nextstrain Auspice JSON files are placed in [./auspice/](auspice), where they can be visualized [as described here](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html).

## Re-running the pipeline as data are updated
As the [pre-built UShER tree](https://hgdownload.gi.ucsc.edu/hubs/GCF/000/864/105/GCF_000864105.1/UShER_NC_007362.1/) is updated, you will want to re-run the pipeline.
To do that, just remove the contents of the [./results/](results), [./docs/](docs), and [./auspice](auspice) subidirectories and then just re-run the pipeline.

The script [check_for_prebuilt_usher_update.py](check_for_prebuilt_usher_update.py) can be set up to run using `crontab -e` to check for such updates.
