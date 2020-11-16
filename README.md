# Snakemake workflow: build_genome_indices

This workflow builds genome indices needed to perform 4sU analysis, which
includes assigning reads to both intronic and exonic regions.
performs analysis of RNA-Seq data from a 4sU treatment time course,
Sequences are downloaded, and STAR and Kallisto indices are built.
The transcriptome index is built with kb_python https://github.com/pachterlab/kb_python,
which extracts the introns plus 30 bp of flanking sequence on either side.

The following repository was used as an initial template:
https://github.com/snakemake-workflows/rna-seq-star-deseq2

## Running the workflow

Copy the files from .example/ to the output directory.

Then run the pipeline:

    snakemake --directory <outdir> --use-conda

## Authors

* Mary Kay Thompson (@marykthompson)
