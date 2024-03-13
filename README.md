# Annotating Human-Primate Ancestral Alleles with T2T Genomes

The purpose of this project is to annotate the ancestral alleles of the human-primate ancestor. We used the parsimony-like method used by the 1000 Genomes Project and Ensembl. See [section 8.3](https://www.nature.com/articles/nature15393) of the supplementary material.

Instead of using a multiple-genome alignment produced by the Enredo-Pecan-Ensembl pipeline, we use Progressive Cactus, which also infers the ancestor's genomes. 

Running `snakemake --cores 22 --configfile configs/8_way_single_copy_hg38.yaml -s calc_ancestor_allele.smk` will create a fasta of each autosome, containing the ancestral allele based in hg38.
