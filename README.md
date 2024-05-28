This repository contains the code and the data to reproduce the bioinformatics analyses from Fingerhut 2024 PLOS Genetic:

Co-transcriptional splicing facilitates transcription of gigantic genes
Jaclyn M. Fingerhut, Romain Lannes, Troy W. Whitfield, Prathapan Thiru, Yukiko M.Yamashita

To reproduce the analysis, please open the jupyter notebook run_everything.ipynb and follow the instructions.

If you want to reproduce the coverage plot you will need to align the data (see manuscript for parameters) and split the bam by strand (see the script: plitting_bam_per_strand.sh).

You need to have Python> 3.8 with Numpy and Pysam installed.

For Convenience, we provide a conda environment .yaml file that comes with all the required dependencies.

In the lab, we developed a new version backed by rust code to make a coverage plot. It is lightning fast; a plot takes less than 30 seconds and does not require pre-splitting the bam files. It is still in beta, but if you are interested in trying it, let us know.




