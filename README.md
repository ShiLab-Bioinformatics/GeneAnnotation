# Background

RNA sequencing is currently the method of choice for genome-wide profiling of gene expression. A popular approach to quantify expression levels of genes from
RNA-seq data is to map reads to a reference genome and then count mapped reads to each gene. Gene annotation data, which include chromosomal coordinates of ex-
ons for tens of thousands of genes, are required for this quantification process. For human and mouse genomes, there are several major sources of gene annotations that can be used for quantification, such as Ensembl, GENCODE, UCSC, and RefSeq databases. However, there is very little understanding of the effect that the choice of annotation has on the quantification of gene expression in a RNA-seq pipeline. In this analysis, we present results from our comparison of Ensembl and RefSeq human annotations on their impact on gene expression quantification using benchmark RNA-seq data generated by the SEquencing Quality Control (SEQC/MAQC III) consortium. We show that the use of RefSeq gene annotation led to better quantification accuracy, based on the correlation with ground truth such as expression data from >800 real-time PCR validated genes.

# Data
The raw data used in this analysis can be downloaded from [here](https://latrobeuni-my.sharepoint.com/:f:/g/personal/dchisanga_ltu_edu_au/EqsmFwoBT21LjnvqwmQNy44BR4wNl0KiomM1MRApyxAI6Q?e=4ivV0a). This consists of the following;

The raw data should be downloaded to the root directory of the analysis pipeline.

**RNA-seq data**\
RNA-seq data consists four reference RNA samples which have been well characterised by the SEQC/MAQC Consortium. The first two are A (Universal Reference RNA) and B (Human Brain Reference RNA) from the SEQC/MAQC Consortium. Other samples include C and D which are derivatives of A and B mixed in the ratios of 3:1 in C and 1:3 in D respectively. Each sample has four replicates and with each replicate having paired-end library reads. Each read was sequenced to a depth of 100bp. 

The files should be saved in the directory "fqs".

**RefSeq NCBI chromosome alias file**\
The annotation from NCBI does not use the UCSC chromosome name format, a chromosome alias file is thus provided for *RefSeq-NCBI release 39* annotation. 

**qRT-PCR gene expression**\
A TaqMan RT-PCR dataset, also, from the SEQC project with expression values measured for over a 1000 genes was used to validate the expression of the RNA-seq
data. The expression values were measured for both the UHRR and HBRR samples together with their respective combinations.

**Microarray BeadChip**\
Microarray data from the SEQC project with samples A to D hybridized to the Illumina Bead arrays

**GRCh38 reference genome**\
The human GRCh38.p12 version 34 reference genome available at [**GENCODE**](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz) is used. 

**Gene annotation files**\

The following gtf files should be downloaded from their respective gene annotation databases and saved in the same directory as the code;

[Ensembl release 100](http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz)

[GENCODE release 34](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz)

[RefSeq-NCBI release 39](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz)

*

# Software

The entire pipeline can be run in [**R version 4.0.2**](https://www.r-project.org/), to replicate the analysis, you will need the following [Bioconductor packages ](https://bioconductor.org)
* [Rsubread_2.2.6](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)
* [limma_3.44.3](https://bioconductor.org/packages/release/bioc/html/limma.html)
* [edgeR_3.30.3](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

Others include;



