---
title: "Impact of gene annotation choice on the quantification of RNA-seq data"
author: "David Chisanga and Wei Shi"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,fig.align = "center")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#check if Rsubread is installed
if(!require(Rsubread))
  BiocManager::install("Rsubread","limma","edgeR","org.Hs.eg.db")
#Set filter to be used to filter lowly expressed genes
cpm.filter=0.5
FC_cutoff=2
FDR_cutoff=0.001
# Helper function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
```

## Background

The aim of this analysis was to compare the effect of genome annotations on gene level quantification results from RNA-seq. To do this, we performed a complete RNA-seq analysis by using two popular genome annotations, NCBI's RefSeq and the Ensembl annotations. RNA-seq data consisted of four reference RNA samples which have been well characterised by the SEQC/MAQC Consortium. The first two are A (Universal Reference RNA) and B (Human Brain Reference RNA) from the SEQC/MAQC Consortium. Other samples included C and D which are derivatives of A and B mixed in the ratios of 3:1 in C and 1:3 in D respectively. Each sample had four replicates and with each replicate having paired-end library reads. Each read was sequenced to a depth of 100bp. In this analysis, ERCC spike-ins were added to the reference before building the index needed for alignment. 

In addition, qRT-PCR as well as BeadChip array gene expression values were used as ground-truth in validating the differential expression calls from each annotation.

## Pre-processing 

It is advisable that the pre-processing which includes downloading the raw files, alignment and read-summarisation be done on 
a linux like server.

### Raw Data

Before running this pipeline, ensure that the raw data is downloaded from [here](https://latrobeuni-my.sharepoint.com/:f:/g/personal/dchisanga_ltu_edu_au/EqsmFwoBT21LjnvqwmQNy44BR4wNl0KiomM1MRApyxAI6Q?e=b1yJV7).

A target dataframe is created from the raw fastq files and is used through out the rest of the analysis.

```{r build.target,eval=F}
#Read fastq files - assumming the fastq files are in the folder fqs
#Check if target file exists
fq.reads1 = list.files(path = "fqs",
                       pattern = "*R1.fastq.gz$",
                       full.names = T)
fq.reads2 = list.files(path = "fqs",
                       pattern = "*R2.fastq.gz$",
                       full.names = T)
#Get replicate sample names
sample.names = unlist(lapply(basename(fq.reads1), function(x)
  gsub("SEQC2012-ILM-AGR-|_R1.fastq.gz", "", x)))
#Get sample group names (A,B,C,D)
sample.groups = unlist(lapply(basename(fq.reads1), function(x)
  strsplit(x, "-")[[1]][4]))
targets = data.frame(
  Reads1 = fq.reads1,
  Reads2 = fq.reads2,
  Sample = sample.names,
  Group = sample.groups,
  stringsAsFactors = FALSE
)
```

```{r}
#load pre-generated targets file
targets<-readRDS("targets.RDS")
```


### Index building and read mapping

The human GRCh38.p12 version 34 reference genome was downloaded from [**GENCODE**](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz). Unzip the file to the same directory as the code and then append the ERCC spikeins to the reference genome annotation as shown in the bash command below.

```{bash appendERCC,eval=F}
#append ERCC spikeins to the reference genome
cat ERCC/ERCC-sequences.fa >> GRCh38.primary_assembly.genome.fa
```

The reference genome together with ERCC RNA spikeins were indexed using the *subread-buildindex* function in **Rsubread**.

```{r buildIndex,eval=F}
#set path to index
index.path = "index/GRCh38.R34"
if (!dir.exists(dirname(index.path)))
  dir.create(dirname(index.path))
#Build index
buildindex(basename = index.path, reference = "GRCh38.primary_assembly.genome.fa")
```


The *align* function within **Rsubread** was used to perform the alignment. Gene annotation files from Ensembl, RefSeq-NCBI and the inbuilt RefSeq annotation in Rsubread were used as part of the alignment. 

The following gtf files should be downloaded from their respective gene annotation databases and saved in the same directory as the code;
[Ensembl release 100](ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz)
[RefSeq-NCBI release 39](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz)


The RefSeq GTF file from NCBI uses gene symbols for the gene_id attribute rather than EntrezIDs. Use the following python script swap gene symbols with EntrezIDs.

```{python,eval=F}
chroms={}
genes={}
count_dbxref=0
with open("GCF_000001405.39_GRCh38.p13_genomic_updated.gtf","w") as dbxref:
    with open("GCF_000001405.39_GRCh38.p13_genomic.gtf") as gtf:
        for row in gtf:
            if "#" in row:
                dbxref.write(row)
                continue
            row=row.strip().split("\t")
            #Check for GeneID in column 9
            col_9=row[8].split(";")
            gene_id=[x.strip() for x in col_9 if "GeneID" in x]
            if len(gene_id)!=0:
                gene_id=gene_id[0].replace("db_xref","").replace("GeneID:","")
                col_9[0]="gene_id "+gene_id
                row[8]=";".join(col_9)
            dbxref.write("\t".join(map(str,row))+"\n")
```


Prior to performing the alignment, ERCC spikeins are appended to the annotations as follows;

```{bash,eval=F}
#append ERCC spikeins to each gtf file
for f in *.gtf;
do
  cat ERCC/ERCC-sequences.gtf >> ${f}
done
```

The ERCC spike-ins are also appended to the inbuilt annotation from Rsubread and saved in SAF format.

```{r}
#Build annotation file for Rsubread inbuilt
inbuilt.annot<- getInBuiltAnnotation("hg38")
ercc <-
  read.delim(
    "ERCC/ERCC-sequences.saf",
    header = F,
    sep = "\t"
  )
colnames(ercc) <- colnames(inbuilt.annot)
inbuilt.annot <- rbind(inbuilt.annot, ercc)
write.table(
  inbuilt.annot,
  file = "Rsubread-inbuilt-hg38-with-ERCC.SAF",
  row.names = F,
  quote = F,
  sep = "\t"
)
```

The following list of annotations are used through out the analysis

```{r}
#build list of annotations
annotation.files <- list(
  Ensembl = "Homo_sapiens.GRCh38.100.gtf",
  RefSeq_NCBI = "GCF_000001405.39_GRCh38.p13_genomic_updated.gtf",
  RefSeq_Rsubread = "Rsubread-inbuilt-hg38-with-ERCC.SAF"
)
```

```{r}
#The annotation from NCBI does not use the UCSC chromosome name format, a chromosome alias file is needed 
NCBI.chrAlias<-"GCF_000001405.39_GRCh38.p13_chromAlias.csv"
```

```{r alignment,eval=F}
#run once
if (!file.exists("Alignment.summary.statistics.RData"))
{
  summary.stats <-
    sapply(names(annotation.files), function(annot) {
      bams <- file.path("bams", annot, paste0(targets$Sample, ".bam"))
      if (!dir.exists(dirname(bams[1])))
        dir.create(dirname(bams[1]), recursive = T)
      #check the annotation name to determine if an alias file is needed
      chrs = NULL
      if (annot == "RefSeq_NCBI")
        chrs = NCBI.chrAlias
      #Check if annotation extension is gtf
      is.gtf <- grepl("gtf$", annotation.files[[annot]])
      stats <- align(
        index = index.path,
        readfile1 = targets$Reads1,
        readfile2 = targets$Reads2,
        type = "rna",
        output_file = bams,
        annot.ext = annotation.files[[annot]],
        useAnnotation = T,
        isGTF = is.gtf,
        chrAliases = chrs,
        nthreads = 4
      )
      return(stats)
    }, USE.NAMES = T, simplify = F)
  #The statistics from the alignment are saved in an R object and can be referenced back
  save(summary.stats, file = "Alignment.summary.statistics.RData")
}
```


### Read summarisation

Read summarization

Gene level read counts for each of the annotations were obtained using **featureCounts** a read count summarization function within *Rsubread*.

```{r featureCounts,eval=F}
#Run once
if (!file.exists("featureCounts.RData"))
{
  #Run featureCounts
  feature.counts <-
    sapply(names(annotation.files), function(annot) {
      bams <- file.path("bams", annot, paste0(targets$Sample, ".bam"))
      #check if chromosome alias file is needed
      chrs = NULL
      if (annot == "RefSeq_NCBI")
        chrs = NCBI.chrAlias
      is.gtf <- grepl("gtf$", annotation.files[[annot]])
      counts <- featureCounts(
        files = bams,
        isPairedEnd = T,
        nthreads = 4,
        isGTFAnnotationFile = is.gtf,
        chrAliases = chrs,
        annot.ext = annotation.files[[annot]]
      )
      return(counts)
    }, USE.NAMES = T, simplify = F)
  #The counts results are saved as an R object
  save(feature.counts, file = "featureCounts.RData")
}
```

```{r}
#load feature counts data
load("featureCounts.RData")
```

## Analysis

### Comparison of genomic features in Ensembl and RefSeq 

The barplot below shows the percentage of mapped fragments for each annotation and summarised in the table below.

```{r}
load("Alignment.summary.statistics.RData")
summary.data <-
  do.call(rbind , lapply(names(summary.stats), function(n) {
    n.stats <- t(summary.stats[[n]])
    xx <- colSums(feature.counts[[n]]$counts)
    n.stats <-
      cbind(
        n.stats,
        Lib.size = xx,
        PropSummarized = 100 * as.numeric(xx) / as.numeric(n.stats[, 1])
      )
    n.stats <- reshape2::melt(n.stats)
    colnames(n.stats) <- c("Library", "Metric", "Value")
    n.stats$Annotation <- n
    return(n.stats)
  }))
```


```{r,fig.height=5,fig.width=8}
prop.counts <-
  subset(summary.data,
         subset = Metric == "PropSummarized")
prop.counts$Value<-as.numeric(prop.counts$Value)
barplot(Value~Annotation+Library,
        prop.counts,
  beside = T,
  las = 2,
  col = gg_color_hue(3),
  xpd = F,
  ylim = c(0,85),
  legend.text = T,
  ylab = "Percentage (%)",
  main="Proportion of summarised fragments",
  args.legend = list(
    xpd = T,
    y = 90,
    horiz = TRUE,
    cex = 0.6,
    bty = "n"
  ),
  cex.lab = 0.8,
  cex.axis = 0.7,
  cex.names = 0.7,
  xlab=NA,
  cex.main=1
)
```

The similarities in the complexity between Ensembl (v100) and the 2 RefSeq (P13 for the public version and P1 for the Rsubread version) annotations are compared below.
```{r}
if(!require(tidyr))
  install.packages("tidyr")
```

```{r,eval=F}
#run once
if (!file.exists("annotations.features.RDS"))
{
  #The exon and gene annot
  annot.features <- sapply(names(annotation.files)[-4], function(x) {
    x.annot <- read.delim(
      annotation.files[[x]],
      sep = "\t",
      header = F,
      stringsAsFactors = F,
      comment.char = "#"
    )
    #retain only exons and genes
    x.annot <- x.annot[x.annot$V3 %in% c("gene", "exon"), ]
    if (sum(grepl("gene_version", head(x.annot$V9, n = 10))) > 0)
    {
      x.annot <- separate(
        x.annot,
        V9,
        into = c("GeneID", "Gene.version"),
        extra = "drop",
        sep = ";"
      )
      x.annot$Gene.version <-
        trimws(gsub("gene_version\\s+", "", x.annot$Gene.version))
      x.annot$GeneID <-
        with(x.annot, paste0(GeneID, ".", Gene.version))
      x.annot$Gene.version <- NULL
    }
    else
      x.annot <-
      separate(x.annot,
               V9,
               into = "GeneID",
               extra = "drop",
               sep = ";")
    x.annot$GeneID <- gsub("gene_id\\s+", "", x.annot$GeneID)
    #exclude ERCCs
    x.annot <- x.annot[!grepl("^ERCC", x.annot$V1), ]
    #exclude duplicates
    x.annot <- x.annot[!duplicated.data.frame(x.annot), ]
    return(x.annot)
  }, USE.NAMES = T, simplify = F)
  #Save pre-processed gene exon annotations for later use
  saveRDS(annot.features,file = "annotations.features.RDS")
}
```


```{r}
annot.features <- readRDS("annotations.features.RDS")
annotations.counts <- do.call(rbind, lapply(annot.features, function(x) {
  data.frame(
    All_genes = length(unique(x$GeneID)),
    Genes_with_exons =
      length(unique(x$GeneID[x$V3 == "exon"])),
    Exons = table(x$V3)[["exon"]],
    stringsAsFactors = F
  )
}))
#Include inbuilt annotation details
RefSeq.subread_df <-
  with(
    inbuilt.annot[!grepl("^ERCC", inbuilt.annot$GeneID), ],
    data.frame(
      All_genes = length(unique(GeneID)),
      Genes_with_exons = length(unique(GeneID)),
      Exons = length(GeneID),
      stringsAsFactors = F
    )
  )
```

```{r}
annot.features <- sapply(names(annot.features), function(x) {
  c.features <-feature.counts[[x]]$annotation
  x.features <- annot.features[[x]]
  if (x == "Ensembl")
    xx <- x.features[gsub("\\.[0-9].*", "", x.features$GeneID) %in%
                c.features$GeneID, ]
  else
    xx <- x.features[x.features$GeneID %in% c.features$GeneID, ]
  #remove ERCC spike-ins
  xx<-xx[!grepl("^ERCC",xx$GeneID),]
  return(xx)
}, simplify = F, USE.NAMES = T)

annotations.counts_all <-
  rbind(annotations.counts, RefSeq.subread_df)
rownames(annotations.counts_all)[4] <- "RefSeq_Rsubread"
knitr::kable(
  t(annotations.counts_all),
  format = "html",
  format.args = list(big.mark = ","),
  align = "c",
  caption = "Comparison of features in Ensembl, RefSeq-NCBI and RefSeq-
Subread"
)
```


```{r}
library(edgeR)
#Create digital gene expression list
dge <- lapply(feature.counts, function(x) {
  x <- DGEList(counts = x$counts,
               genes = x$annotation)
  colnames(x) <- targets$Sample
  x$samples$group = targets$Group
  x$genes$Chr = unlist(lapply(strsplit(x$genes$Chr, ";"),
                              function(y)
                                paste(unique(y), collapse = "|")))
  return(x)
})
```


```{r}
library(org.Hs.eg.db,quietly = T)
#Create a mapping of RefSeq and Ensembl IDs
geneIDs_maps <- list(
  Ensembl.2RefSeq = suppressMessages(as.list(
    mapIds(
      x = org.Hs.eg.db,
      keys = unique(
        gsub("\\.[1-9].*", "", annot.features$Ensembl$GeneID)
      ),
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "list"
    )
  )),
RefSeq_NCBI.2Ensembl = suppressMessages(as.list(
    mapIds(
      x = org.Hs.eg.db,
      keys = unique(annot.features$RefSeq_NCBI$GeneID),
      column = "ENSEMBL",
      keytype = "ENTREZID",
      multiVals = "list"
    )
  )),
  RefSeq_Rsubread.2Ensembl = suppressMessages(as.list(
    mapIds(
      x = org.Hs.eg.db,
      keys = dge$RefSeq_Rsubread$genes$GeneID,
      column = "ENSEMBL",
      keytype = "ENTREZID",
      multiVals = "list"
    )
  ))
)
```


```{r}
#remove unmapped IDs in ref.2ensembl
RefSeq_Rsubread.2ensembl <-
  geneIDs_maps$RefSeq_Rsubread.2Ensembl[
    sapply(geneIDs_maps$RefSeq_Rsubread.2Ensembl, function(x) {
    keep = T
    if (length(x) == 1)
      if (is.na(x))
        keep = F
      return(keep)
  }, simplify = T, USE.NAMES = F)]

RefSeq_NCBI.2ensembl <- geneIDs_maps$RefSeq_NCBI.2Ensembl[
  sapply(geneIDs_maps$RefSeq_NCBI.2Ensembl, function(x) {
  keep = T
  if (length(x) == 1)
    if (is.na(x))
      keep = F
    return(keep)
}, simplify = T, USE.NAMES = F)]

ensembl.2RefSeq_Rsubread <-
  ensembl.2RefSeq_NCBI <- geneIDs_maps$Ensembl.2RefSeq[!is.na(geneIDs_maps$Ensembl.2RefSeq)]

#remove multi-mapping IDs in ref.2ensembl
RefSeq_Rsubread.2ensembl <-
  RefSeq_Rsubread.2ensembl[sapply(RefSeq_Rsubread.2ensembl, length) == 1]
length(RefSeq_Rsubread.2ensembl)

RefSeq_NCBI.2ensembl <-
  RefSeq_NCBI.2ensembl[sapply(RefSeq_NCBI.2ensembl, length) == 1]
length(RefSeq_NCBI.2ensembl)

#remove multi-mapping IDs in ensembl.2ref
ensembl.2RefSeq_Rsubread <-
  ensembl.2RefSeq_Rsubread[unlist(lapply(ensembl.2RefSeq_Rsubread, function(x)
    length(x) == 1))]
length(ensembl.2RefSeq_Rsubread)
ensembl.2RefSeq_NCBI <-
  ensembl.2RefSeq_NCBI[unlist(lapply(ensembl.2RefSeq_NCBI, function(x)
    length(x) == 1))]
length(ensembl.2RefSeq_NCBI)

#keep common IDs
RefSeq_Rsubread.2ensembl <-
  RefSeq_Rsubread.2ensembl[intersect(names(RefSeq_Rsubread.2ensembl),
                                     unlist(ensembl.2RefSeq_Rsubread, use.names = F))]
length(RefSeq_Rsubread.2ensembl)
ensembl.2RefSeq_Rsubread<-ensembl.2RefSeq_Rsubread[unlist(RefSeq_Rsubread.2ensembl)]
length(ensembl.2RefSeq_Rsubread)

#keep common IDs
RefSeq_NCBI.2ensembl <-
  RefSeq_NCBI.2ensembl[intersect(names(RefSeq_NCBI.2ensembl),
                      unlist(ensembl.2RefSeq_NCBI, use.names = F))]
length(RefSeq_NCBI.2ensembl)
ensembl.2RefSeq_NCBI <-
  ensembl.2RefSeq_NCBI[unlist(RefSeq_NCBI.2ensembl)]
length(ensembl.2RefSeq_NCBI)
```

### Comparison of genes common between annotations

```{r,fig.width=5,fig.height=5}
genes.ensembl <-
  unique(gsub("\\.[0-9].*", "", annot.features$Ensembl$GeneID))
genes.ensembl[genes.ensembl %in% names(ensembl.2RefSeq_Rsubread)] <-
  unlist(ensembl.2RefSeq_Rsubread[match(
    genes.ensembl, names(ensembl.2RefSeq_Rsubread), nomatch = 0)], use.names = F)
genes.ensembl[genes.ensembl %in% names(ensembl.2RefSeq_NCBI)] <-
  unlist(ensembl.2RefSeq_NCBI[match(
    genes.ensembl, names(ensembl.2RefSeq_NCBI), nomatch = 0)], use.names = F)
genes.RefSeq_Rsubread <- dge$RefSeq_Rsubread$genes$GeneID
genes.RefSeq_NCBI <- dge$RefSeq_NCBI$genes$GeneID
#Exclude ERCC spike-ins
genes.ensembl <- genes.ensembl[!grepl("^ERCC", genes.ensembl)]
genes.RefSeq_Rsubread <-
  genes.RefSeq_Rsubread[!grepl("^ERCC", genes.RefSeq_Rsubread)]
genes.RefSeq_NCBI <-
  genes.RefSeq_NCBI[!grepl("ERCC", genes.RefSeq_NCBI)]
venn::venn(
  list(
    "Ensembl" = genes.ensembl,
    "RefSeq-RSubread" = genes.RefSeq_Rsubread,
    "RefSeq-NCBI" = genes.RefSeq_NCBI
  ),
  cexsn = 0.8
)
```


### Total effective lengths

The barplot below compares the total transcriptom size for the 3 annotations.

```{r,fig.height=6,fig.width=4}
transcript.size <- do.call(rbind, lapply(dge[c("Ensembl", "RefSeq_NCBI", "RefSeq_Rsubread")],
                                         function(x)
                                           sum(x$genes$Length))) / 1000000
par(mai=c(1.8,1,0.5,0.5))
barplot(
  t(transcript.size),
  las = 2,
  cex.axis = 0.7,
  col = "lightgray",
  ylab = "Total transcript size (x10^6)",
  cex.lab = 0.9,
  cex.main=1
) 
```

### Differences in effective gene lengths

```{r,fig.height=7,fig.width=4}
par(mai=c(1.8,1,0.5,0.5))
boxplot(
  lapply(dge, function(x)
    log2(x$genes$Length)),
  main = "Comparison of effective gene lengths",
  las = 2,
  range = 0,
  ylab = "Log2 effective lengths",
  cex.main = 1,
  cex.lab = 0.8,
  cex.names = 0.5
)
```

```{r}
#Get differences in gene lengths between annotations
#Ensembl vs RefSeq-NCB
ensembl_refseq_ncbi <- geneIDs_maps$RefSeq_NCBI.2Ensembl
ensembl_refseq_ncbi <-
  ensembl_refseq_ncbi[!is.na(ensembl_refseq_ncbi) &
                        unlist(lapply(ensembl_refseq_ncbi, function(x)
                          length(x) == 1))]
ensembl_refseq_ncbi <-
  log2(dge$Ensembl$genes[unlist(ensembl_refseq_ncbi, use.names = F), "Length"]) -
  log2(dge$RefSeq_NCBI$genes[names(ensembl_refseq_ncbi), "Length"])
#Ensembl vs RefSeq-Rsubread
ensembl_refseq_Rsubread <- geneIDs_maps$RefSeq_Rsubread.2Ensembl
ensembl_refseq_Rsubread <-
  ensembl_refseq_Rsubread[!is.na(ensembl_refseq_Rsubread) &
                            unlist(lapply(ensembl_refseq_Rsubread, function(x)
                              length(x) == 1))]
ensembl_refseq_Rsubread <-
  log2(dge$Ensembl$genes[unlist(ensembl_refseq_Rsubread, use.names = F), "Length"]) -
  log2(dge$RefSeq_Rsubread$genes[names(ensembl_refseq_Rsubread), "Length"])
#RefSeq-NCB vs RefSeq-Rsubread
RefSeq.genes <- intersect(dge$RefSeq_Rsubread$genes$GeneID,
                          dge$RefSeq_NCBI$genes$GeneID)
RefSeq.genes <-
  with(dge, log2(RefSeq_NCBI$genes[RefSeq.genes, "Length"]) -
         log2(RefSeq_Rsubread$genes[RefSeq.genes, "Length"]))
```

```{r,fig.height=7,fig.width=4}
par(mai=c(2.5,1,0.5,0.5))
boxplot(
  list(
    "Ensembl vs RefSeq-NCBI" = ensembl_refseq_ncbi,
    "Ensembl vs RefSeq-Rsubread" = ensembl_refseq_Rsubread,
    "RefSeq-NCBI vs RefSeq-Rsubread" = RefSeq.genes
  ),
  las = 2,
  range = 0,
  cex.axis = 0.8,
  cex.lab = 0.7,
  cex.names=0.5,
  ylab = "Absolute difference in log2 effective gene lengths",
  main = "Differences between effective lengths",
  cex.main = 1
)
```

### Experiment design matrix

The experiment design matrix was created using the Group attribute from the targets file

```{r design_matrix}
design = model.matrix( ~ 0 + Group, targets)
colnames(design) <- gsub("Group", "", colnames(design))
t(design)
```

```{r helperFunction1}
boxplot.groups<-function(num.groups,labels,x,dist.groups=3,cex.ax=1,groups,
                         mai=c(0.8,1,0.5,1),axis.las=2,legend.show=T,
                         legend.factor=1.05,...)
{
    labels.l<-length(labels)
    box.pos<-sort(unlist(lapply(1:num.groups,function(x)
      seq(x,by=dist.groups,length.out = labels.l))))
    box.txt<-seq(mean(1:num.groups),by=dist.groups,length.out = labels.l)
    box.cols<-gg_color_hue(num.groups)#rep(rev(gray.colors(num.groups)),labels.l) #
    par(mai=mai,xpd=T)
    xx<-boxplot(x,col=box.cols,xaxt="n",xlab=NA,at=box.pos,las=1,...)
     axis(1,at=box.txt,labels = labels,tick = T,las=axis.las,cex.axis=cex.ax)
    if(legend.show)
    legend(x=max(box.pos)*legend.factor,y=max(xx$stats),legend = groups,
           fill = unique(box.cols),cex=0.8,bty = "n",title = "Key",
           title.adj = 0.2)
}
```

### Gene filtering

For each annotation,genes with read counts >= `r cpm.filter` in at least 4 libraries were ratained for further analysis.

```{r gene_filters}
annotations<-names(annotation.files)
#Check distribution of reads before filtering via a boxplot
rpkm.data<-do.call(rbind,lapply(annotations,function(x){
  y<-as.data.frame(rpkm(dge[[x]],log = T))
  y<-gather(y,key = "Samples",value="log2RPKM")
  y<-cbind(y,Annotation=x)
  return(y)
}))
rpkm.data<-rpkm.data[order(rpkm.data$Samples,decreasing = F),]
```
```{r,fig.width=10,fig.height=7}
boxplot.groups(num.groups = 3,dist.groups = 4,labels = targets$Sample,cex.main=1,
         main="Distribution of intensity range before filtering",
         x=log2RPKM~Annotation+Samples,data=rpkm.data,range=0,mai=c(1,1,0.5,1.6),
         cex.axis=0.8,groups = unique(rpkm.data$Annotation),
         ylim=c(-10,15)) 
```


### Number of filtered genes

```{r filtercpms}
dgl.filt<-sapply(names(dge), function(x){
  exp<-rowSums(cpm(dge[[x]],normalized.lib.sizes = T)>= cpm.filter) >=4 
  cat(x,"")
  print(table(exp))
  return(dge[[x]][exp,,keep.lib.sizes=FALSE])
},simplify = F,USE.NAMES = T)

count.dt=do.call(cbind,sapply(annotations,function(x) 
  c(nrow(dge[[x]]),nrow(dgl.filt[[x]])),simplify = F,USE.NAMES = T))
rownames(count.dt)=c("Before","After")
```


```{r,fig.width=6,fig.height=5}
p <- barplot(
  count.dt,
  beside = T,
  legend.text = T,
  col = gg_color_hue(2),
  args.legend = list(cex = 1, bty = "n"),
  cex.axis = 0.8,
  cex.names = 0.8,
  cex.lab = 0.9,
  cex.main = 1,
  ylab = "Gene count",
  main = "Gene count before & after filtering for low counts",
  ylim = c(0, max(count.dt) + 5000)
)
count.vals <- gather(as.data.frame(count.dt))$value
text(
  p,
  y = count.vals + 1400,
  labels = prettyNum(count.vals, big.mark = ","),
  cex = 0.8
)
```

Venn diagram of genes common 

```{r}
ensembl.genes <-
  c(ensembl.2RefSeq_Rsubread[!names(
    ensembl.2RefSeq_Rsubread) %in% ensembl.2RefSeq_NCBI],
    ensembl.2RefSeq_NCBI[!names(
      ensembl.2RefSeq_NCBI) %in% ensembl.2RefSeq_Rsubread])
ensembl.genes <-
  ensembl.genes[!grepl("^ERCC", names(ensembl.genes))]
ensembl.genes <-
  c(dgl.filt$Ensembl$genes$GeneID[
    !dgl.filt$Ensembl$genes$GeneID %in% names(ensembl.genes)],
    ensembl.genes[dgl.filt$Ensembl$genes$GeneID])
```

```{r,fig.width=8,fig.height=5}
par(omi = c(0, 0, 0.5, 0))
venn::venn(
  x = list(
    Ensembl = ensembl.genes,
    "RefSeq-NCBI" = dgl.filt$RefSeq_NCBI$genes$GeneID,
    "RefSeq-Rsubread" = dgl.filt$RefSeq_Rsubread$genes$GeneID
  ),
)
title(main = "Number of common genes after filtering for lowly expressed genes",
      outer = T,
      cex.main = 1)
```

```{r}
#Boxplot of RPKMs for common genes
com.ercc<-Reduce(intersect,lapply(annotations,function(x) 
  dgl.filt[[x]]$genes$GeneID))
com.ercc.lst=as.list(com.ercc)
names(com.ercc.lst)=com.ercc
#get common genes
dgl.com.genes=ensembl.2RefSeq_Rsubread[names(ensembl.2RefSeq_Rsubread) %in% 
              intersect(intersect(names(ensembl.2RefSeq_Rsubread),
                                  names(ensembl.2RefSeq_NCBI)),
                        dgl.filt$Ensembl$genes$GeneID)]
length(dgl.com.genes)
dgl.com.genes=c(dgl.com.genes,com.ercc.lst)
ref.common_genes<-Reduce(intersect,lapply(
  annotations[grepl("Ref",annotations)],
  function(x) dgl.filt[[x]]$genes$GeneID))
dgl.com.genes=dgl.com.genes[unlist(dgl.com.genes) %in% ref.common_genes]
cat("Number of common genes ",length(dgl.com.genes),"\n")
ensembl.r<-dgl.filt$Ensembl$genes$GeneID[!dgl.filt$Ensembl$genes$GeneID %in% 
                                                names(dgl.com.genes)]
cat("\n","Number of genes with no (unique) mapping in Ensembl:",
    length(ensembl.r),"\n")
```


### Comparison of effective gene lengths, Ensembl vs RefSeq

#### With all genes in each annotation

```{r fig.width=4,fig.height=7}
#Boxplot of all gene lengths in each annotation after filtering for lowly expressed genes
anots.gene_lengths <- sapply(annotations, function(x) {
  x.genes <- names(dgl.com.genes)
  if (grepl("Ref", x))
    x.genes <- unlist(dgl.com.genes, use.names = F)
  return(log2(dgl.filt[[x]]$genes[, "Length"]))
}, USE.NAMES = T, simplify = F)
par(mai=c(1.4,1,0.5,0.5))
boxplot(
  anots.gene_lengths,
  range = 0,
  ylab = "log2 gene-length",
  cex.main = 1,
  cex.lab = 0.9,
  cex.axis = 0.8,
  las=2,
  main = "Boxplot of effective lengths after filters"
)
```



```{r,fig.height=8,fig.width=4}
#Boxplot of Ensembl vs RefSeq_public effective lengths
#Get common genes between annotations
Ensembl.vs.RefSeq.NCBI <-
  RefSeq_NCBI.2ensembl[dgl.filt$RefSeq_NCBI$genes$GeneID]
Ensembl.vs.RefSeq.NCBI <-
  Ensembl.vs.RefSeq.NCBI[!is.na(names(Ensembl.vs.RefSeq.NCBI))]
Ensembl.vs.RefSeq.Rsubread <-
  RefSeq_Rsubread.2ensembl[dgl.filt$RefSeq_Rsubread$genes$GeneID]
Ensembl.vs.RefSeq.Rsubread <-
  Ensembl.vs.RefSeq.Rsubread[!is.na(names(Ensembl.vs.RefSeq.Rsubread))]
RefSeq.NCBI.vs.RefSeq.Rsubread <-
  intersect(rownames(dgl.filt$RefSeq_Rsubread),
            rownames(dgl.filt$RefSeq_NCBI))
par(mai=c(2.9,1,0.5,0.5))
boxplot(
  list(
    "Ensembl vs RefSeq-NCBI" = (
      log2(dgl.filt$Ensembl$genes[unlist(Ensembl.vs.RefSeq.NCBI, use.names = F), "Length"]) -
        log2(dgl.filt$RefSeq_NCBI$genes[names(Ensembl.vs.RefSeq.NCBI), "Length"])
    ),
    "Ensembl vs RefSeq-Rsubread" = (
      log2(dgl.filt$Ensembl$genes[unlist(Ensembl.vs.RefSeq.Rsubread, use.names = F), "Length"]) -
        log2(dgl.filt$RefSeq_Rsubread$genes[names(Ensembl.vs.RefSeq.Rsubread), "Length"])
    ),
    "RefSeq-NCBI vs RefSeq-Rsubread" = (
      log2(dgl.filt$RefSeq_NCBI$genes[RefSeq.NCBI.vs.RefSeq.Rsubread, "Length"]) -
        log2(dgl.filt$RefSeq_Rsubread$genes[RefSeq.NCBI.vs.RefSeq.Rsubread, "Length"])
    )
  ),
  ylab = "Difference of log2-effective gene length",
  cex.main = 0.8,
  cex.lab = 0.8,
  range = 0,
  las = 2,
  cex.names = 0.5,
  main = "Comparison of the log2 difference of lengths"
)
```


#### With common genes between annotations

```{r fig.width=5,fig.height=7}
anots.gene_lengths <- sapply(annotations, function(x) {
  x.genes <- names(dgl.com.genes)
  if (grepl("Ref", x))
    x.genes <- unlist(dgl.com.genes, use.names = F)
  return(log2(dgl.filt[[x]]$genes[x.genes, "Length"]))
}, USE.NAMES = T, simplify = F)
boxplot(
  anots.gene_lengths,
  range = 0,
  ylab = "log2 gene-length",
  cex.main = 1,
  cex.lab = 0.9,
  cex.axis = 0.8,
  main = "Boxplot of effective lengths after filters"
)
```

```{r,fig.height=8,fig.width=5}
#Boxplot of Ensembl vs RefSeq_public effective lengths
#Get common genes between annotations
par(mai=c(2.8,1,0.5,0.5))
boxplot(
  list(
    "Ensembl vs RefSeq-NCBI" = (
      log2(dgl.filt$Ensembl$genes[names(dgl.com.genes), "Length"]) -
        log2(dgl.filt$RefSeq_NCBI$genes[unlist(dgl.com.genes, use.names = F), "Length"])
    ),
    "Ensembl vs RefSeq-Rsubread" = (
      log2(dgl.filt$Ensembl$genes[names(dgl.com.genes), "Length"]) -
        log2(dgl.filt$RefSeq_Rsubread$genes[unlist(dgl.com.genes, use.names = F), "Length"])
    ),
    "RefSeq-NCBI vs RefSeq-Rsubread" = (
      log2(dgl.filt$RefSeq_NCBI$genes[unlist(dgl.com.genes), "Length"]) -
        log2(dgl.filt$RefSeq_Rsubread$genes[unlist(dgl.com.genes, use.names = F), "Length"])
    )
  ),
  ylab = "Difference of log2-effective gene length",
  cex.main = 1,
  cex.lab = 0.8,
  range = 0,
  las = 2,
  cex.names = 0.5,
  main = "Comparison of the log2 difference of lengths"
)
```


```{r,fig.width=10,fig.height=7}
#Plot boxplot of logRPKMs accross all samples
RPKM.df<-do.call(rbind,lapply(annotations,function(x){
  y<-as.data.frame(rpkm(dgl.filt[[x]],log = T))
  y<-gather(y,key = "Samples",value="log2RPKM")
  y<-cbind(y,Annotation=x)
  return(y)
}))
boxplot.groups(num.groups = 3,dist.groups = 4,labels = targets$Sample,cex.main=1,ylim=c(-10,15),
     main="Boxplot of log2RPKM values using all genes after filters",
     x=log2RPKM~Annotation+Samples,data=RPKM.df,range=0,mai = c(1,1,0.5,1.5),
     cex.ax = 0.8,cex.axis=0.8,groups =unique(RPKM.df$Annotation))
```


### Normalization

The read summarisations from each of the two annotations were transformed using voom and then nomarlised using the Library size, quantile and TMM respectively in readness for linear modelling as shown by the mean-variance trend relationships estimated from the data in the figures below


#### TMM normalization

TMM normalization is applied to account for compositional differences between libraries

```{r tmm_norm, results="asis",fig.width=10,fig.height=7,results="asis"}
dgl.filt.TMM=lapply(dgl.filt,function(x) 
  calcNormFactors(x,method = "TMM"))
tmm.rst=do.call(cbind,lapply(dgl.filt.TMM, function(x)
  x$samples$norm.factors))
rownames(tmm.rst)=targets$Sample
```


#### Voom transformation

The read summarisations from each of the two annotations were transformed using voom. Three normalisation methods; quantile, library size and TMM with library sizes were applied to the read counts for each of the annotation options.

```{r voom_setup}
#Perform voom transformation
norm.methods=c("none","quantile","TMM")
norm.methods_lst=paste(c("Library size", "Quantile","TMM"),"normalization")
names(norm.methods_lst)=norm.methods
```

```{r voom_transf,results="asis",fig.height=4,fig.width=12}
y.voom.rpkm<-sapply(norm.methods,function(m){
  sapply(annotations,function(a){
    dgl.data=dgl.filt[[a]]
    mm=m
    if(m=="TMM")
    {
      mm="none"
      dgl.data<-dgl.filt.TMM[[a]]
    }
    voom.rst<-voom(dgl.data,design = design,plot = F,normalize.method = mm)
    voom.rst$E<-voom.rst$E -log2(voom.rst$genes$Length /1000)
  },USE.NAMES = T,simplify = F)
},USE.NAMES = T,simplify = F)
names(y.voom.rpkm)<-norm.methods_lst
```

### Comparison of normalised gene expression

```{r}
fit.all=lapply(y.voom.rpkm, function(x){return(
  lapply(x,function(y) return(lmFit(y,design = design))))})
```

#### With all genes

```{r ,results="asis",fig.height=5,fig.width=12}
#Plot boxplot of fitted values for all genes
par(mfrow=c(1,3),oma=c(0,0,4,0))
for(n in names(fit.all)) {
  exp <- do.call(rbind, lapply(annotations, function(a) {
    y <- as.data.frame(fit.all[[n]][[a]]$coefficients)
    y <- gather(y, key = "Samples", value = "log2RPKM")
    y <- cbind(y, Annotation = a)
    return(y)
  }))
  boxplot.groups(
    num.groups = 3,
    dist.groups = 4,
    labels = unique(targets$Group),
    cex.main = 1,
    main = n,
    ylim = c(-10, 15),
    ylab = "log2RPKM",
    x = log2RPKM ~ Annotation + Samples,
    data = exp,
    range = 0,
    legend.factor = 1.08,
    cex.ax = 0.8,
    cex.axis = 0.8,
    groups = unique(exp$Annotation),
    mai = c(0.4, 0.8, 0.5, 1.2)
  )
}
title("RPKMs after normalization with all genes",outer = 1)
```


#### With common genes

```{r,fig.height=5,fig.width=12}
#Plot boxplot of fitted values for common genes
par(mfrow=c(1,3),oma=c(0,0,4,0))
for(n in names(fit.all)) {
  exp <- do.call(rbind, lapply(annotations, function(a) {
    genes<-unlist(dgl.com.genes)
    if(a=="Ensembl") genes<-names(dgl.com.genes)
    y <- as.data.frame(fit.all[[n]][[a]]$coefficients[genes,])
    y <- gather(y, key = "Samples", value = "log2RPKM")
    y <- cbind(y, Annotation = a)
    return(y)
  }))
  boxplot.groups(
    num.groups = 3,
    dist.groups = 4,
    labels = unique(targets$Group),
    cex.main = 1,
    main = n,
    ylim = c(-10, 15),
    ylab = "log2RPKM",
    x = log2RPKM ~ Annotation + Samples,
    data = exp,
    range = 0,
    legend.factor = 1.08,
    cex.ax = 0.8,
    cex.axis = 0.8,
    groups = unique(exp$Annotation),
    mai = c(0.4, 0.8, 0.5, 1.2)
  )
}
title("RPKMs after normalization with common genes",outer = 1)

```


### Titration monotonicty

First, Empirical Bayes moderated t-statistics were used to assess differences in expression for each gene across sample conditions. Genes were called as DE if they achieved a FDR of < `r FDR_cutoff` and a fold-change of `r FC_cutoff`..

```{r deanalysis,results="asis"}
#Create contrasts for DE analysis
contrast.design = makeContrasts(
  A.vs.B = A - B,
  A.vs.C = A - C,
  A.vs.D = A - D,
  B.vs.C = B - C,
  B.vs.D = B - D,
  C.vs.D = C - D,
  levels = colnames(design)
)
fit.treat.all <- lapply(fit.all, function(x) {
  lapply(x, function(y) {
    treat(contrasts.fit(y, contrasts = contrast.design),
          lfc = log2(FC_cutoff))
  })
})
```

```{r}
de.dt.all<-lapply(fit.treat.all,function(x)
  lapply(x,function(y) decideTests(y,p.value = FDR_cutoff)))
```


```{r helperFunTitration,echo=F}
## This function is used to determine the number of genes following the titration monotonicity
follow_mono <- function(fold_df)
{
  return(sum((
    fold_df[, "A"] > fold_df[, "B"] &
      fold_df[, "A"] >= fold_df[, "C"] & fold_df[, "C"] >=
      fold_df[, "D"] &
      fold_df[, "D"] >= fold_df[, "B"]
  ) | (
    fold_df[, "A"] < fold_df[, "B"] & fold_df[, "A"] <= fold_df[, "C"] &
      fold_df[, "C"] <= fold_df[, "D"] & fold_df[, "D"] <= fold_df[, "B"]
  )
  ))
}

#Function implements the titration monotonicity formula
titration <- function(x.vals)
{
  names(x.vals) = NULL
  x.vals = sort(x.vals, decreasing = F)
  y = log2(((3 * 2 ^ x.vals) + 1) / ((2 ^ x.vals) + 3))
  return(list(x = x.vals, y = y))
}

fit.titration <- function(x, y)
{
  l.model <- loessFit(y, x)
  j <- order(x)
  return(list(x = x[j], y = l.model$fitted[j]))
}

plot.titration_plots <-
  function(x,
           y,
           z = NULL,
           main,
           pch = 16,
           cex.axis = 1.5,
           cex.lab = 1.8,
           cex = 0.2,
           cols = c("green", "red"),
           cex.main = 2,
           tit.vals,
           cex.leg = 1.5)
  {
    plot(
      x,
      y,
      pch = 16,
      main = gsub("subread", "Rsubread", gsub("public", "NCBI", gsub("_", "-", main))),
      xlab = "A.vs.B",
      ylab = "C.vs.D",
      ylim = c(-4, 4),
      cex = cex,
      cex.axis = cex.axis,
      cex.lab = cex.lab,
      cex.main = cex.main
    )
    lines(tit.vals, col = cols[2])
    pred = fit.titration(x, y)
    yy <- y[order(x, decreasing = F)]
    mse = Metrics::mse(actual = yy, predicted = tit.vals$y)
    legend(
      "topleft",
      legend = paste0("MSE: ", round(mse, 3)),
      cex = cex.leg,
      bty = "n"
    )
    legend(
      "topright",
      legend = paste0("Gene count:", format(length(x), big.mark = ",")),
      cex = cex.leg,
      bty = "n"
    )
    pred$mse = mse
    return(pred)
  }
```

#### With Library size normalization

```{r}
titration.b4filt <- sapply(c("all", "common"), function(g) {
  tit_treat <- sapply(annotations, function(m) {
    tit <- dge[[m]]
    tit <- tit[!grepl("ERCC", rownames(tit)), , keep.lib.size = F]
    if (g == "common")
    {
      if (m == "Ensembl")
        tit <-
          tit[rownames(tit) %in% names(ensembl.2RefSeq_Rsubread), , keep.lib.sizes =
                FALSE]
      else
        tit <-
          tit[rownames(tit) %in% unlist(ensembl.2RefSeq_Rsubread, use.names = F), ,
              keep.lib.sizes = FALSE]
    }
    tt <-
      treat(contrasts.fit(lmFit(
        voom(tit, design = design, normalize.method = "none"),
        design = design
      ), contrasts = contrast.design),
      lfc = FC_cutoff)
    tit.vals = titration(tt$coefficients[, "A.vs.B"])
    return(list(
      fold.change = tt$coefficients[, c("A.vs.B", "C.vs.D")],
      titration.vals = tit.vals
    ))
  }, USE.NAMES = T, simplify = F)
}, USE.NAMES = T, simplify = F)
```

```{r}
#get titration values for annotations after filtering for lowly expressed genes
titration.afterfilter <- sapply(names(fit.treat.all), function(meth) {
  sapply(c("all", "common"), function(x) {
    sapply(annotations, function(annot) {
      n.fittreat = fit.treat.all[[meth]][[annot]]
      #n.fittreat<-n.fittreat[!grepl("^ERCC",rownames(n.fittreat)),]
      if (x == "common")
      {
        genes = names(dgl.com.genes)
        if (annot != "Ensembl")
          genes = unlist(dgl.com.genes, use.names = F)
        # print(paste0(length(genes),"-",sum(genes %in% rownames(n.fittreat))))
        n.fittreat <- n.fittreat[genes, ]
      }
      x.vals <- n.fittreat$coefficients[, "A.vs.B"]
      tit.vals = titration(x.vals)
      return(list(
        fold.change = n.fittreat$coefficients[, c("A.vs.B", "C.vs.D")],
        titration.vals = tit.vals
      ))
    }, simplify = F, USE.NAMES = T)
  }, simplify = F, USE.NAMES = T)
}, USE.NAMES = T, simplify = F)
```


```{r titration_after_filters,comparisonPlots,fig.height=12,fig.width=14}
par(
  mfrow = c(3, 3),
  oma = c(0, 3, 0, 0),
  mai = c(0.7, 0.8, 0.5, 0.2)
)
for (plt.type in names(titration.b4filt)[1])
{
  for (anot in names(titration.b4filt[[plt.type]]))
  {
    ppt <- titration.b4filt[[plt.type]][[anot]]
    plot.titration_plots(
      x = ppt$fold.change[, "A.vs.B"],
      y = ppt$fold.change[, "C.vs.D"],
      tit.vals =  ppt$titration.vals,
      main = paste0(anot, "")
    )
    if (anot == "Ensembl")
      mtext("Before filtering", 2, line = 6, cex = 1.5)
  }
}
par(mai = c(0.7, 0.8, 0.1, 0.2))
for (plt.type in names(titration.afterfilter$`Library size normalization`))
{
  for (annot in names(titration.afterfilter$`Library size normalization`[[plt.type]]))
  {
    anot.data <- titration.afterfilter$`Library size normalization`[[plt.type]][[annot]]
    plot.titration_plots(
      x = anot.data$fold.change[, "A.vs.B"],
      y = anot.data$fold.change[, "C.vs.D"],
      tit.vals =  anot.data$titration.vals,
      main = NULL
    )
    if (annot == "Ensembl")
      mtext(
        ifelse(plt.type == "all", "After filtering", "Common genes"),
        2,
        line = 6,
        cex = 1.5
      )
  }
} 
```
#### With Quantile normalization

```{r,fig.height=8,fig.width=14}
par(mfrow=c(2,3),oma=c(0,3,0,0),mai=c(0.6,0.8,0.5,0.2))
for(plt.type in names(titration.afterfilter$`Quantile normalization`))
{
  for(annot in names(titration.afterfilter$`Quantile normalization`[[plt.type]]))
  {
   show.main=NULL
   if(plt.type=="all") show.main=annot
   anot.data<-titration.afterfilter$`Quantile normalization`[[plt.type]][[annot]]
   plot.titration_plots(x = anot.data$fold.change[,"A.vs.B"],y = anot.data$fold.change[,"C.vs.D"],
                         tit.vals =  anot.data$titration.vals,
                        main=show.main) 
   if(annot=="Ensembl")
    mtext(ifelse(plt.type=="all","All genes","Common genes"),2, line=6,cex = 1.5)
  }
}
```

#### With TMM normalization

```{r,fig.height=8,fig.width=14}
par(mfrow=c(2,3),oma=c(0,3,0,0),mai=c(0.6,0.8,0.5,0.2))
for(plt.type in names(titration.afterfilter$TMM))
{
  for(annot in names(titration.afterfilter$TMM[[plt.type]]))
  {
   show.main=NULL
   if(plt.type=="all") show.main=annot
   anot.data<-titration.afterfilter$TMM[[plt.type]][[annot]]
   plot.titration_plots(x = anot.data$fold.change[,"A.vs.B"],y = anot.data$fold.change[,"C.vs.D"],
                         tit.vals =  anot.data$titration.vals,main=show.main)  
   if(annot=="Ensembl")
    mtext(ifelse(plt.type=="all","All genes","Common genes"),2, line=6,cex = 1.5)
  }
}
```

## Validation of RNA-Seq against qRT-PCR dataset

The differential expression data from the RNA-seq is validated against the qRT-PCR dataset. The raw gene expression values are first log2 transformed before being compared against the RNA-seq dataset. For genes with dudplicate qRT-PCR values, only the row with the highest mean qRT-PCR value across all samples is retained. 
```{r, results="asis", warning=F}
#Read in qRT-PCR data
rtpcr.data = read.delim(
  "Taqman-raw.txt",
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
cat("\nNumber of RT-PCR genes: ",
    dim(rtpcr.data)[1], "\n")
#Remove genes without EntrezIDs
rtpcr.data = rtpcr.data[!is.na(rtpcr.data$EntrezID), ]
cat(
  "\nNumber of RT-PCR genes after removing genes without EntrezIDs: ",
  dim(rtpcr.data)[1],
  "\n"
)
rtpcr.value.cols = colnames(rtpcr.data)[grep(pattern = "*value$",
                                             colnames(rtpcr.data))]
rtpcr.detect.cols = colnames(rtpcr.data)[grep(pattern = "*detection$",
                                              colnames(rtpcr.data))]
#Take care of duplicates by getting the genes with the highest row means
rtpcr.data = do.call(rbind, lapply(split(rtpcr.data, rtpcr.data$EntrezID), function(x) {
  if (nrow(x) > 1)
  {
    m = rowMeans(x[, rtpcr.value.cols])
    x = x[m == max(m), ]
  }
  return(x)
}))
rownames(rtpcr.data) = rtpcr.data$EntrezID
cat("\nNumber of RT-PCR genes after removing duplicates: ",
    dim(rtpcr.data)[1],
    "\n")
```

```{r}
rtpcr.data.log2=rtpcr.data
#Log2 transform the qRT-PCR values
rtpcr.data.log2[,rtpcr.value.cols]  <- log2(rtpcr.data.log2[,rtpcr.value.cols])
group.samples=c("A","B","C","D")
#Fit a linear model
rtpcr.fit <- lmFit(rtpcr.data.log2[,rtpcr.value.cols], design)
```

```{r,fig.width=12,fig.height=8}
rtpcr.com.genes<-dgl.com.genes[unlist(dgl.com.genes) %in% intersect(
  unlist(dgl.com.genes),rtpcr.data.log2$EntrezID)]
cat("Number of common genes between RT-PCR and RNA-seq:",length(rtpcr.com.genes),"\n")
```

### Correlation of groups between RNA-seq and qRT-PCR data


```{r, results="asis",fig.width=10,fig.height=6}
rtpcr.fit.data=rtpcr.fit$coefficients
rtpcr.fit.cdata<-rtpcr.fit.data[unlist(rtpcr.com.genes,use.names = F),
                               group.samples]
```



```{r, results="asis"}
rtpcr.ens.genes <-
  suppressMessages(as.list(
    mapIds(
      x = org.Hs.eg.db,
      keys = rownames(rtpcr.fit.data),
      column = "ENSEMBL",
      keytype = "ENTREZID",
      multiVals = "filter"
    )
  ))
rtpcr.ens.genes <-
  rtpcr.ens.genes[unlist(rtpcr.ens.genes) %in% dgl.filt$Ensembl$genes$GeneID]
cat("Number of common genes between PCR and Ensembl: ",
    length(rtpcr.ens.genes),
    "\n")
ref.rtpcr.genes <-
  sapply(annotations[grepl("RefSeq", annotations)], function(x)
    intersect(dgl.filt[[x]]$genes$GeneID, rownames(rtpcr.fit.data)), simplify = F, USE.NAMES = T)
pcr_ref_com <- sapply(ref.rtpcr.genes, function(x)
  length(x))
cat("Number of common genes between PCR and RefSeq")
print(pcr_ref_com)
rna_seq.rtpcr <- NULL
```

```{r}
pcr_cor.data <- list(
  all = sapply(names(fit.all), function(meth)
    do.call(rbind, lapply(annotations, function(anot)
    {
      genes.anot <- unlist(rtpcr.ens.genes, use.names = F)
      genes.rtpcr <- names(rtpcr.ens.genes)
      if (grepl("RefSeq", anot))
        genes.anot <- genes.rtpcr <- ref.rtpcr.genes[[anot]]
      anot.exp <- fit.all[[meth]][[anot]]$coefficients[genes.anot, ]
      rtpcr.exp <- rtpcr.fit.data[genes.rtpcr, group.samples]
      cor.d <- diag(cor(x = anot.exp, y = rtpcr.exp))
      return(data.frame(
        Correlation = cor.d,
        Samples = names(cor.d),
        Annotation = anot
      ))
    })), USE.NAMES = T, simplify = F),
  common = sapply(names(fit.all), function(meth)
    do.call(rbind, lapply(annotations, function(anot)
    {
      genes.anot <- names(rtpcr.com.genes)
      if (grepl("RefSeq", anot))
        genes.anot <- unlist(rtpcr.com.genes, use.names = F)
      anot.exp <- fit.all[[meth]][[anot]]$coefficients[genes.anot, ]
      cor.d <- diag(cor(x = anot.exp, y = rtpcr.fit.cdata))
      return(data.frame(
        Correlation = cor.d,
        Samples = names(cor.d),
        Annotation = anot
      ))
    })), USE.NAMES = T, simplify = F)
)
```


```{r,fig.width=7,fig.height=12}
axis.cex = 1.5
par(mfrow = c(3, 2), oma = c(0, 1, 0, 0))
for (meth in names(fit.all))
{
  show.leg = F
  com.main = NULL
  all.main = NULL
  if (meth == "Library size normalization")
  {
    show.leg = c("Ensembl", "RefSeq-NCBI", "RefSeq-Rsubread")
    com.main = "Common genes"
    all.main = "All genes"
  }
  par(mai = c(0, 1, 0.8, 0.1))
  if (meth == "TMM")
    par(mai = c(0.4, 1, 0.5, 0.1))
  barplot(
    formula = Correlation ~ Annotation + Samples,
    data = pcr_cor.data$all[[meth]],
    beside = T,
    xpd = F,
    ylim = c(0.70, 0.88),
    legend.text = NULL,
    cex.main = 2,
    col = gg_color_hue(3),
    main = all.main,
    ylab = "Correlation coefficient",
    cex.names = 1.5,
    cex.lab = 2,
    cex.axis = axis.cex,
    xlab = NA
  )
  mtext(
    paste0(meth),
    side = 2,
    line = 6,
    cex = 1.5
  )
  if (meth == "Library size normalization")
    legend(
      x = 0.2,
      y = 0.89,
      xpd = T,
      legend = c(
        paste0("Ensembl:", length(rtpcr.ens.genes)),
        paste0("RefSeq-NCBI:", pcr_ref_com[["RefSeq_NCBI"]]),
        paste0("RefSeq-Rsubread:", pcr_ref_com[["RefSeq_Rsubread"]])
      ),
      title = "Gene count",
      title.adj = 0.1,
      cex = 1.2,
      bty = "n",
      pt.cex = 1.8,
      pch = 15,
      col = gg_color_hue(3)
    )
  par(mai = c(0, 0.5, 0.8, 0.5))
  if (meth ==  "TMM normalization")
    par(mai = c(0.4, 0.5, 0.5, 0.5))
  barplot(
    formula = Correlation ~ Annotation + Samples,
    data = pcr_cor.data$common[[meth]],
    beside = T,
    xpd = F,
    ylim = c(0.70, 0.88),
    cex.main = 2,
    legend.text = show.leg,
    args.legend = list(
      xpd = T,
      x = 20,
      y = 0.89,
      bty = "n",
      cex = 1.2,
      title = "Key",
      title.adj = 0.1,
      horiz = F
    ),
    cex.axis = axis.cex,
    col = gg_color_hue(3),
    main = com.main,
    xlab = NA,
    cex.names = 1.5,
    cex.lab = 2,
    ylab = NA
  )
  if (meth == "Library size normalization")
    legend(
      x = -2,
      y = 0.89,
      xpd = T,
      legend = paste0("Gene count:", length(rtpcr.com.genes)),
      title.adj = 0.1,
      cex = 1.2,
      bty = "n",
      pt.cex = 1.8
    )
}
```

## Validation of RNA-seq against Microarray dataset 

The gene expression results from the annotations are compared against Microarray gene expression results.

```{r,readBead}
bead.data=read.table("BeadChip-Log2.table",header = T,
                     stringsAsFactors = F,sep = "\t")
#Pick max rowMean for duplicate EntrezIDs
dim(bead.data)
bead.data=do.call(rbind,lapply(split(bead.data,bead.data$EntrezID),function(x){
  m=rowMeans(x[,group.samples])
  return(x[max(m)==m,])}))
dim(bead.data)
row.names(bead.data)=bead.data[,"EntrezID"]
```

A linear model was fitted to the gene expression dataset followed by a bayesian inference of differences between groups using a trend

```{r beadDE,fig.height=5,fig.width=7}
bead.samples=colnames(bead.data)[6:ncol(bead.data)]
bead.target=data.frame(Samples=bead.samples,
                       Group=unlist(lapply(bead.samples,function(x)
                         strsplit(x,"\\.")[[1]][1])))
bead.design=model.matrix(~0+Group,bead.target)
colnames(bead.design)=unlist(lapply(colnames(bead.design), function(x) gsub("Group","",x)))
bead.fit=lmFit(object = bead.data[,bead.samples],design = bead.design)
```

### Correlations between Microarray and RNA-seq datasets
Here, the correlation of group expression values for each annotation with the microarray data are compared.

```{r corplotsMicro, results="hold",fig.width=7,fig.height=6,warning=F}
i=1
bead.rna.cor.all=NULL
ref.bead.genes<-sapply(annotations[grepl("RefSeq",annotations)],function(x)
      as.character(intersect(bead.data[,"EntrezID"],dgl.filt[[x]]$genes$GeneID)),
      simplify = F,USE.NAMES = T)
ens.bead.genes.only<-geneIDs_maps$Ensembl.2RefSeq[unlist(
  lapply(geneIDs_maps$Ensembl.2RefSeq,function(x) length(x)==1))]
ens.bead.genes.only<-ens.bead.genes.only[ens.bead.genes.only %in% 
                                           as.character(bead.data$EntrezID)&
                                           names(ens.bead.genes.only) %in%
                                           dgl.filt$Ensembl$genes$GeneID]
cat("Number of Bead genes mapping Ensembl ",length(ens.bead.genes.only),"\n")
for(n in names(ref.bead.genes))
  cat("Number of Bead genes mapping ",n,length(ref.bead.genes[[n]]),"\n")
bead.rna.cor_com=NULL
cex.p=0.4
ref_ens_bead.com_genes=dgl.com.genes[
  unlist(dgl.com.genes,use.names = F) %in% bead.data$EntrezID]
cat("Number of common genes RNA-seq and microarray: ",length(ref_ens_bead.com_genes),"\n")

array_cor.data<-list(all=sapply(names(fit.all),function(meth){
  fit=fit.all[[meth]]
  do.call(rbind,lapply(annotations,function(annot)
  {
    genes.anot<-names(ens.bead.genes.only)
    genes.mic<-unlist(ens.bead.genes.only,use.names = F)
    if(grepl("RefSeq",annot)) genes.anot<-genes.mic<-ref.bead.genes[[annot]]
    anot.exp<-fit.all[[meth]][[annot]]$coefficients[genes.anot,]
    mic.exp<-bead.fit$coefficients[genes.mic,]
    cor.d<-diag(cor(x=anot.exp,y=mic.exp))
    return(data.frame(Correlation=cor.d,Samples=names(cor.d),
                       Annotation=annot))}))
},simplify = F,USE.NAMES = T),
common=sapply(names(fit.all),function(meth){
  do.call(rbind,lapply(annotations,function(anot)
  {
    genes.anot<-names(ref_ens_bead.com_genes)
    genes.mic<-unlist(ref_ens_bead.com_genes,use.names = F)
    if(grepl("RefSeq",anot)) genes.anot<-genes.mic
    anot.exp<-fit.all[[meth]][[anot]]$coefficients[genes.anot,]
    mic.exp<-bead.fit$coefficients[genes.mic,]
    cor.d<-diag(cor(x=anot.exp,y=mic.exp))
    return(data.frame(Correlation=cor.d,Samples=names(cor.d),
                       Annotation=anot))}))},
  simplify = F,USE.NAMES = T))
```


```{r,fig.height=12,fig.width=7}
bar.ylim<-c(0.60,0.80)
pos.leg<-0.82
par(mfrow = c(3, 2), oma = c(0, 1, 0, 0))
for (meth in names(fit.all))
{
  show.leg = F
  com.main = NULL
  all.main = NULL
  if (meth == "Library size normalization")
  {
    show.leg = c("Ensembl", "RefSeq-NCBI", "RefSeq-Rsubread")
    com.main = "Common genes"
    all.main = "All genes"
  }
  par(mai = c(0, 1, 0.8, 0.1))
  if (meth == "TMM normalization")
    par(mai = c(0.4, 1, 0.5, 0.1))
  barplot(
    formula = Correlation ~ Annotation + Samples,
    data = array_cor.data$all[[meth]],
    beside = T,
    xpd = F,
    ylim = bar.ylim,
    legend.text = NULL,
    cex.main = 2,
    col = gg_color_hue(3),
    main = all.main,
    ylab = "Correlation coefficient",
    cex.names = axis.cex,
    cex.lab = 2,
    cex.axis = axis.cex,
    xlab = NA
  )
  mtext(meth,
        side = 2,
        line = 6,
        cex = 1.5)
  if (meth == "Library size normalization")
    legend(
      x = 0.2,
      y = pos.leg,
      xpd = T,
      legend = c(
        paste0("Ensembl:",
               format(
                 length(ens.bead.genes.only), big.mark = ","
               )),
        paste0("RefSeq-NCBI:", format(
          length(ref.bead.genes$RefSeq_NCBI), big.mark = ","
        )),
        paste0("RefSeq-Rsubread:", format(
          length(ref.bead.genes$RefSeq_Rsubread), big.mark = ","
        ))
      ),
      title = "Gene count",
      title.adj = 0.1,
      cex = 1.2,
      bty = "n",
      pt.cex = 1.8,
      pch = 15,
      col = gg_color_hue(3)
    )
  
  par(mai = c(0, 0.5, 0.8, 0.5))
  if (meth == "TMM normalization")
    par(mai = c(0.4, 0.5, 0.5, 0.5))
  barplot(
    formula = Correlation ~ Annotation + Samples,
    data = array_cor.data$common[[meth]],
    beside = T,
    xpd = F,
    ylim = bar.ylim,
    cex.main = 2,
    legend.text = show.leg,
    args.legend = list(
      xpd = T,
      x = 20,
      y = pos.leg,
      bty = "n",
      cex = 1.2,
      title = "Key",
      title.adj = 0.1,
      horiz = F
    ),
    cex.axis = axis.cex,
    col = gg_color_hue(3),
    main = com.main,
    xlab = NA,
    cex.names = axis.cex,
    cex.lab = 2,
    ylab = NA
  )
  if (meth == "Library size normalization")
    legend(
      x = -2,
      y = pos.leg,
      xpd = T,
      legend = paste0("Gene count:", format(
        length(ref_ens_bead.com_genes), big.mark = ","
      )),
      title.adj = 0.1,
      cex = 1.2,
      bty = "n",
      pt.cex = 1.8
    )
}
```


## Session information

```{r sessionInfo}
sessionInfo()
```
