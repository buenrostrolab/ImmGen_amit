
# Load Libraries
library(readr)
library(Matrix)
library(GenomicRanges)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(cqn)
library(stringi)
library(ggplot2)
library(reshape2)


## ATAC DATA ##

# Register Paralleization; change 2 to 1 if multiple cores are not available
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = TRUE))  # Update this with more cores if appropriate

bed <- "../immgen_dat/ImmGenATAC1219.peak.bed"   # Point to new bed file with peaks
bamdir <- "../bamFiles/amit_atac_bam/"  # Point to directory with bam files

bamfiles <- list.files(bamdir, full.names = TRUE, pattern = "\\.bam$")
peaks <- get_peaks(bed, sort_peaks = FALSE)
counts <- get_counts(bamfiles, peaks, paired = FALSE)  # Takes a while to execute
counts <- add_gc_bias(counts, genome = BSgenome.Mmusculus.UCSC.mm10)
xx <- data.matrix(assays(counts)[["counts"]])

# Sum over cell type; not replicates but multiple sequencing runs

HSC <- rowSums(xx[,c("SRR1533862.bam","SRR1533863.bam","SRR1533864.bam")])
CMP <- rowSums(xx[,c("SRR1533865.bam","SRR1533866.bam","SRR1533867.bam")])
GMP <- rowSums(xx[,c("SRR1533850.bam","SRR1533851.bam")])
B   <- rowSums(xx[,c("SRR1533847.bam","SRR1533848.bam","SRR1533849.bam")])
Mono<-         xx[,c("SRR1533852.bam")]
MEP <- rowSums(xx[,c("SRR1533868.bam", "SRR1533870.bam")]) # SRR1533869.bam is empty
EryA<-         xx[,c("SRR3001797.bam")]
NK  <- rowSums(xx[,c("SRR1533859.bam","SRR1533860.bam")])
CD8 <-         xx[,c("SRR1533861.bam")]                
CD4 <- rowSums(xx[,c("SRR1533853.bam","SRR1533854.bam", "SRR1533855.bam","SRR1533856.bam")])
Gran<- rowSums(xx[,c("SRR1533857.bam","SRR1533858.bam")])               

collapsed_counts <- data.matrix(cbind(HSC, CMP, GMP, B, Mono, MEP, EryA, NK, CD8, CD4, Gran))
write.table(collapsed_counts, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE, file = "processed/rawATAC.csv")


# CQN Normalization
#cqno <- cqn(collapsed_counts, x = counts@rowRanges@elementMetadata@listData$bias, lengthMethod = "fixed", lengths = rep(1000, dim(collapsed_counts)[1]))
#cqn_counts <- cqno$y + cqno$offset
#write.table(counts, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE, file = "../data/ImmGen_ATAC_Counts.csv")

# RNA-Seq Collapse
rna_fileList <- list.files("RNASEQ/counts/", full.names = TRUE)
rc <- sapply(rna_fileList, function(f){ read.table(f)[,2]})
colnames(rc) <- stri_sub(colnames(rc), 16, 25)

# Sum over cell type; not replicates but multiple sequencing runs
HSC <- rowSums(rc[,c("SRR1536379","SRR1536380","SRR1536381","SRR1536382")])
CMP <- rowSums(rc[,c("SRR1536389","SRR1536390","SRR1536391","SRR1536392")])
GMP <- rowSums(rc[,c("SRR1536393","SRR1536394","SRR1536395","SRR1536396")])
B   <- rowSums(rc[,c("SRR1536411","SRR1536412")])
Mono<- rowSums(rc[,c("SRR1536407","SRR1536408","SRR1536409","SRR1536410")])
MEP <- rowSums(rc[,c("SRR1536423","SRR1536424","SRR1536425","SRR1536426")])
EryA<- rowSums(rc[,c("SRR1536427","SRR1536428")])
NK  <- rowSums(rc[,c("SRR1536421","SRR1536422")])
CD8 <- rowSums(rc[,c("SRR1536417","SRR1536418","SRR1536419","SRR1536420")])
CD4 <-rowSums(rc[,c("SRR1536413","SRR1536414","SRR1536415","SRR1536416")])
Gran<- rowSums(rc[,c("SRR1536401","SRR1536402","SRR1536403","SRR1536404","SRR1536405","SRR1536406")])

# Read in one file to get first column
genenames <- read.table(rna_fileList[1], stringsAsFactors = FALSE)[,1]
collapsed_rna_counts <- data.matrix(cbind(HSC, CMP, GMP, B, Mono, MEP, EryA, NK, CD8, CD4, Gran))
rownames(collapsed_rna_counts) <- genenames
collapsed_rna_counts <- head(collapsed_rna_counts, -5)

pdf <- data.frame(colSums(collapsed_counts), colSums(collapsed_rna_counts))
colnames(pdf) <- c("ATAC_Counts", "RNA_Counts")
ggplot(pdf, aes(ATAC_Counts, RNA_Counts, label = rownames(pdf))) + geom_text() + theme_bw()+
  labs(title = "Amit Samples Counts")
write.table(collapsed_rna_counts, sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE, file = "processed/rawRNA.csv")


