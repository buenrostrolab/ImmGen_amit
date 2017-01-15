#!/bin/bash

# Original Code by Hide
# Modified by Caleb for MGH Erisone cluster
# And for single-end reads

module load hisat
module load python/2.7.3
module load htseq/0.6.1

fqs="SRR1536417_1.fastq.gz"

for fq in $fqs
do

echo "mapping by HISAT2"
hisat2 -q --transcriptome-mapping-only --no-discordant -x /data/aryee/caleb/hemeATAC/amit/rnaseq/HISAT_2_mm10_ImmGen/mm10.refGene -U $1 -S $1.sam 

echo "all reads in : " $1.sam
samtools view -Sc $1.sam

echo "removing unmapped but keep low-MAPQ read "
samtools view -Sb -F 4 $1.sam > $1.MAPQ0.bam	# keep only mapped reads
echo "all.mapped.reads in : " $1.MAPQ0.bam
samtools view -c $1.MAPQ0.bam


echo "removing low-MAPQ read "
samtools view -b -q 5 $1.MAPQ0.bam > $1.MAPQ5.bam	# keep only mapped reads

rm $1.sam


echo "mapped reads in : " $1.MAPQ5.bam
samtools view -c $1.MAPQ5.bam


#echo "converting into bam " 
#samtools view -Sb $1.MAPQ5.sam > $1.MAPQ5.bam
#rm $1.MAPQ5.sam

echo "removing duplicate "
### Remove duplicates
echo "soring by picard..."
java -Xms16000m -jar picard.jar SortSam INPUT=$1.MAPQ5.bam OUTPUT=sorted.$1.bam SORT_ORDER=coordinate

echo "Removing duplicates..."
java -Xms16000m -jar picard.jar MarkDuplicates INPUT=sorted.$1.bam OUTPUT=dup.removed.$1.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true METRICS_FILE=$1.dup.metrics
echo "reads_in_dup.removed.bam"
samtools view -c dup.removed.$1.bam 


rm $1.MAPQ5.bam
rm sorted.$1.bam

	
echo "checking_flagstat_in_dup.removed.bam"
samtools flagstat dup.removed.$1.bam


echo "sorting by read name"
samtools sort dup.removed.$1.bam -o namesorted.$1.bam

rm dup.removed.$1.bam

echo "checking_flagstat_in_only.paired.bam"
samtools flagstat namesorted.$1.bam  

echo "checking insert length"
echo "soring by picard coordinate before check insert length"
java -Xms16000m -jar picard.jar SortSam INPUT=namesorted.$1.bam OUTPUT=sorted.$1.only.paired.bam SORT_ORDER=coordinate
java -Xms16000m -jar picard.jar CollectInsertSizeMetrics I=sorted.$1.only.paired.bam O=$1.size_metrics.txt H=$1.size_histogram.pdf

samtools index sorted.$1.only.paired.bam	

##read counts on genes by htseq
htseq-count -f bam -r name -s no sorted.$1.only.paired.bam /data/aryee/caleb/hemeATAC/amit/rnaseq/2016.5.mm10.refGene.gtf > $1.count

done

