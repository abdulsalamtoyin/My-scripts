#!/bin/bash

##############################################################################################
### Created By Abdulsalam Toyin, IITA 2016						######
### Striga RNA-Seq Analysis								######
###											######
##############################################################################################

WD=/syn02iscsi/toyin/STRIGA-ANALYSIS/
DIR=/data/biodata/Zmays/StrigaMaizeRNA-seq
anno=/data/biodata/Zmays/Zmays_284_6a.gene_exons.gff3
GEN=/data/biodata/Zmays/Zmays_chr1-10.fa
GenDir=/syn02iscsi/toyin/STRIGA-ANALYSIS

#####Quality control##########

##read quality control

#Make sure w are in the right directory

# FastQC
##mkdir results-fastqc

##sample_file=$(sed -n "$SLURM_ARRAY_TASK_ID"p fastq.list)


##fastqc -o results-fastqc -t 32 ${sample_file}


mkdir $WD/trim_seq
fastq=/data/biodata/Zmays/StrigaMaizeRNA-seq/
Dir=/syn02iscsi/toyin/STRIGA-ANALYSIS/trim_reads

for i in $(ls $fastq/*_001.fastq.gz | rev | cut -c 16- | rev | uniq)
do
java -jar /opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 30 ${i}R1_001.fastq.gz ${i}R2_001.fastq.gz \
${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" \
ILLUMINACLIP:/opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:14
mv ${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" $Dir
done

##for file in $DIR/*_R1_001.fastq.gz 
##do
##base=`basename $file _R1_001.fastq.gz`
##if [ ! -f $DIR/$base"_R1_001_trim.fastq.gz" ]; then
##sickle pe -f $DIR/$base"_R1_001.fastq.gz" -r $DIR/$base"_R2_001.fastq.gz" -o $WD/trim_seq/$base"_R1_001_trim.fastq.gz" -p $WD/trim_seq/$base"_R2_001_trim.fastq.gz" -s $WD/trim_seq/$base"_R1_001_unpair.trim.fq.gz" -t sanger 
##fi
##done

##cd $WD/trim_seq

##sample_file=$(sed -n "$SLURM_ARRAY_TASK_ID"p fastq.list)

# FastQC
##fastqc -o results-fastqc ${sample_file}

##cd ..

#index genome
##
##for file in $GEN;
##if [ ! -f $GEN.bt2 ];
##then
##bowtie2-build $GEN.fa $GEN
##fi
##done

perl $WD/runTuxedo.pl 

##STAR  --runMode genomeGenerate --runThreadN 30 --genomeDir $GenDir/Genome --genomeFastaFiles $GEN.fasta

##for i in $(ls results-QC/*_001_trim.fastq | rev | cut -c 19- | rev | uniq)
##do
##STAR --genomeDir $GenDir/Genome  --runThreadN 30 --readFilesIn ${i}_R1_001_trim.fastq ${i}_R2_001_trim.fastq --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outWigType bedGraph read2 --outWigStrand Stranded --outWigNorm None --outFileNamePrefix ${i}Star
##done

==================================================
tuxedo workflow script ends here
==================================================

R CMD BATCH cummeRbund.R
Rscript cummeRbund.R
