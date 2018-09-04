#!/bin/bash

#################
#       Abdulsalam Toyin 2014
#
#       Command line: sh rna-seq.sh <(1)new-dir> <(2)read-seq1> <(3)read-seq2> <(4)ref-seq> <(5)ref-index> <(6)ref-annot> <(7)tophat_out> <(8)Cufflinks_out> <(9)cuffmerge> <(10)cuffdiff_out>
#
###################
Gen=/data/biodata/Zmays/Zmays_chr1-10.fa
Index=/data/biodata/Zmays/Zmays_chr1-10
anno=/data/biodata/Zmays/Zmays_284_6a.gene_exons.gff3
reads=/syn02iscsi/toyin/STRIGA-ANALYSIS/trim_reads

echo "$(date +"%m-%d-%Y_%H:%M") Starting pipeline"

dir=`pwd`

##read quality control

##mkdir $dir/$1

##cp $2 $1/$2
##cp $3 $1/$3

##cd $1

# read quality control

##for i in $(ls $fastq/*_001.fastq.gz | rev | cut -c 16- | rev | uniq)

##do

##java -jar /opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 30 ${i}R1_001.fastq.gz ${i}R2_001.fastq.gz \
##${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" \
##ILLUMINACLIP:/opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:14

##mv ${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" $Dir

##done


## Mapping based transcript reconstruction

##cd /data/biodata/Cassava/genome/V6/masked
##bowtie2-build cassava_chr_1-18.fa cassava_chr_1-18
##samtools faidx cassava_chr_1-18.fa

# Prepare directories

echo "$(date +"%m-%d-%Y_%H:%M") Mapping with STAR... "        

##STAR  --runMode genomeGenerate --runThreadN 30 --genomeDir $dir/Genome --genomeFastaFiles $Gen

for i in $(ls $reads/*.trimmed_PE.fastq | rev | cut -c 20- | rev | uniq)
do
STAR --genomeDir $dir/Genome  --runThreadN 30 --readFilesIn ${i}"R1.trimmed_PE.fastq" ${i}"R2.trimmed_PE.fastq" --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outWigType bedGraph read2 --outWigStrand Stranded --outWigNorm None --outFileNamePrefix ${i}Star
done

for i in $(ls $reads/*Aligned.sortedByCoord.out.bam*)
do
cufflinks -p 30 -o ${i}_cufflinks --no-update-check -N -g $anno ${i}Star"Aligned.sortedByCoord.out.bam"
done

##tophat2 -i 12000 --no-mixed --no-coverage-search --max-multihits 1 -o $6 -p 30 ${INDEX} $2 $3

## Merge cufflinks outputs using Cuffmerge.
echo "====              Now, starting Cuffmerge...            ===="

cd $reads

## List all transcript outputs from Cufflinks in a text file.
ls -d1 cufflinks*/ | sed 's/^/.\//' | sed 's/$/transcripts.gtf/' > assemblies.txt

# Vars to identify Log files.
today=`date +'%Y-%m-%d'`
user=`pwd | cut -d "/" -f 3`

cuffmerge -p 30 -o CUFFMERGE -g $anno -s  assemblies.txt &> Log_Merge_${today}_${user}.txt

echo "Cuffmerge DONE."


## Calculate transcript expression, splicing, and promoter use with Cuffdiff.

echo "====             Now, starting CuffDiff...            ===="

# Generate CuffDiff args:
minimalnamesSep=`echo $minimalnames | sed 's/,/ /g'`
toplist=''

for miniID in $minimalnamesSep
    do
        temp=`grep "$miniID" replicateSamples.txt| cut -d " " -f 2 | while read -r line; do echo -n "${i}Star"Aligned.sortedByCoord.out.bam","; done | sed 's/,$//'`
        toplist="$toplist $temp"
    done


# Vars to identify Log files.
today=`date +'%Y-%m-%d'`
user=`pwd | cut -d "/" -f 3`

cuffmerge -p 30 -o CUFFMERGE -g $anno -s  assemblies.txt &> Log_Merge_${today}_${user}.txt

echo "Cuffmerge DONE."


## Calculate transcript expression, splicing, and promoter use with Cuffdiff.

echo "====             Now, starting CuffDiff...            ===="

# Generate CuffDiff args:
minimalnamesSep=`echo $minimalnames | sed 's/,/ /g'`
toplist=''

for miniID in $minimalnamesSep
    do
        temp=`grep "$miniID" replicateSamples.txt| cut -d " " -f 2 | while read -r line; do echo -n "tophat_out_$line/accepted_hits.bam,"; done | sed 's/,$//'`
        toplist="$toplist $temp"
    done

cuffdiff -o CUFFDIFF --no-update-check -b $Gen -p 30 -L $minimalnames -u CUFFMERGE/merged.gtf $toplist &> Log_Cuffdiff_${today}_${user}.txt




