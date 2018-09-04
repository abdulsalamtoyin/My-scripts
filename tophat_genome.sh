#!/usr/bin/env bash

###############################################################################
# Basic mapping pipeline for paired end data using tophat2 then preparing bam files for IGV
###############################################################################

###############################################################################
# BEFORE RUNNING SCRIPT DO THE FOLLOWING:
###############################################################################

# Make sure modules are loaded:
#       tophat
#       samtools
#       subread

# Check the location of the hard variables

###############################################################################
# Hard variables
###############################################################################

echo "$(date +"%m-%d-%Y_%H:%M") Starting pipeline"

# read quality control

##for i in $(ls $fastq/*_001.fastq.gz | rev | cut -c 16- | rev | uniq)

##do

##java -jar /opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 30 ${i}R1_001.fastq.gz ${i}R2_001.fastq.gz \
##${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" \
##ILLUMINACLIP:/opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:14

##mv ${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" $Dir

##done


# 	Relevant PATHS
#	Reference Genome
Gen="/data/biodata/Zmays/Zmays_chr1-10.fa"
#	trimed Sequencing Reads
reads="/syn02iscsi/toyin/STRIGA-ANALYSIS/trim_reads"
#       bowtie2 indices
INDEX="/data/biodata/Zmays/Zmays_chr1-10"
#       GTF annotation file
GTF="/data/biodata/Zmays/Zmays_284_6a.gene_exons.gff3"
GenDir="/syn02iscsi/toyin/STRIGA-ANALYSIS"
#       command to display software versions used during the run
##modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)

###############################################################################
###############################################################################


usage="
    USAGE
       step1:   load the following modules: bbmap bowtie samtools python (default version)
       step2:   bash bowtie2_pipeline.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing read files (fastq.gz format)

        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_r1' or '_r2' depending on respective orientation.
    "

# Set default parameters
PAIRED=true


# Parse command line parameters

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi


while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--dir)
        DIR="$2"
        shift
        ;;
        -p|--paired)
        PAIRED="true"
        ;;
    esac
shift
done

# Remove trailing "/" from input directory if present

if [[ ${DIR:(-1)} == "/" ]]; then
    DIR=${DIR::${#DIR}-1}
fi

# Print out loaded modules to keep a record of which software versions were used in this run

echo "$modules"

###############################################################################
###############################################################################

# Prepare directories

if [ ! -d "star_out" ]; then
    mkdir star_out
fi

if [ ! -d "bam" ]; then
    mkdir bam
fi

if [ ! -d "count" ]; then
    mkdir count
fi

##STAR  --runMode genomeGenerate --runThreadN 30 --genomeDir $GenDir/Genome --genomeFastaFiles $Gen



echo "$(date +"%m-%d-%Y_%H:%M") Starting pipeline"

for file in ${DIR}/*.fastq; do

    SKIPFILE=false

    if [[ $PAIRED == "true" ]]; then

        # Paired end

        if [[ ${file:(-19)} == "R1.trimmed_PE.fastq" ]]; then

            FBASE=$(basename $file .trimmed_PE.fastq)
            BASE=${FBASE%_R1}

            # Map reads using STAR

            if [ ! -d ./star_out/${BASE} ]; then
                mkdir ./star_out/${BASE}
            fi

            echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${BASE} with STAR... "

            STAR --genomeDir $GenDir/Genome --runThreadN 32 —readFilesIn ${DIR}/${BASE}_R1.trimmed_PE.fastq ${DIR}/${BASE}_R2.trimmed_PE.fastq --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outWigType bedGraph read2 --outWigStrand Stranded --outWigNorm None --outFileNamePrefix ./star_out/${BASE}
        else

            # Avoid double mapping by skipping the r2 read file

            SKIPFILE=true
        fi
    else

        # Single end

        BASE=$(basename $file .trimmed_PE.fastq)
    fi

        # Map reads using Tophat

        if [ ! -d ./star_out/${BASE} ]; then
            mkdir ./star_out/${BASE}
        fi

        echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${BASE} with STAR... "

        STAR --genomeDir $GenDir/Genome --runThreadN 32 —readFilesIn ${DIR}/${BASE}_R1.trimmed_PE.fastq ${DIR}/${BASE}_R2.trimmed_PE.fastq --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outWigType bedGraph read2 --outWigStrand Stranded --outWigNorm None --outFileNamePrefix ./star_out/${BASE}
    if [[ $SKIPFILE = false ]]; then

        echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${BASE}"

        echo "$(date +"%m-%d-%Y_%H:%M") Sorting and indexing ${BASE}.bam"

        samtools sort -o ./bam/${BASE}_sorted.bam ./tophat_out/${BASE}/accepted_hits.bam

        samtools index ./bam/${BASE}_sorted.bam
    fi
done

echo $(date +"%m-%d-%Y_%H:%M")" Counting reads with featureCounts... "

ARRAY=()

for file in ./bam/*_sorted.bam
do

ARRAY+=" "${file}

done

##featureCounts -a ${GTF} -o ./count/counts.txt -T 4 -t exon -g transcript_id${ARRAY}



	##Performing Transcript Recunstruction with Cufflinks#####
echo $(date +"%m-%d-%Y_%H:%M")" Performing Transcript Recunstruction with Cufflinks......."

cufflinks -p 30 -o ${BASE}_cufflinks --no-update-check -N -g $anno ./star_out/${BASE}/accepted_hits.bam

echo "Tophat and Cufflinks DONE."

## Merge cufflinks outputs using Cuffmerge.
echo "====              Now, starting Cuffmerge...            ===="

## List all transcript outputs from Cufflinks in a text file.
ls -d1 ${BASE}_cufflinks/ | sed 's/^/.\//' | sed 's/$/transcripts.gtf/' > assemblies.txt

# Vars to identify Log files.
today=`date +'%Y-%m-%d'`
user=`pwd | cut -d "/" -f 3`

cuffmerge -p 32 -o CUFFMERGE -g $anno -s  assemblies.txt &> Log_Merge_${today}_${user}.txt

echo "Cuffmerge DONE."

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

cuffdiff -o CUFFDIFF --no-update-check -b $Gen -p 32 -L $minimalnames -u CUFFMERGE/merged.gtf $toplist &> Log_Cuffdiff_${today}_${user}.txt

