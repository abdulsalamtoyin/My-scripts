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

#       bowtie2 indices
INDEX="/nas02/home/s/f/sfrenk/proj/seq/WS251/genome/bowtie2/genome"

#       GTF annotation file
GTF="/nas02/home/s/f/sfrenk/proj/seq/WS251/genes.gtf"

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
PAIRED=false


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

if [ ! -d "tophat_out" ]; then
    mkdir tophat_out
fi

if [ ! -d "bam" ]; then
    mkdir bam
fi

if [ ! -d "count" ]; then
    mkdir count
fi

echo "$(date +"%m-%d-%Y_%H:%M") Starting pipeline"

for file in ${DIR}/*.fastq.gz; do
    
    SKIPFILE=false

    if [[ $PAIRED == "true" ]]; then
            
        # Paired end

        if [[ ${file:(-11)} == "r1.fastq.gz" ]]; then
        
            FBASE=$(basename $file .fastq.gz)
            BASE=${FBASE%_r1}

            # Map reads using Tophat
                
            if [ ! -d ./tophat_out/${BASE} ]; then
                mkdir ./tophat_out/${BASE}
            fi

            echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${BASE} with Tophat... "        
            tophat -i 12000 --no-mixed --no-coverage-search --max-multihits 1 -o ./tophat_out/${BASE} -p 4 ${INDEX} ${DIR}/${BASE}_r1.fastq.gz ${DIR}/${BASE}_r2.fastq.gz

        else

            # Avoid double mapping by skipping the r2 read file
                
            SKIPFILE=true
        fi
    else

        # Single end

        BASE=$(basename $file .fastq.gz)
    fi

        # Map reads using Tophat
                
        if [ ! -d ./tophat_out/${BASE} ]; then
            mkdir ./tophat_out/${BASE}
        fi

        echo "$(date +"%m-%d-%Y_%H:%M") Mapping ${BASE} with Tophat... "        
        tophat -i 12000 --no-mixed --no-coverage-search --max-multihits 1 -o ./tophat_out/${BASE} -p 4 ${INDEX} ${DIR}/${BASE}_r1.fastq.gz ./${DIR}/${BASE}.fastq.gz

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

featureCounts -a ${GTF} -o ./count/counts.txt -T 4 -t exon -g transcript_id${ARRAY}
