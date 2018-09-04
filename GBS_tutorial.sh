#!/bin/bash

# Install gedit and dos2unix
sudo pip install gedit dos2unix

# Create a working directory, make writable, and download the data and documentation:
sudo mkdir /mnt/GBS
sudo chmod 777 /mnt/GBS
cd /mnt/GBS
wget http://ncsubit815.s3-website-us-east-1.amazonaws.com/SRR072188.fastq
wget http://ncsubit815.s3-website-us-east-1.amazonaws.com/GBSbarcodes
wget http://www.maizegenetics.net/tassel/docs/TasselPipelineGBS.pdf

# Downloading files can take a while, so be patient.
# Spend the time reading the rest of the script functions and 
# understanding what each command-line does, using the GBS pipeline
# document as a guide. 
# The GBS documentation recommends the following sub-directories
# for data management: 

mkdir fastq              #  for raw data files, one file per flowcell lane
mkdir tagCounts          # for output from FastqToTagCountPlugin
mkdir mergedTagCounts    # for output from MergeMultipleTagCountPlugin
mkdir tbt                # for output from QseqToTBTPlugin OR FastqToTBTPlugin
mkdir mergedTBT          # for output from MergeTagsByTaxaFilesPlugin

# The path to the perlscript is /usr/share/java/tassel-3.0/run_pipeline.pl
# A description of the GBS experiment is available at 
# http://sra.dnanexus.com/experiments/SRX030801

# The required format for the barcode file is 8 tab-delimited columns:
# flowcell lane barcode sampleID plateID row col blank
# The eighth column must be text, but can contain the word 
# blank repeated as necessary - this change is not yet 
# documented in the literature, but was explained to me in an email
# exchange with Jeff Glaubitz and James Harriman.

# The fastq file we are using contains headers with SRA format, 
# starting with @SRR072188.[read#], space, then
# instrument name_runname:lane:tile:x:y length=86. 
# The TASSEL pipeline will not accept this format - it
# expects header lines to be in standard Illumina format,
# either pre-CASAVA v1.8 or CASAVA v1.8 (see Wikipedia page
# on Fastq format for details of these header line formats).

# The pipeline also expects files to be named 
# flowcellID_laneID_fastq.txt or ..._fastq.gz,
# and the flowcell ID in the sequence header 
# lines has to match the file name and the flowcell
# ID given in the barcode file.

# Use sed to change the sequence header lines to 
# match the example flowcell ID in the GBS document.

sed -r 's/SRR072188.[0-9]{1,} HWUSI-EAS690_0010/EAS690:10:42A87AAXX/g' SRR072188.fastq > sequences.fq

# This editing takes a while also. Be patient...

# Compress the sequence file with gzip, then move to 
# the fastq folder and re-name:

gzip sequences.fq
mv sequences.fq.gz fastq/42A87AAXX_7_fastq.gz

######   increase memory allocated to Java by changing
#run_pipeline.pl script using sed, after re-setting dir permissions:
sudo chmod 777 /Users/Toyin/Desktop/training/tassel

# change -Xms512m -Xmx1536m in last line to -Xms15g -Xmx30g
sed 's/-Xms512m -Xmx1536m/-Xms10g -Xmx15g/' /Users/Toyin/Desktop/training/tassel/run_pipeline.pl > /Users/Toyin/Desktop/training/tassel/run_pipeline2.pl

# Make new file executable; convert to unix format
chmod 777 /Users/Toyin/Desktop/training/tassel/run_pipeline2.pl

dos2unix /Users/Toyin/Desktop/training/tassel/run_pipeline2.pl


# Run the TASSEL pipeline plugins to process the data

# The -c 5 requires a tag to be detected 5 times before 
# it is kept - this helps filter out sequencing errors.
# If merging data from multiple lanes, where the same 
# sample is present in different files, this can be
# reduced to avoid losing data - a tag that appears 
# once in one lane may appear more times in another lane.

/Users/Toyin/Desktop/training/tassel/run_pipeline2.pl -fork1 -FastqToTagCountPlugin -i fastq -k GBSbarcodes -e ApeKI -c 5 -o tagCounts -endPlugin -runfork1

# The MergeMultipleTagCount plugin merges data from multiple
#  flowcells or lanes into one file
# We have only a single file, but we need the data in the format
# produced by this plugin, so we run the plugin anyway.
#  The -c option can be used here instead of
# in the FastqToTagCountPlugin step,
# but for a single file, it is not necessary to use it both places.

/Users/Toyin/Desktop/training/tassel/run_pipeline2.pl -fork1 -MergeMultipleTagCountPlugin -i ./tagCounts -o ./mergedTagCounts/MasterTags.cnt -c 48 -endPlugin -runfork1

# The -y option exports tag counts from 0 to 127, instead of the
# default 0 (absent) or 1 (present)
# The -c option at this step is the number of individuals in which a
# tag must be detected to be output
# -c 20 specifies that only tags found in at least 20 individuals are exported.

/Users/Toyin/Desktop/training/tassel/run_pipeline2.pl -fork1 -FastqToTBTPlugin -i fastq -k GBSbarcodes -e ApeKI -o tbt -t mergedTagCounts/MasterTags.cnt -y -endPlugin -runfork1

# The output from the FastqToTBTPlugin is in binary format -
# convert it to a text format for analysis by other programs.

/Users/Toyin/Desktop/training/tassel/run_pipeline2.pl -fork1 -BinaryToTextPlugin -i tbt/42A87AAXX_7.tbt.byte -o tbt/TagsByTaxa.txt -t TBTByte -endPlugin -runfork1

# The TagsByTaxa.txt file can be loaded into R and summarized using
# various tools in R.
