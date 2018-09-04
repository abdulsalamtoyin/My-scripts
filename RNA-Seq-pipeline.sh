#!/bin/bash
##
##      Pipeline for RNA-seq using TopHat/Cufflinks protocol.
##      Aligned to NCBIm37 / mm9 with Ensembl annotation.
##
##                                              Cynthia Alexander Rascon 2014-06
##
################################################################################

usage() { echo "Usage: $0 [-help]";}

help="
# Run within a directory containing (only) the folders for all samples to be 
# compared. Each folder containing zip files of reads in *fastq format.
#
# Obtain TopHat alignments to Ensembl mm9, BAM (+ index) of sorted-
# non redundant tags and TDF file, as well as CuffLinks, CuffMerge and CuffDiff
# outputs.
#
# Needs: Tophat, Cufflinks, Samtools and igvtools. 
#        Tophat reference genome, Cufflinks GTF guide, Cuffmerge
# 	 annotation GTF and genomic seqs references, and Chr sizes file. 
#
# Infiles: (no args) Will assume each sample fastq files are located in 
# 	separated folders. Samples will be renamed and organized through
# 	stdin. 
# Outfiles: 
#	RAW: (dir) *.fastq files after unzip and concatenation.
#	TOPHAT (dir) all TopHat outputs.	
#	BAM: (dir) Tophat alignment in  BAM format of sorted, non redundant tags
#		   *.mm9Ensemblalign.sorted.nonred.bam
#		   BAM indexes
#		   *.mm9Ensemblalign.sorted.nonred.bam.bai
#	CUFFLINKS: (dir) all Cufflinks outputs.
#	CUFFMERGE: (dir) all Cuffmerge outputs.
#	CUFFDIFF: (dir) all Cuffdiff outputs.
#	TDF: tdf file ready for igv visualization.
#	LOGS: (dir) Log files generated for this pipeline.  
#	"


################################################################################
#			Rename and Organize samples
################################################################################

# Get arguments
while getopts ":h" o; do
	case "${o}" in
                h)
                 echo "$help"
                 ;;
	esac
done

## List all samples to be compared, change names by user input.

[ -f  samples.txt ] && rm samples.txt   # Just as precaution
[ -f  renamedSamples.txt ] && rm renamedSamples.txt
[ -f  replicateSamples.txt ] && rm replicateSamples.txt
ls -d */ | sed 's/\/$//' > samples.txt

echo -e "\nTo rename sample, type in new name + [ENTER], else just hit [ENTER].
     *Only [a-z][A-Z][0-9][.][-][_] allowed.
     *Consider a nomenclature like yyyymm_sample_repNum_ownerInitials
     i.e. 201401_sampleA_rep1_CA\n"

Samples="./samples.txt"
while IFS= read -r name
        do
                echo  -en "Current name : $name  \nNew name : "
                read  newname </dev/tty
                echo -e "$name ${newname:=$name}" >> renamedSamples.txt
        done <"$Samples"

## Get sample names and assign replicates
echo -ne "\nMinimal sample names (disregard replicates) separated by comma w/o
 spaces (i.e. HFSCs,TA,DP): "
read minimalnames

echo -ne "\nAssign each replicate to a sample [$minimalnames]
 (i.e. TA for Sample-TA-Rep001).\n"
Samples="./renamedSamples.txt"
while IFS= read -r line
        do
                echo -n "$line corresponds to : "
                read rep  </dev/tty
                echo "$line $rep" >> replicateSamples.txt
        done <"$Samples"

## Unzip and concatenate samples each in a single .fastq file.

echo "====             Uncompressing files...             ===="
gunzip */*

echo "====             Concatenating files...             ===="
Samples="./renamedSamples.txt"
while IFS=' ' read -r name newname
	do 
		cat $name/*fastq > $newname.fastq
	done < "$Samples"

################################################################################
#	       		    TopHat and CuffLinks
################################################################################

## Align and estimate transcripts abundance for each sample using Tophat and 
## Cufflinks.

echo "====          Start Tophat and Cufflinks...         ===="

for sample in *fastq*;
	do
		echo "Processing $sample ..."
		tophat -p4 -o tophat_out_${sample/.fastq/} ~/Vol/RefGenomes/Bowtie2Index_mm9/genome $sample
		cufflinks -p4 -o cufflinks_${sample/.fastq/} --no-update-check -N -g ~/Vol/RefGenomes/RNAref_mm9_Ensembl_genes.gtf tophat_out_${sample/.fastq/}/accepted_hits.bam
	done

echo "Tophat and Cufflinks DONE."

## Merge cufflinks outputs using Cuffmerge.
echo "==== 		Now, starting Cuffmerge...	      ===="

## List all transcript outputs from Cufflinks in a text file.
ls -d1 cufflinks*/ | sed 's/^/.\//' | sed 's/$/transcripts.gtf/' > assemblies.txt

# Vars to identify Log files.
today=`date +'%Y-%m-%d'`
user=`pwd | cut -d "/" -f 3`

cuffmerge -p4 -o CUFFMERGE -g ~/Vol/RefGenomes/RNAref_mm9_Ensembl_genes.gtf -s ~/Vol/RefGenomes/DNAref_mm9_Ensembl_Chromosomes assemblies.txt &> Log_Merge_${today}_${user}.txt

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

cuffdiff -o CUFFDIFF --no-update-check -b ~/Vol/RefGenomes/WholeGenome_mm9_Ensembl.fa -p 4 -L $minimalnames -u CUFFMERGE/merged.gtf $toplist &> Log_Cuffdiff_${today}_${user}.txt

################################################################################
#		       BAM, bai, and TDF files + Clean up.
################################################################################

## Sort, filter and index tophat aligned reads. Plus generate TDF files (for IgV).
ls -d1 tophat*/ | sed 's/$/accepted_hits.bam/' > aligned.txt

TophatOuts="./aligned.txt"
while IFS= read -r line
	do
		temp=${line#.\/tophat_out_}
		sample=${temp%\/accepted_hits.bam}
		samtools sort $line temp_${sample}.mm9Ensemblalign.sorted
        samtools rmdup -s temp_${sample}.mm9Ensemblalign.sorted.bam ${sample}.mm9Ensemblalign.sorted.nonred.bam
        samtools index ${sample}.mm9Ensemblalign.sorted.nonred.bam
		igvtools count ${sample}.mm9Ensemblalign.sorted.nonred.bam ${sample}.mm9Ensemblalign.sorted.nonred.tdf ~/Vol/RefGenomes/mm9.chrom.sizes
	done <"$TophatOuts"

## Clean your mess.
rm temp*mm9Ensemblalign.sorted.bam
mkdir BAM; mv *.bam BAM; mv *.bai BAM
mkdir RAW; mv *fastq RAW
mkdir TOPHAT; mv tophat* TOPHAT
mkdir CUFFLINKS; mv cuff* CUFFLINKS
mkdir LOGS; mv samples.txt LOGS; mv renamedSamples.txt LOGS
mv *.txt LOGS; mv igv.log LOGS; mv replicateSamples.txt LOGS
mkdir TDF; mv *.tdf TDF

echo "====			FINISHED	            ===="
