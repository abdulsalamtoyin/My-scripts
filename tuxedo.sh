## Align and estimate transcripts abundance for each sample using Tophat and 
## Cufflinks.

Gen=/data/biodata/Zmays/Zmays_chr1-10.fa
Index=/data/biodata/Zmays/Zmays_chr1-10
anno=/data/biodata/Zmays/Zmays_284_6a.gene_exons.gff3
reads=/syn02iscsi/toyin/STRIGA-ANALYSIS/trim_reads

echo "====          Start Tophat and Cufflinks...         ===="

for sample in $reads/*fastq*;
	do
		echo "Processing $sample ..."
		tophat2 -p 30 -o tophat_out_${sample} $Index $sample
		cufflinks -p 30 -o ${sample}_cufflinks --no-update-check -N -g $anno ${sample}_tophat_out/accepted_hits.bam
	done

echo "Tophat and Cufflinks DONE."

## Merge cufflinks outputs using Cuffmerge.
echo "==== 		Now, starting Cuffmerge...	      ===="

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
        temp=`grep "$miniID" replicateSamples.txt| cut -d " " -f 2 | while read -r line; do echo -n "tophat_out_$line/accepted_hits.bam,"; done | sed 's/,$//'`
        toplist="$toplist $temp"
    done

cuffdiff -o CUFFDIFF --no-update-check -b $Gen -p 30 -L $minimalnames -u CUFFMERGE/merged.gtf $toplist &> Log_Cuffdiff_${today}_${user}.txt

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

