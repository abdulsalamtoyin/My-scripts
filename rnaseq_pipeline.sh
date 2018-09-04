#!/bin/bash

echo "
 #==========================================================================#
 #     RNA-SEQ PIPELINE FOR THE BIOSCIENCE CENTER,							#
 #	   INTERNATIONAL INSTITIUTE OF TROPICAL AGRICULTURE						#
 #	    				AUGUST	2016											#
 #				      ABDULSALAM TOYIN										#
 #    																		#
 #   																		#
 #  																		#
 #  																   		#
 #   																	    #
 #                                                                          #
 #==========================================================================#
"

##########################################################

# initialize some variables:

echo "Checking dependencies..."

#check for java
if ! which java ; then
	echo "Could not find java in current directory or in PATH"
	exit 1
fi

#check for Rscript
if ! which Rscript ; then
	echo "Could not locate the Rscript engine in current directory or PATH"
	exit 1
fi

PYTHON='/cccbstore-rc/projects/cccb/apps/bin/python2.7'
#check for python- for regex syntax need 2.7 or greater!
if ! which $PYTHON ; then
	echo "Could not access python located at $PYTHON.  Require version 2.7 or greater"
	exit 1
fi


#  a function that prints the proper usage syntax
function usage
{
	echo "**************************************************************************************************"
	echo "usage: 
		-d | --dir sample_directory 
		-s | --samples samples_file 
		-g | --genome
		-c | --contrasts contrast_file (optional) 
		-config (optional, configuration file-- if not given, use default)
                -noalign (optional, default behavior is to align) 
                -paired (optional, default= single-end) 
		-a | --aligner (optional, default is STAR)
		-test (optional, for simple test)"
	echo "**************************************************************************************************"
}


#  expects the following args:
#  $1: a file containing the sample and condition (tab-separated) (one per line)
function print_sample_report
{
	if [ -e $1 ]; then
		echo ""
		printf "%s\t%s\n" Sample Condition
		 while read line; do
			printf "%s\t%s\n" $(echo $line | awk '{print $1}') $(echo $line | awk '{print $2}')
		done < $1
		echo ""
		echo ""
	else
		echo "Sample file ("$1") was not found."
	fi
}


#  expects the following args:
#  $1: a file containing the base/control and experimental/case condition (tab-separated) (one per line)
function print_contrast_report
{
	if [ -e $1 ]; then
		echo ""
		printf "%s\t%s\n" Base/Control Case/Condition
	        while read line; do
			printf "%s\t%s\n" $(echo $line | awk '{print $1}') $(echo $line | awk '{print $2}')
		done < $1
		echo ""
		echo ""
	else
		echo "Contrast file ("$1") was not found."
	fi
}



##########################################################

#read input from commandline:
while [ "$1" != "" ]; do
	case $1 in
		-c | --contrasts )
			shift
			CONTRAST_FILE=$1
			;;
		-d | --dir )
			shift
			PROJECT_DIR=$1
			;;
		-g | --genome )
			shift
			ASSEMBLY=$1
			;;
		-s | --samples )
			shift
			SAMPLES_FILE=$1
			;;
                -config )
                        shift
                        CONFIG=$1
                        ;;
		-a | --aligner )
			shift
			ALIGNER=$1
			;;
		-noalign )
			ALN=0
			;;
		-paired )
			PAIRED_READS=1
			;;
		-h | --help )
			usage
			exit
			;;
		-test )
			TEST=1
			;;
		* )
			usage
			exit 1
	esac
	shift
done

############################################################


############################################################

#check that we have all the required input:

if [ "$PROJECT_DIR" == "" ]; then
    echo ""
    echo "ERROR: Missing the project directory.  Please try again."
    usage
    exit 1
fi

if [ "$SAMPLES_FILE" == "" ]; then
    echo ""
    echo "ERROR: Missing the samples file.  Please try again."
    usage
    exit 1
fi


if [ "$ASSEMBLY" == "" ]; then
    echo ""
    echo "ERROR: Missing the genome.  Please try again."
    usage
    exit 1
fi

if [ "$PAIRED_READS" == "" ]; then
    echo ""
    echo "Single-end sequencing protocol specified by default."
    PAIRED_READS=0
fi


# Set some default parameters if they were not explicitly set in the input args:

#if ALN was not set, then -noalign flag was NOT invoked, meaning we DO align
if [ "$ALN" == "" ]; then
    ALN=1       
fi

#if TEST was not set, then do NOT test
if [ "$TEST" == "" ]; then
    TEST=0
fi

#if the aligner was not explicitly set, default to STAR
if [ "$ALIGNER" == "" ]; then
    ALIGNER=STAR
fi


#if no configuration file was given, then use teh default one
if [ "$CONFIG" == "" ]; then
    CONFIG=/cccbstore-rc/projects/cccb/pipelines/RNA_Seq_pipeline/config.txt
    # double check that the configuration file exists:
    if [[ ! -f "$CONFIG" ]]; then
        echo "Could not locate a configuration file at "$CONFIG
        exit 1
    fi 
    echo ""
    echo "Using the default configuration file located at: "$CONFIG
fi


############################################################

#read-in the non-dynamic configuration parameters (and export via set to have these as environment variables):
# !!! important-- import the configuration file !!!
set -a
source $CONFIG
set +a


############################################################
# After inputs have been read, proceed with setting up parameters based on these inputs:

# construct the full paths to some files by prepending the project directory:
export VALID_SAMPLE_FILE=$PROJECT_DIR'/'$VALID_SAMPLE_FILE
export DESIGN_MTX_FILE=$PROJECT_DIR'/'$DESIGN_MTX_FILE
export COUNTS_DIR=$PROJECT_DIR'/'$COUNTS_DIR
export ASSEMBLY
export PROJECT_DIR
export SAMPLES_FILE
export PAIRED_READS
export ALIGNER

#############################################################

#identify the correct genome files to use
if [[ "$ASSEMBLY" == hg19 ]]; then
    GTF=/cccbstore-rc/projects/db/genomes/Human/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.gtf
    GENOMEFASTA=/cccbstore-rc/projects/db/genomes/Human/GRCh37.75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
    SNAPR_GENOME_INDEX=/cccbstore-rc/projects/db/genomes/Human/GRCh37.75/SNAPR/index-dir
    SNAPR_TRANSCRIPTOME_INDEX=/cccbstore-rc/projects/db/genomes/Human/GRCh37.75/SNAPR/transcriptome-dir
    STAR_GENOME_INDEX=/cccbstore-rc/projects/db/genomes/Human/GRCh37.75/STAR_INDEX  
    GTF_FOR_RNASEQC=/cccbstore-rc/projects/db/genomes/Human/GRCh37.75/GTF/Homo_sapiens.GRCh37.75.transcript_id.chr_trimmed.gtf
elif [[ "$ASSEMBLY" == mm10 ]]; then
    GTF=/cccbstore-rc/projects/db/genomes/Mm/build38/Mus_musculus.GRCm38.75.chr_trimmed.gtf
    GENOMEFASTA=/cccbstore-rc/projects/db/genomes/Mm/build38/mm10.fa
    SNAPR_GENOME_INDEX=/cccbstore-rc/projects/db/genomes/Mm/build38/snapr_genome_index
    SNAPR_TRANSCRIPTOME_INDEX=/cccbstore-rc/projects/db/genomes/Mm/build38/snapr_transcriptome_index
    STAR_GENOME_INDEX=/cccbstore-rc/projects/db/genomes/Mm/build38/STAR_INDEX
    GTF_FOR_RNASEQC=/cccbstore-rc/projects/db/genomes/Mm/build38/Mus_musculus.GRCm38.75.transcript_id.chr_trimmed.gtf
else
    echo "Unknown or un-indexed genome."
    exit 1
fi

export GTF
#############################################################

# parameters based on the aligner selected:
if [[ $ALIGNER == $STAR  ]]; then
    ALN_DIR_NAME=$STAR_ALIGN_DIR
    ALIGN_SCRIPT=$STAR_ALIGN_SCRIPT
    GENOME_INDEX=$STAR_GENOME_INDEX
    TRANSCRIPTOME_INDEX=    #nothing
elif [[ $ALIGNER == $SNAPR ]]; then
    ALN_DIR_NAME=$SNAPR_ALIGN_DIR
    ALIGN_SCRIPT=$SNAPR_ALIGN_SCRIPT
    GENOME_INDEX=$SNAPR_GENOME_INDEX
    TRANSCRIPTOME_INDEX=$SNAPR_TRANSCRIPTOME_INDEX    
else
    echo ""
    echo "ERROR: Unrecognized aligner.  Exiting"
    exit 1
fi

export ALIGNER
export ALN_DIR_NAME
export ALIGN_SCRIPT
export GENOME_INDEX
export TRANSCRIPTOME_INDEX
############################################################

############################################################

#print out the parameters for logging:

echo ""
echo "Will attempt to perform analysis on samples (from "$SAMPLES_FILE"):"
print_sample_report $SAMPLES_FILE
echo ""
if [ "$CONTRAST_FILE" == "" ]; then
	echo "Will NOT perform differential analysis since no contrast file supplied."
else
	echo "Will attempt to perform the following contrasts (from "$CONTRAST_FILE"):"
	print_contrast_report $CONTRAST_FILE
fi
echo ""
echo "Project home directory: "$PROJECT_DIR
if [ $ALN -eq $NUM1 ]; then
	echo "Perform alignment with "$ALIGNER
	echo "Alignment will be performed against: "$ASSEMBLY
fi
echo ""
echo ""

############################################################


############################################################
#check if alignment is needed
# if yes, perform alignment

if [ $ALN -eq $NUM1 ]; then

    #call a python script that scans the sample directory, checks for the correct files,
    # and injects the proper parameters into the alignment shell script
    $PYTHON $PREPARE_ALIGN_SCRIPT || { echo "Something went wrong in preparing the alignment scripts.  Exiting"; exit 1; }

    echo "After examining project structure, will attempt to align on the following samples:"
    print_sample_report $VALID_SAMPLE_FILE

    #given the valid samples (determined by the python script), run the alignment
    # note that this is NOT done in parallel!  SNAPR and STAR are VERY memory intensive

    for sample in $( cut -f1 $VALID_SAMPLE_FILE ); do
		FLAG=0
		while [ $FLAG -eq 0 ]; do
			#check if other SNAPR or STAR processes are running first:
			if [ "$(pgrep $SNAPR_PROCESS)" == "" ] && [ "$(pgrep $STAR)" == "" ]; then
			        echo "Run alignment with script at: "
				echo $PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$sample'/'$sample$FORMATTED_ALIGN_SCRIPT_NAMETAG   
				date
				echo ""
                                chmod a+x $PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$sample'/'$sample$FORMATTED_ALIGN_SCRIPT_NAMETAG
				
				#kickoff the script and wait until completion:
			        if [ $TEST -eq $NUM0 ]; then
					$PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$sample'/'$sample$FORMATTED_ALIGN_SCRIPT_NAMETAG                
				else
					echo "...[Mock alignment]..."
					mkdir $PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$sample'/'$ALN_DIR_NAME
					touch $PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$sample'/'$ALN_DIR_NAME'/'$sample$SORTED_TAG$BAM_EXTENSION
				fi
				echo "Alignment on sample $sample completed at: "
				date
				FLAG=1
			else
				if [ $ATTEMPTS < $MAX_ALIGN_ATTEMPTS ]; then
					echo "Another snapr or star process is running.  Will sleep for $SLEEP_TIME seconds and try again. (Max attempts=$MAX_ALIGN_ATTEMPTS)"
					sleep $SLEEP_TIME
					let ATTEMPTS+=1
				else
					echo "Reached the maximum amount of attempts.  Exiting with task incomplete."
					exit 1
				fi
			fi
		done
	
    done

else

    echo "Skipping alignment based on input parameters (-noalign).  Locating BAM files..."
   
    #since we did not align here using STAR, etc., change the ALN_DIR_NAME to something generic:
    ALN_DIR_NAME="bam_file"

    #given bam files contained anywhere in PROJECT_HOME, construct the assumed project
    #hierarchy and create symbolic links to the bam files

    while read line; do
    	SAMPLE=$(echo $line | awk '{print $1}')

	#the name of the link (in our convention) which will link to the original bam file
	BASE_BAM_FILE=$SAMPLE$SORTED_TAG$BAM_EXTENSION
	BAM_FILE=$(find -L $PROJECT_DIR -type f -name $SAMPLE*$BAM_EXTENSION)
	if [ "$BAM_FILE" != "" ]; then
		SAMPLE_ALN_DIR=$PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$SAMPLE'/'$ALN_DIR_NAME
		mkdir -p $SAMPLE_ALN_DIR
		ln -s $BAM_FILE $SAMPLE_ALN_DIR'/'$BASE_BAM_FILE
#		echo $line >> $VALID_SAMPLE_FILE
		printf "%s\t%s\n" $(echo $SAMPLE) $(echo $line | awk '{print $2}') >> $VALID_SAMPLE_FILE
	else
		echo "Could not locate a properly named BAM file for sample "$SAMPLE
	fi
    done < $SAMPLES_FILE

    echo "Found BAM files for the following samples:"
    print_sample_report $VALID_SAMPLE_FILE
    
fi


############################################################

# check for the appropriate bam files and update the valid sample file accordingly:
$PYTHON $CHECK_BAM_SCRIPT || { echo "Error.  Exiting"; exit 1; }

############################################################

############################################################
# create or move the count files to prepare for differential analysis:

mkdir -p $COUNTS_DIR

if [ ! -d "$COUNTS_DIR" ]; then
    echo "Could not create the count directory (permissions?).  Exiting"
    exit 1
fi

if [ $TEST -eq $NUM0 ]; then

	#if used STAR for alignment or was directly provided with BAM files, get the read counts: 
	if [[ $ALIGNER == $STAR  ]] || [[ $ALN -eq $NUM0 ]] ; then

		#a tag for the temporary count file:
		TMP='.tmp'

    		#create the count files from the BAM files
    		for sample in $( cut -f1 $VALID_SAMPLE_FILE ); do
		    COUNT_FILE=$COUNTS_DIR'/'$sample$COUNTFILE_SUFFIX

    		    featureCounts -a $GTF \
				  -o $COUNT_FILE$TMP \
				  -t exon \
				  -g gene_name $PROJECT_DIR"/"$SAMPLE_DIR_PREFIX$sample"/"$ALN_DIR_NAME"/"$sample$SORTED_TAG$BAM_EXTENSION

		    Rscript \
			$PROCESS_COUNT_FILE_SCRIPT \
			$COUNT_FILE$TMP \
			$COUNT_FILE || { echo "Parsing of the raw count files failed.  Exiting."; exit 1; }
		    
                    rm $COUNT_FILE$TMP &

    		done
	elif [[ $ALIGNER == $SNAPR ]]; then
	    #move the count files
	    for sample in $( cut -f1 $VALID_SAMPLE_FILE ); do
	    	mv $PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$sample'/'$ALN_DIR_NAME'/'$sample$SORTED_TAG$COUNTFILE_SUFFIX $COUNTS_DIR'/'$sample$COUNTFILE_SUFFIX        
            done
	fi

#if test case (create dummy files in counts directory):
else 
	for sample in $( cut -f1 $VALID_SAMPLE_FILE ); do
		touch $COUNTS_DIR'/'$sample$COUNTFILE_SUFFIX
	done
fi

#with the count files moved, create a design matrix for DESeq:
$PYTHON $CREATE_DESIGN_MATRIX_SCRIPT || { echo "Failed in creating design matrix.  Exiting"; exit 1; }

############################################################


############################################################

#create a report directory to hold the report and the output analysis:
REPORT_DIR=$PROJECT_DIR'/'$REPORT_DIR
mkdir $REPORT_DIR

############################################################



############################################################
#get the normalized counts via DESEQ:

#first create an output directory for the deseq scripts that will be run:
DESEQ_RESULT_DIR=$REPORT_DIR'/'$DESEQ_RESULT_DIR
mkdir $DESEQ_RESULT_DIR

NORMALIZED_COUNTS_FILE=$DESEQ_RESULT_DIR'/'$NORMALIZED_COUNTS_FILE
export NORMALIZED_COUNTS_FILE

if [ $TEST -eq $NUM0 ]; then
	Rscript $NORMALIZED_COUNTS_SCRIPT $DESIGN_MTX_FILE $NORMALIZED_COUNTS_FILE
else
	echo "Perform mock normalized counts with DESeq."
fi

if [ -e "$CONTRAST_FILE" ]; then
	echo "Run differential expression with DESeq"
	while read contrast; do
    		conditionA=$(echo $contrast | awk '{print $1}')
    		conditionB=$(echo $contrast | awk '{print $2}')

    		if [ $TEST -eq $NUM0 ]; then
    			echo Rscript $DESEQ_SCRIPT $DESEQ_RESULT_DIR $DESIGN_MTX_FILE $DESEQ_OUTFILE_TAG $conditionA $conditionB $HEATMAP_FILE $HEATMAP_GENE_COUNT
    			Rscript $DESEQ_SCRIPT $DESEQ_RESULT_DIR $DESIGN_MTX_FILE $DESEQ_OUTFILE_TAG $conditionA $conditionB $HEATMAP_FILE $HEATMAP_GENE_COUNT
		else
			echo "Perform mock DESeq step on contrast between "$conditionA "and" $conditionB
    		fi
	done < $CONTRAST_FILE
else
	echo "Skipping differential analysis since no contrast file was specified."
fi
############################################################


############################################################
#Report creation:

#create the reports with RNA-seQC:
while read line; do
 	SAMPLE=$(echo $line | awk '{print $1}')
	SAMPLE_BAM=$PROJECT_DIR'/'$SAMPLE_DIR_PREFIX$SAMPLE'/'$ALN_DIR_NAME'/'$SAMPLE$SORTED_TAG$BAM_EXTENSION

        #create output directory
        SAMPLE_QC_DIR=$REPORT_DIR'/'$RNA_SEQC_DIR'/'$SAMPLE
        mkdir -p $SAMPLE_QC_DIR
 
    	if [ $TEST -eq $NUM0 ]; then
		java -jar $RNA_SEQC_JAR \
		-o $SAMPLE_QC_DIR \
		-r $GENOMEFASTA \
		-s "$SAMPLE|$SAMPLE_BAM|-" \
		-t $GTF_FOR_RNASEQC || { echo "Something failed on performing QC step.  Check the output for guidance."; }
	else
		echo "Perform mock QC analysis, etc. on "$SAMPLE
	fi
done < $VALID_SAMPLE_FILE

############################################################
# GSEA:

#create the GSEA directory:
GSEA_OUTPUT_DIR=$REPORT_DIR'/'$GSEA_OUTPUT_DIR
mkdir $GSEA_OUTPUT_DIR
export GSEA_OUTPUT_DIR

#create the formatted GSEA input files
GSEA_CLS_FILE=$GSEA_OUTPUT_DIR'/'$GSEA_CLS_FILE
GSEA_GCT_FILE=$GSEA_OUTPUT_DIR'/'$GSEA_GCT_FILE

Rscript $CREATE_GSEA_CLS_SCRIPT $VALID_SAMPLE_FILE $GSEA_CLS_FILE
Rscript $CREATE_GSEA_GCT_SCRIPT $NORMALIZED_COUNTS_FILE $VALID_SAMPLE_FILE $GSEA_GCT_FILE


if [ -e "$CONTRAST_FILE" ]; then
        while read contrast; do
                conditionA=$(echo $contrast | awk '{print $1}')
                conditionB=$(echo $contrast | awk '{print $2}')

                if [ $TEST -eq $NUM0 ]; then
			$RUN_GSEA_SCRIPT \
				$GSEA_JAR \
				$GSEA_ANALYSIS \
				$GSEA_GCT_FILE \
				$GSEA_CLS_FILE \
				$conditionA'_versus_'$conditionB \
				$DEFAULT_GMX_FILE \
				$NUM_GSEA_PERMUTATIONS \
				$conditionA'_versus_'$conditionB \
				$DEFAULT_CHIP_FILE \
				$GSEA_OUTPUT_DIR || { echo "Error occurred in running GSEA.  Exiting."; exit 1; }
                else
                    	echo "Perform mock GSEA step on contrast between "$conditionA "and" $conditionB
                fi
        done < $CONTRAST_FILE
else
    	echo "Skipping GSEA analysis since no contrast file was specified."
fi



############################################################


#Report creation:

#copy the necessary libraries to go with the html report:
cp -r $REPORT_TEMPLATE_LIBRARIES $REPORT_DIR

#run the injection script to create the report:
if [ $TEST -eq $NUM0 ]; then
	$PYTHON $CREATE_REPORT_SCRIPT
else
	echo "Perform mock creation of output report."
fi
############################################################


############################################################
#cleanup

#rm $VALID_SAMPLE_FILE


