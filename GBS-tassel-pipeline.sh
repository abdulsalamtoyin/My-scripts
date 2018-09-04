mkdir my-gbs-project
mkdir my-gbs-project/analysis
mkdir my-gbs-project/analysis/fastq
mkdir my-gbs-project/analysis/discovery
mkdir my-gbs-project/analysis/discovery/db
mkdir my-gbs-project/analysis/discovery/tagexport
mkdir my-gbs-project/analysis/discovery/production


#generate tags from fastq file
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -GBSSeqToTagDBPlugin -e HindIII \
-i fastq \
-db db/GBS_all.db \
-k keysv2_16.txt \
-c 10 \
-kmerLength 80 \
-minKmerL 20 \
-mnQS 20 \
-mxKmerNum 400000000 \
-batchSize 8 \
-endPlugin -runfork1 | tee GBSSeqToTagDB-GBS_all.log

#extract tags from db for alignment to reference
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -fork1 -TagExportToFastqPlugin \
-db db/GBS_all.db \
-o db/tagexport/tagsForAlign_all.fa.gz \
-c 1 \
-endPlugin -runfork1 | tee TagExportToFastq_all.log


#map tags to reference                   
bowtie2 -p 6 \
--very-sensitive-local \
-x /run/media/diersadmin/Soybean/GBS_analysis/Gmax2.0/Gmax2.0_bowtie_index/Gmax_2.0_bowtie \
-U db/tagexport/tagsForAlign_all.fa.gz \
-S db/tagexport/tagsForAlignFullvs_all.sam

#store tag physical map position in database
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -fork1 -SAMToGBSdbPlugin \
-i db/tagexport/tagsForAlignFullvs_all.sam \
-db db/GBS_all.db \
-aProp 0 \
-aLen 0 \
-endPlugin -runfork1 | tee SAMToGBSdb_all.log

#run SNP discovery plugin
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -DiscoverySNPCallerPluginV2 \
-db db/GBS_all.db \
-ref /run/media/diersadmin/Soybean/GBS_analysis/Gmax2.0/Gmax_275_v2.0_renamed.fa \
-mnLCov 0.05 \
-mnMAF 0.01 \
-deleteOldData true
-endPlugin -runfork1 | tee DiscoverySNPCaller_all.log

#check SNP quality parameters (Marco's samples)
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -SNPQualityProfilerPlugin \
-db db/GBS_all.db \
-taxa "taxa_marco.txt" \
-statFile "outputStatsMarco.txt" \
-endPlugin -runfork1 | tee SNPQualityProfilerMarco.log

#check SNP quality parameters (Eric's samples)
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -SNPQualityProfilerPlugin \
-db db/GBS_all.db \
-taxa "taxa_eric.txt" \
-statFile "outputStatsEric.txt" \
-deleteOldData true -endPlugin -runfork1 | tee SNPQualityProfilerEric.log

#calculate quality scores based on parameters estimated from previous plugin (used R script for this purpose, please read instructions in snp-quality.R script file). Runnin first for Marco's samples.
Rscript ../programs/snp-quality-marco.R

#calculate quality scores based on parameters estimated from previous plugin (used R script for this purpose, please read instructions in snp-quality.R script file). Now for Eric's samples.
Rscript ../programs/snp-quality-eric.R

#update database with quality scores (Marco's samples)
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -UpdateSNPPositionQualityPlugin \
-db db/GBS_all.db \
-qsFile myQsFileMarco.txt \
-endPlugin -runfork1 | tee UpdateSNPPositionQualityMarco.log

#update database with quality scores (Eric's samples)
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -UpdateSNPPositionQualityPlugin \
-db db/GBS_all.db \
-qsFile myQsFileEric.txt \
-endPlugin -runfork1 | tee UpdateSNPPositionQualityEric.log


#run production pipeline (Marco's samples)
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -ProductionSNPCallerPluginV2 \
-db db/GBS_all.db \
-e HindIII \
-i fastq \
-k keysv2_production_marco.txt \
-kmerLength 80 \
-minPosQS 3 \
-o discovery/production/filtered-snps-production-pipeline-Marco.vcf \
-endPlugin -runfork1 | tee ProductionSNPCallerV2Marco.log

#run production pipeline (Eric's samples)
perl '/run/media/diersadmin/Soybean/GBS_analysis/programs/tassel5/run_pipeline.pl' -Xmx50g -fork1 -ProductionSNPCallerPluginV2 \
-db db/GBS_all.db \
-e HindIII \
-i fastq \
-k keysv2_production_eric.txt \
-kmerLength 80 \
-minPosQS 3 \
-o discovery/production/filtered-snps-production-pipeline-Ericall0.05.vcf \
-endPlugin -runfork1 | tee ProductionSNPCallerV2Eric.log