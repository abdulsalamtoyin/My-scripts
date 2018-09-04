fastq=/data/biodata/Zmays/StrigaMaizeRNA-seq/
Dir=/syn02iscsi/toyin/STRIGA-ANALYSIS/trim_reads

for i in $(ls $fastq/*_001.fastq.gz | rev | cut -c 16- | rev | uniq)
do
java -jar /opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 30 ${i}R1_001.fastq.gz ${i}R2_001.fastq.gz ${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" ILLUMINACLIP:/opt/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:14
mv ${i}"R1.trimmed_PE.fastq.gz" ${i}"R1.trimmed_SE.fastq.gz" ${i}"R2.trimmed_PE.fastq.gz" ${i}"R2.trimmed_SE.fastq.gz" $fastq
done
