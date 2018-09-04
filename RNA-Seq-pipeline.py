#!/usr/bin/env python
# -*- coding: utf-8 -*-
__metaclass__ = type
# Author: Li Guochao

import os
from time import strftime as time
from optparse import OptionParser
from operator import itemgetter
from math import log

parser = OptionParser(usage="%prog [-g] Reference.gtf [-a] Reference.fa [-t] num_of_threads [-e] run enrichment analysis on the Internet or not (default is not) [-i] control,fastq1,fastq2:treat,fastq1,fastq2 [-o] output_path", version="Version: 0.8")

parser.add_option("-g", 
                  "--gtf",
                  dest      = "reference_gtf", 
                  type      = "string",
                  help      = "The annotation GTF file of reference. (eg: /leofs/sunyl_group/yaolsh/ref/hg19/hg.gtf or ./hg.gtf)") 

parser.add_option("-a", 
                  "--fa",
                  dest      = "reference_fa", 
                  type      = "string",
                  help      = "The reference name in .fa format. (eg: /leofs/sunyl_group/yaolsh/ref/hg19/hg.fa or ./hg.fa)") 

parser.add_option("-t", 
                  "--thread",
                  dest      = "thread", 
                  type      = "int",
                  default   = 8,
                  help      = "The number of threads, default is 8.") 

parser.add_option("-e", 
                  "--enrichment", 
                  dest      = "enrichment", 
                  action    = "store_true",
                  default   = False,
                  help      = "Tell Analyzer_RNA-Seq to run GO and KEGG enrichment analysis on the Internet or not. Default is False. If '-e' is given, this step will run.") 

parser.add_option("-i", 
                  "--input",
                  "--fq",
                  dest      = "input_groups_and_filenames", 
                  type      = "string",
                  help      = "The names of fastq files for analysis. (eg: 231-2,/leofs/sunyl_group/ligch/YuHui/data/RNA-Seq/231-2/231-2_1.fastq,/leofs/sunyl_group/ligch/YuHui/data/RNA-Seq/231-2/231-2_2.fastq:231-1,/leofs/sunyl_group/ligch/YuHui/data/RNA-Seq/231-1/231-1_1.fastq,/leofs/sunyl_group/ligch/YuHui/data/RNA-Seq/231-1/231-1_2.fastq or NCSC,./231-2/231-2_1.fastq,./231-2/231-2_1.fastq:CSC,./231-1/231-1_1.fastq,./231-1/231-1_2.fastq)") 

parser.add_option("-o", 
                  "--output_path", 
                  dest      = "output_path", 
                  type      = "string",
                  help      = "The path of output. (eg: /leofs/sunyl_group/ligch/YuHui/analysis/RNA-Seq or ./)") 

(options, args) = parser.parse_args() 


### Preparation
if options.output_path and options.reference_gtf and options.reference_fa and options.input_groups_and_filenames and options.thread:

    gtf             = os.path.abspath(options.reference_gtf)
    fa              = os.path.abspath(options.reference_fa)
    thread          = options.thread
    enrichment      = options.enrichment
    groups_and_fqs  = options.input_groups_and_filenames.split(":")
    output_path     = os.path.abspath(options.output_path)
    
    group_1         = groups_and_fqs[0].split(",")[0]
    fq_group_1      = [os.path.abspath(i) for i in groups_and_fqs[0].split(",")[1:]]
    group_2         = groups_and_fqs[1].split(",")[0]
    fq_group_2      = [os.path.abspath(i) for i in groups_and_fqs[1].split(",")[1:]]

    os.chdir(output_path)


    print 
    print "=" * 40 + " Analyzer_RNA-Seq starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 40 + "\n"
    with open("Run_time", "w") as run_time:
        run_time.write("Analyzer_RNA-Seq starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n\n")


### Step 1. QC
    print
    print "=" * 40 + " 1. QC starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 40 + "\n"
    with open("Run_time", "a") as run_time:
        run_time.write("1. QC starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n")

    def qc(qc_path, group, fq_1, fq_2):
        qc_group_path = os.path.join(qc_path, group)
        os.mkdir(qc_group_path)
        qc_group_path_prefix = os.path.join(qc_path, group, group)
        qc_cmdline = """piplineForQC preprocessing \
                -f %(fq_1)s \
                -r %(fq_2)s \
                -o %(qc_group_path_prefix)s \
                -s 1,2,3,4,6
        """ % (
                {
                    "fq_1"                  : fq_1,
                    "fq_2"                  : fq_2,
                    "qc_group_path_prefix"  : qc_group_path_prefix,
                }
              )
        os.system(qc_cmdline)

    def fastqc(qc_path, group, fq_1, fq_2, before_or_after_QC):
        fastqc_before_or_after_QC_path = os.path.join(qc_path, group, "fastqc", before_or_after_QC)
        os.mkdir(fastqc_before_or_after_QC_path)
        fastqc_cmdline = """fastqc \
                -t %(thread)d \
                -o %(fastqc_before_or_after_QC_path)s \
                %(fq_1)s \
                %(fq_2)s
        """ % (
                {
                    "thread"                            : thread,
                    "fq_1"                              : fq_1,
                    "fq_2"                              : fq_2,
                    "fastqc_before_or_after_QC_path"    : fastqc_before_or_after_QC_path,
                }
              )
        os.system(fastqc_cmdline)

    qc_path = os.path.join(output_path, "1_QC")
    os.mkdir(qc_path)
    qc(qc_path, group_1, fq_group_1[0], fq_group_1[1])
    qc(qc_path, group_2, fq_group_2[0], fq_group_2[1])

    os.mkdir( os.path.join(qc_path, group_1, "fastqc") )
    os.mkdir( os.path.join(qc_path, group_2, "fastqc") )

    fq_group_1_after_QC = []
    fq_group_1_after_QC.append(os.path.join(qc_path, group_1, group_1 + "_1.fastq"))
    fq_group_1_after_QC.append(os.path.join(qc_path, group_1, group_1 + "_2.fastq"))
    fq_group_2_after_QC = []
    fq_group_2_after_QC.append(os.path.join(qc_path, group_2, group_2 + "_1.fastq"))
    fq_group_2_after_QC.append(os.path.join(qc_path, group_2, group_2 + "_2.fastq"))

    fastqc(qc_path, group_1, fq_group_1[0]         , fq_group_1[1]         , "before_QC")
    fastqc(qc_path, group_1, fq_group_1_after_QC[0], fq_group_1_after_QC[1], "after_QC" )
    fastqc(qc_path, group_2, fq_group_2[0]         , fq_group_2[1]         , "before_QC")
    fastqc(qc_path, group_2, fq_group_2_after_QC[0], fq_group_2_after_QC[1], "after_QC" )


### Step 2. tophat
    print
    print "=" * 40 + " 2. Tophat starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 40 + "\n"
    with open("Run_time", "a") as run_time:
        run_time.write("2. Tophat starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n")

    def tophat(output_path, gtf, fa, group, fq_1, fq_2, thread):
        tophat_group_path = os.path.join(output_path, "2_tophat", group)
        os.mkdir(tophat_group_path)
        tophat_cmdline = """tophat \
            -p %(thread)d \
            -G %(gtf)s \
            -o %(tophat_group_path)s \
            %(fa)s \
            %(fq_1)s \
            %(fq_2)s
        """ % (
                {
                    "thread"            : thread,
                    "gtf"               : gtf,
                    "tophat_group_path" : tophat_group_path,
                    "fa"                : fa,
                    "fq_1"              : fq_1,
                    "fq_2"              : fq_2,
                }
              )
        os.system(tophat_cmdline)

    os.mkdir( os.path.join(output_path, "2_tophat") )
    tophat(output_path, gtf, fa, group_1, fq_group_1_after_QC[0], fq_group_1_after_QC[1], thread)
    tophat(output_path, gtf, fa, group_2, fq_group_2_after_QC[0], fq_group_2_after_QC[1], thread)


### Step 3. cufflinks
    print
    print "=" * 40 + " 3. Cufflinks starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 40 + "\n"
    with open("Run_time", "a") as run_time:
        run_time.write("3. Cufflinks starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n")

    def cufflinks(output_path, group, thread):
        cufflinks_group_path =  os.path.join(output_path, "3-5_cufflinks", "3_cufflinks", group)
        os.mkdir(cufflinks_group_path)
        cufflinks_cmdline = """cufflinks \
            -p %(thread)d \
            -o %(cufflinks_group_path)s \
            %(tophat_bam_path)s
        """ % (
                {
                    "thread"                : thread,
                    "cufflinks_group_path"  : cufflinks_group_path,
                    "tophat_bam_path"       : os.path.join(output_path, "2_tophat", group, "accepted_hits.bam"),
                }
              )
        os.system(cufflinks_cmdline)

    os.mkdir( os.path.join(output_path, "3-5_cufflinks") )
    os.mkdir( os.path.join(output_path, "3-5_cufflinks", "3_cufflinks") )
    cufflinks(output_path, group_1, thread)
    cufflinks(output_path, group_2, thread)


### Step 4. cuffmerge
    print
    print "=" * 40 + " 4. Cuffmerge starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 40 + "\n"
    with open("Run_time", "a") as run_time:
        run_time.write("4. Cuffmerge starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n")

    def cuffmerge(cuffmerge_output_path, gtf, fa, thread, assemblies):
        cuffmerge_cmdline = """cuffmerge \
            -o %(cuffmerge_output_path)s \
            -g %(gtf)s \
            -s %(fa)s \
            -p %(thread)d \
            %(assemblies.txt)s
        """ % (
                {
                    "cuffmerge_output_path" : cuffmerge_output_path,
                    "gtf"                   : gtf,
                    "fa"                    : fa,
                    "thread"                : thread,
                    "assemblies.txt"        : assemblies,
                }
              )
        os.system(cuffmerge_cmdline)

    cuffmerge_output_path = os.path.join(output_path, "3-5_cufflinks", "4_cuffmerge") 
    os.mkdir(cuffmerge_output_path)
    with open(os.path.join(cuffmerge_output_path, "assemblies.txt"), "w") as cuffmerge_assemblies:
        out_line_1 = os.path.join(output_path, "3-5_cufflinks", "3_cufflinks", group_1, "transcripts.gtf")
        out_line_2 = os.path.join(output_path, "3-5_cufflinks", "3_cufflinks", group_2, "transcripts.gtf")
        cuffmerge_assemblies.write(out_line_1 + "\n" + out_line_2 + "\n")
    assemblies = os.path.join(cuffmerge_output_path, "assemblies.txt")
    cuffmerge(cuffmerge_output_path, gtf, fa, thread, assemblies)


### Step 5. cuffdiff
    print
    print "=" * 40 + " 5. Cuffdiff starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 40 + "\n"
    with open("Run_time", "a") as run_time:
        run_time.write("5. Cuffdiff starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n")

    def cuffdiff(cuffdiff_output_path, thread, merged_gtf, tophat_group_1_bam, tophat_group_2_bam):
        cuffdiff_cmdline = """cuffdiff \
            -o %(cuffdiff_output_path)s \
            -b %(fa)s \
            -p %(thread)d \
            -L N,T \
            -u %(merged_gtf)s \
            %(tophat_group_1_bam)s \
            %(tophat_group_2_bam)s 
        """ % (
                {
                    "cuffdiff_output_path"  : cuffdiff_output_path,
                    "fa"                    : fa,
                    "thread"                : thread,
                    "merged_gtf"            : merged_gtf,
                    "tophat_group_1_bam"    : tophat_group_1_bam,
                    "tophat_group_2_bam"    : tophat_group_2_bam,
                }
              )
        os.system(cuffdiff_cmdline)

    cuffdiff_output_path =  os.path.join(output_path, "3-5_cufflinks", "5_cuffdiff")
    os.mkdir(cuffdiff_output_path)
    merged_gtf = os.path.join(cuffmerge_output_path, "merged.gtf")
    tophat_group_1_bam = os.path.join(output_path, "2_tophat", group_1, "accepted_hits.bam")
    tophat_group_2_bam = os.path.join(output_path, "2_tophat", group_2, "accepted_hits.bam")
    cuffdiff(cuffdiff_output_path, thread, merged_gtf, tophat_group_1_bam, tophat_group_2_bam)



### Step 6. Find genes with different expression level (up-regulated or down-regulated) and filter them by p or q value. Output results in text format. In addition, output another list of all the diff-expression genes for valcano plot. Draw density, scatter and volcano plot and heatmap based on the results above. Finally, various heatmaps are drawn by different methods of clustering and measurements of distance. 
    print
    print "=" * 35 + " 6. Find_diff_exp_genes starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 35 + "\n"
    with open("Run_time", "a") as run_time:
        run_time.write("6. Find_diff_exp_genes starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n")

    find_diff_exp_gene_output_path =  os.path.join(output_path, "6_find_diff_exp_gene")
    os.mkdir(find_diff_exp_gene_output_path)

    cuffdiff_gene_exp_diff            = open(os.path.join(cuffdiff_output_path          , "gene_exp.diff"               ), "r")
    diffgene_p_up                     = open(os.path.join(find_diff_exp_gene_output_path, "diffgene_p_up"               ), "w")
    diffgene_p_down                   = open(os.path.join(find_diff_exp_gene_output_path, "diffgene_p_down"             ), "w")
    diffgene_q_up                     = open(os.path.join(find_diff_exp_gene_output_path, "diffgene_q_up"               ), "w")
    diffgene_q_down                   = open(os.path.join(find_diff_exp_gene_output_path, "diffgene_q_down"             ), "w")
    gene_log2FC_p_q_for_volcano       = open(os.path.join(find_diff_exp_gene_output_path, "gene_log2FC_p_q_for_volcano" ), "w")

    diffgene_path_p_for_heatmap       =      os.path.join(find_diff_exp_gene_output_path, "diffgene_p_for_heatmap"      )    
    diffgene_path_q_for_heatmap       =      os.path.join(find_diff_exp_gene_output_path, "diffgene_q_for_heatmap"      )    
    diffgene_p_for_heatmap            = open(diffgene_path_p_for_heatmap                                                 , "w")
    diffgene_q_for_heatmap            = open(diffgene_path_q_for_heatmap                                                 , "w")

    diffgene_genenamelist_path_p_up   =      os.path.join(find_diff_exp_gene_output_path, "diffgene_genenamelist_p_up"  )
    diffgene_genenamelist_path_p_down =      os.path.join(find_diff_exp_gene_output_path, "diffgene_genenamelist_p_down")
    diffgene_genenamelist_path_q_up   =      os.path.join(find_diff_exp_gene_output_path, "diffgene_genenamelist_q_up"  )
    diffgene_genenamelist_path_q_down =      os.path.join(find_diff_exp_gene_output_path, "diffgene_genenamelist_q_down")
    diffgene_genenamelist_p_up        = open(diffgene_genenamelist_path_p_up                                             , "w")
    diffgene_genenamelist_p_down      = open(diffgene_genenamelist_path_p_down                                           , "w")
    diffgene_genenamelist_q_up        = open(diffgene_genenamelist_path_q_up                                             , "w")
    diffgene_genenamelist_q_down      = open(diffgene_genenamelist_path_q_down                                           , "w")

    cuffdiff_gene_exp_diff_lines = cuffdiff_gene_exp_diff.readlines()
    cuffdiff_gene_exp_diff_title = cuffdiff_gene_exp_diff_lines[0].strip().split("\t")

    diffgene_p_up.write("\t".join(cuffdiff_gene_exp_diff_title) + "\n")
    diffgene_p_down.write("\t".join(cuffdiff_gene_exp_diff_title) + "\n")
    diffgene_p_for_heatmap.write("\t".join([cuffdiff_gene_exp_diff_title[2]] + ["log2(N_FPKM)"] + ["log2(T_FPKM)"]) + "\n")
    diffgene_q_up.write("\t".join(cuffdiff_gene_exp_diff_title) + "\n")
    diffgene_q_down.write("\t".join(cuffdiff_gene_exp_diff_title) + "\n")
    diffgene_q_for_heatmap.write("\t".join([cuffdiff_gene_exp_diff_title[2]] + ["log2(N_FPKM)"] + ["log2(T_FPKM)"]) + "\n")
    gene_log2FC_p_q_for_volcano.write("\t".join([cuffdiff_gene_exp_diff_title[2]] \
                                            + [cuffdiff_gene_exp_diff_title[9]] \
                                            + [cuffdiff_gene_exp_diff_title[11]] \
                                            + [cuffdiff_gene_exp_diff_title[12]]) \
                                            + "\n")

    diffgene_gene_line = {}
    for line in cuffdiff_gene_exp_diff_lines[1:]:
        line = line.strip().split("\t")
        if line[2] == "-":                      # Go to the next line if genename in this line is "-" (none).
            continue

        if line[2] not in diffgene_gene_line:
            diffgene_gene_line[line[2]] = line 
        else:
            diffgene_gene_line[line[2]] += line

    diffgene_p_list = []
    diffgene_q_list = []
    for line in diffgene_gene_line.values():
        if len(line) > 14:                      
            continue                            # There are some lines whose genenames are same and locuses are different. However, their expression patterns might be completely contrast, in other words, a gene might be both up-regulated and down-regulated in the same sample. In this version, Analyzer_RNA-Seq just simply deletes these lines. A better solution might come in the future.
        else:
            if float(line[11]) < 0.05:
                if float(line[9]) > 1:
                    diffgene_p_up.write("\t".join(line) + "\n")
                    diffgene_genenamelist_p_up.write(line[2] + "\n")
                    if float(line[7]) > 0.1 and float(line[8]) > 0.1:
                        diffgene_p_list.append([line[2]] + [str(log(float(line[7]),2))] + [str(log(float(line[8]),2))] + [float(line[11])])
                elif float(line[9]) < -1:
                    diffgene_p_down.write("\t".join(line) + "\n")
                    diffgene_genenamelist_p_down.write(line[2] + "\n")
                    if float(line[7]) > 0.1 and float(line[8]) > 0.1:
                        diffgene_p_list.append([line[2]] + [str(log(float(line[7]),2))] + [str(log(float(line[8]),2))] + [float(line[11])])
                        
            if float(line[12]) < 0.05:
                if float(line[9]) > 1:
                    diffgene_q_up.write("\t".join(line) + "\n")
                    diffgene_genenamelist_q_up.write(line[2] + "\n")
                    if float(line[7]) > 0.1 and float(line[8]) > 0.1:
                        diffgene_q_list.append([line[2]] + [str(log(float(line[7]),2))] + [str(log(float(line[8]),2))] + [float(line[11])])
                elif float(line[9]) < -1:
                    diffgene_q_down.write("\t".join(line) + "\n")
                    diffgene_genenamelist_q_down.write(line[2] + "\n")
                    if float(line[7]) > 0.1 and float(line[8]) > 0.1:
                        diffgene_q_list.append([line[2]] + [str(log(float(line[7]),2))] + [str(log(float(line[8]),2))] + [float(line[11])])

            gene_log2FC_p_q_for_volcano.write("\t".join([line[2]] + [line[9]] + [line[11]] + [line[12]]) + "\n")

    diffgene_p_sort = sorted(diffgene_p_list, key = itemgetter(3))
    diffgene_q_sort = sorted(diffgene_q_list, key = itemgetter(3))

    for line in diffgene_p_sort:
        diffgene_p_for_heatmap.write("\t".join(line[:3]) + "\n")

    for line in diffgene_q_sort:
        diffgene_q_for_heatmap.write("\t".join(line[:3]) + "\n")

    cuffdiff_gene_exp_diff.close()
    diffgene_p_up.close()
    diffgene_p_down.close()
    diffgene_p_for_heatmap.close()
    diffgene_q_up.close()
    diffgene_q_down.close()
    diffgene_q_for_heatmap.close()
    gene_log2FC_p_q_for_volcano.close()
    diffgene_genenamelist_p_up.close()
    diffgene_genenamelist_p_down.close()
    diffgene_genenamelist_q_up.close()
    diffgene_genenamelist_q_down.close()

    r_script_summary_density_scatter = """
# set work directory
setwd("%(find_diff_exp_gene_output_path)s")
library(cummeRbund)
cuff_data <- readCufflinks("%(cuffdiff_output_path)s")
sink("summary.txt")
print(cuff_data)
sink()
tiff(file = "Density.tif", res = 600, width = 3260, height = 2000, compression = "lzw")
csDensity(genes(cuff_data))
dev.off()
tiff(file = "Scatter.tif", res = 600, width = 2000, height = 2000, compression = "lzw")
csScatter(genes(cuff_data), "N", "T")
dev.off()
    """ % ( 
            {
            "cuffdiff_output_path"              : cuffdiff_output_path,
            "find_diff_exp_gene_output_path"    : find_diff_exp_gene_output_path
            }
        )

    r_script_summary_density_scatter_path_filename = os.path.join(find_diff_exp_gene_output_path, "r_script_summary_density_scatter")
    with open(r_script_summary_density_scatter_path_filename, "w") as r_script_summary_density_scatter_filehand:
        r_script_summary_density_scatter_filehand.write(r_script_summary_density_scatter)
    os.system("Rscript " + r_script_summary_density_scatter_path_filename)
    os.system("rm " + r_script_summary_density_scatter_path_filename)

    r_script_volcano = """
# set work directory
setwd("%(find_diff_exp_gene_output_path)s")
library(ggplot2)
volcano_by_p_and_q<-read.table("gene_log2FC_p_q_for_volcano",head = T, row.names = 1)
# Data preparation.
df<-data.frame(volcano_by_p_and_q)
log_FC<-c(volcano_by_p_and_q[,1])
p<-c(volcano_by_p_and_q[,2])
q<-c(volcano_by_p_and_q[,3])
# Draw volcano plot of diffgene by p_value.
df.G<-subset(df,log_FC <= -1 & p < 0.05)
df.G<-cbind(df.G,rep(1,nrow(df.G)))
colnames(df.G)[4]<-"Color"
df.B<-subset(df,log_FC > -1 & log_FC < 1 | p > 0.05)
df.B<-cbind(df.B,rep(2,nrow(df.B)))
colnames(df.B)[4]<-"Color"
df.R<-subset(df,log_FC >= 1 & p < 0.05)
df.R<-cbind(df.R,rep(3,nrow(df.R)))
colnames(df.R)[4]<-"Color"
df.t<-rbind(df.G,df.B,df.R)
df.t$Color <- as.factor(df.t$Color)
tiff(file = paste0("Volcano_p.tif"), res = 600, width = 3236, height = 2000, compression = "lzw")
ggplot(data = df.t, aes(x = df.t[,1], y = -log10(df.t[,2]), color= Color )) + 
    geom_point(alpha = 0.5, size = 2) + 
    theme(legend.position = "none") + 
    xlim(c(-5, 5)) + 
    ylim(c(0, 20)) + 
    scale_color_manual(values = c("green", "grey", "red")) + 
    labs(x=expression(log[2](Fold_Change)), y=expression(-log[10](p_value))) + 
    theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=8)) + 
    theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=8))
dev.off()
# Draw volcano plot of diffgene by q_value.
df.G<-subset(df,log_FC <= -1 & q < 0.05)
df.G<-cbind(df.G,rep(1,nrow(df.G)))
colnames(df.G)[4]<-"Color"
df.B<-subset(df,log_FC > -1 & log_FC < 1 | q > 0.05)
df.B<-cbind(df.B,rep(2,nrow(df.B)))
colnames(df.B)[4]<-"Color"
df.R<-subset(df,log_FC >= 1 & q < 0.05)
df.R<-cbind(df.R,rep(3,nrow(df.R)))
colnames(df.R)[4]<-"Color"
df.t<-rbind(df.G,df.B,df.R)
df.t$Color <- as.factor(df.t$Color)
tiff(file = paste0("Volcano_q.tif"), res = 600, width = 3236, height = 2000, compression = "lzw")
ggplot(data = df.t, aes(x = df.t[,1], y = -log10(df.t[,3]), color= Color )) + 
    geom_point(alpha = 0.5, size = 1) + 
    theme( legend.position = "none") + 
    xlim(c(-5, 5)) + 
    ylim(c(0, 20)) + 
    scale_color_manual(values = c("green", "grey", "red")) + 
    labs(x=expression(log[2](Fold_Change)), y=expression(-log[10](q_value))) + 
    theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=8)) + 
    theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=8))
dev.off()
    """ % ( 
            {
                "find_diff_exp_gene_output_path"    : find_diff_exp_gene_output_path
            }
        )
        
    r_script_volcano_path_filename = os.path.join(find_diff_exp_gene_output_path, "r_script_volcano")
    with open(r_script_volcano_path_filename, "w") as r_script_volcano_filehand:
        r_script_volcano_filehand.write(r_script_volcano)
    os.system("Rscript " + r_script_volcano_path_filename)
    os.system("rm " + r_script_volcano_path_filename)
    os.system("rm " + os.path.join(find_diff_exp_gene_output_path , "gene_log2FC_p_q_for_volcano"))

    r_script_heatmap = """
# set work directory
work_directory <- "%(find_diff_exp_gene_output_path)s"
setwd(work_directory)
pheatmap_directory <- paste0(work_directory, "/heatmap")
dir.create(pheatmap_directory)
pheatmap_by_p_or_q <- function(data_directory, p_or_q) {
    pheatmap_directory_p_or_q <- paste0(work_directory, "/heatmap/", p_or_q)
    dir.create(pheatmap_directory_p_or_q)
    setwd(pheatmap_directory_p_or_q)
    diffgene_for_heatmap<-read.table(data_directory, head =TRUE, row.names = 1)
    names(diffgene_for_heatmap) <- c("log2(N_FPKM)" , "log2(T_FPKM)")
    num_of_rows <- nrow(diffgene_for_heatmap)
    if ( 200 <= num_of_rows ) {
        x = 0.9
        y = x - 0.1
        z = head(diffgene_for_heatmap, n = 200L)
    }
    if ( 25 < num_of_rows & num_of_rows < 200 ) {
        x = -0.009415 * num_of_rows + 2.608537
        y = x - 0.1
        z = diffgene_for_heatmap
    }
    if ( num_of_rows <= 25 ) {
        x = 2.5
        y = x - 0.1
        z = diffgene_for_heatmap
    }
    library(pheatmap)
    for (m in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
        for (n in c("correlation","euclidean", "maximum", "manhattan", "canberra", "minkowski")){
            tiff(file = paste0("heatmap_", m, "_", n, ".tif"), res = 600, width = 1800, height = 2600, compression = "lzw")
            pheatmap(z, color=colorRampPalette(c("white","yellow","red"))(299), show_rownames = TRUE, fontsize = 6, fontsize_col = 5, fontsize_row = y, cellwidth = 10, cellheight = x, clustering_method = m, clustering_distance_cols = n, clustering_distance_rows = n, main = paste0("Top 200 (by ", p_or_q," value) genes with different expression\nlevels (represnted by log2(FPKM)) between N and T"))
            dev.off()
        }
    }
}
pheatmap_by_p_or_q("%(diffgene_path_p_for_heatmap)s","p")
pheatmap_by_p_or_q("%(diffgene_path_q_for_heatmap)s","q")
    """ % ( 
            {
                "find_diff_exp_gene_output_path"    : find_diff_exp_gene_output_path,
                "diffgene_path_p_for_heatmap"       : diffgene_path_p_for_heatmap,
                "diffgene_path_q_for_heatmap"       : diffgene_path_q_for_heatmap
            }
        )
        
    r_script_heatmap_path_filename = os.path.join(find_diff_exp_gene_output_path, "r_script_heatmap")
    with open(r_script_heatmap_path_filename, "w") as r_script_heatmap_filehand:
        r_script_heatmap_filehand.write(r_script_heatmap)
    os.system("Rscript " + r_script_heatmap_path_filename)
    os.system("rm " + r_script_heatmap_path_filename)


### Step 7. GO and KEGG enrichmant analysis by R. If '-e' is not given to Analyzer_RNA-Seq, this step will be passed.
    if enrichment == True:
        print
        print "=" * 30 + " 7. GO and KEGG  enrichmant analysis (by R) starts at " + time('%Y-%m-%d %H:%M:%S') + " " + "=" * 30 + "\n"
        with open("Run_time", "a") as run_time:
            run_time.write("7. GO_and_KEGG_enrichment_analysis starts at " + time('%Y-%m-%d %H:%M:%S') + ".\n\n")

        GO_and_KEGG_enrichment_output_path          =  os.path.join(output_path, "7_GO_and_KEGG_enrichment_analysis"      )
        GO_and_KEGG_enrichment_output_path_p_up     =  os.path.join(GO_and_KEGG_enrichment_output_path, "p_up"   )
        GO_and_KEGG_enrichment_output_path_p_down   =  os.path.join(GO_and_KEGG_enrichment_output_path, "p_down" )
        GO_and_KEGG_enrichment_output_path_q_up     =  os.path.join(GO_and_KEGG_enrichment_output_path, "q_up"   )
        GO_and_KEGG_enrichment_output_path_q_down   =  os.path.join(GO_and_KEGG_enrichment_output_path, "q_down" )
        os.mkdir(GO_and_KEGG_enrichment_output_path)
        os.mkdir(GO_and_KEGG_enrichment_output_path_p_up)
        os.mkdir(GO_and_KEGG_enrichment_output_path_p_down)
        os.mkdir(GO_and_KEGG_enrichment_output_path_q_up)
        os.mkdir(GO_and_KEGG_enrichment_output_path_q_down)

        print GO_and_KEGG_enrichment_output_path_p_up
        print diffgene_genenamelist_path_p_up

        r_script_GO_and_KEGG_enrichment = """
# Use "GOstats" package to do GO or KEGG enrichment analysis of given gene symbols.
# library necessary packages
library("org.Hs.eg.db")
library("GSEABase")
library("GOstats")
library(Category)
library("pathview")
GO_and_KEGG_enrichment <- function(work_directory, data_directory) {
    # setwd and get current_wd
    setwd(work_directory)
    wd <- getwd()
    
    # read a gene name list for analysis
    genes <- read.table(data_directory)
    genes <- as.vector(unlist(genes))
    
    
    
    # GO and KEGG enrichment analysis
    
    # 1.mkdir
    GO_result_directory <- paste0(wd, "/GO_enrichment_result")
    dir.create(GO_result_directory)
    KEGG_result_directory <- paste0(wd, "/KEGG_enrichment_result")
    dir.create(KEGG_result_directory)
    # 2.GO enrichment analysis
    # (1) cd
    setwd(GO_result_directory)
    
    # (2) analysis
    goAnn <- get("org.Hs.egGO")
    universe <- Lkeys(goAnn)
    entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
    entrezIDs <- as.character(entrezIDs)
    for (m in c("BP","CC","MF")) {
        params <- new("GOHyperGParams", geneIds=entrezIDs, universeGeneIds=universe, annotation="org.Hs.eg.db", ontology = m, pvalueCutoff=0.05, conditional=FALSE, testDirection="over")
        over <- hyperGTest(params)
        glist <- geneIdsByCategory(over)
        glist <- sapply(glist, function(.ids) {
        .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
        .sym[is.na(.sym)] <- .ids[is.na(.sym)]
        paste(.sym, collapse=";")
        })
        gene_ontology <- summary(over)
        gene_ontology$Symbols <- glist[as.character(unlist(gene_ontology[1]))]
        write.table(gene_ontology,paste0("GO_enrichment_result_", m, ".txt"),sep = "\t",row.names = FALSE)
    }
    
    # (3) combination of 3 ontologys
    bp <- read.table(paste0(work_directory, "/GO_enrichment_result/","GO_enrichment_result_BP.txt"), sep = "\t", header = TRUE)
    bp_ParentTerm <- rep("Biological Process", length(bp$Term))
    go_enrichment_analysis_reuslt_BP <- cbind(bp[1:6],bp_ParentTerm,bp[7:8])
    names(go_enrichment_analysis_reuslt_BP)<-c("GO_ID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Sub_Term","Symbols")
    cc <- read.table(paste0(work_directory, "/GO_enrichment_result/","GO_enrichment_result_CC.txt"), sep = "\t", header = TRUE)
    cc_ParentTerm <- rep("Cellular Components", length(cc$Term))
    go_enrichment_analysis_reuslt_CC <- cbind(cc[1:6],cc_ParentTerm,cc[7:8])
    names(go_enrichment_analysis_reuslt_CC)<-c("GO_ID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Sub_Term","Symbols")
    mf <- read.table(paste0(work_directory, "/GO_enrichment_result/","GO_enrichment_result_MF.txt"), sep = "\t", header = TRUE)
    mf_ParentTerm <- rep("Molecullar Function", length(mf$Term))
    go_enrichment_analysis_reuslt_MF <- cbind(mf[1:6],mf_ParentTerm,mf[7:8])
    names(go_enrichment_analysis_reuslt_MF)<-c("GO_ID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Sub_Term","Symbols")
    go_enrichment_analysis_reuslt_all <- rbind(go_enrichment_analysis_reuslt_BP,go_enrichment_analysis_reuslt_CC,go_enrichment_analysis_reuslt_MF)
    write.table(go_enrichment_analysis_reuslt_all, paste0("GO_enrichment_result_ALL.txt"), sep = "\t", row.names = FALSE)
    
    
    
    # 3. KEGG enrichment analysis
    # (1) cd
    setwd(KEGG_result_directory)
    
    # (2) analysis
    keggAnn <- get("org.Hs.egPATH")
    universe <- Lkeys(keggAnn)
    params <- new("KEGGHyperGParams", geneIds=entrezIDs, universeGeneIds=universe, annotation="org.Hs.eg.db", categoryName="KEGG", pvalueCutoff=0.05,testDirection="over")
    over <- hyperGTest(params)
    kegg <- summary(over)
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
        .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
        .sym[is.na(.sym)] <- .ids[is.na(.sym)]
        paste(.sym, collapse=";")
    })
    kegg$Symbols <- glist[as.character(kegg$KEGGID)]
    gIds <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
    gEns <- unlist(gIds)
    gene.data <- rep(1, length(gEns))
    names(gene.data) <- gEns
    for(i in 1:length(kegg$KEGGID)){
        pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[i], species="hsa", out.suffix="pathview", kegg.native=T)
    }
}
    GO_and_KEGG_enrichment( "%(GO_and_KEGG_enrichment_output_path_p_up)s"   , "%(diffgene_genenamelist_path_p_up)s"   )
    GO_and_KEGG_enrichment( "%(GO_and_KEGG_enrichment_output_path_p_down)s" , "%(diffgene_genenamelist_path_p_down)s" )
    GO_and_KEGG_enrichment( "%(GO_and_KEGG_enrichment_output_path_q_up)s"   , "%(diffgene_genenamelist_path_q_up)s"   )
    GO_and_KEGG_enrichment( "%(GO_and_KEGG_enrichment_output_path_q_down)s" , "%(diffgene_genenamelist_path_q_down)s" )
        """ % ( 
                {
                    "GO_and_KEGG_enrichment_output_path_p_up"       :  GO_and_KEGG_enrichment_output_path_p_up,
                    "diffgene_genenamelist_path_p_up"               :  diffgene_genenamelist_path_p_up,
                    "GO_and_KEGG_enrichment_output_path_p_down"     :  GO_and_KEGG_enrichment_output_path_p_down,
                    "diffgene_genenamelist_path_p_down"             :  diffgene_genenamelist_path_p_down,
                    "GO_and_KEGG_enrichment_output_path_q_up"       :  GO_and_KEGG_enrichment_output_path_q_up,
                    "diffgene_genenamelist_path_q_up"               :  diffgene_genenamelist_path_q_up,
                    "GO_and_KEGG_enrichment_output_path_q_down"     :  GO_and_KEGG_enrichment_output_path_q_down,
                    "diffgene_genenamelist_path_q_down"             :  diffgene_genenamelist_path_q_down
                }
            )

        r_script_GO_and_KEGG_enrichment_path_filename = os.path.join(GO_and_KEGG_enrichment_output_path, "r_script_GO_and_KEGG_enrichment")
        with open(r_script_GO_and_KEGG_enrichment_path_filename, "w") as r_script_GO_and_KEGG_enrichment_filehand:
            r_script_GO_and_KEGG_enrichment_filehand.write(r_script_GO_and_KEGG_enrichment)
        os.system("Rscript " + r_script_GO_and_KEGG_enrichment_path_filename)
        os.system("rm " + r_script_GO_and_KEGG_enrichment_path_filename)
    

    else:
        print
        print "=" * 35 + " 7. GO and KEGG  enrichmant analysis (by R) does not run " + "=" * 35 + "\n"
        with open("Run_time", "a") as run_time:
            run_time.write("7. GO_and_KEGG_enrichment_analysis does not run." + "\n\n")



    with open("Run_time", "a") as run_time:
        run_time.write("\n" + "Analyzer_RNA-Seq ends at " + time('%Y-%m-%d %H:%M:%S') + ".\n")



else:
print "Error: There is not enough parameters!"
