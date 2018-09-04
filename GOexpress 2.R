## R code from vignette source 'GOexpress-UsersGuide.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: GOexpress-UsersGuide.Rnw:24-25
###################################################
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 3: GOexpress-UsersGuide.Rnw:106-108 (eval = FALSE)
###################################################
source("http://bioconductor.org/biocLite.R")
biocLite("GOexpress")
library("GOexpress")


###################################################
### code chunk number 4: GOexpress-UsersGuide.Rnw:119-120
###################################################
maintainer("GOexpress")


###################################################
### code chunk number 5: GOexpress-UsersGuide.Rnw:151-152
###################################################
citation(package="GOexpress")


###################################################
### code chunk number 6: GOexpress-UsersGuide.Rnw:207-209
###################################################
library(GOexpress) # load the GOexpress package
data()
View(genes)
data <- as.matrix(genes)
# import the training dataset
View(data)

###################################################
### code chunk number 7: GOexpress-UsersGuide.Rnw:216-218
###################################################
exprs(data)[1:5,1:5] # Subset of the expression data
testhead(pData(data)) # Subset of the phenotypic information


###################################################
### code chunk number 8: GOexpress-UsersGuide.Rnw:229-230
###################################################
head(rownames(exprs(test))) # Subset of gene identifiers


###################################################
### code chunk number 9: GOexpress-UsersGuide.Rnw:256-258
###################################################
is.factor(AlvMac$Treatment) # assertion test
AlvMac$Treatment # visual inspection


###################################################
### code chunk number 10: GOexpress-UsersGuide.Rnw:268-269 (eval = FALSE)
###################################################
AlvMac$Treatment <- factor(AlvMac$Treatment)


###################################################
### code chunk number 11: GOexpress-UsersGuide.Rnw:291-295
###################################################
set.seed(4543) # set random seed for reproducibility
AlvMac_results <- GO_analyse(
  eSet = AlvMac, f = "Treatment",
  GO_genes=AlvMac_GOgenes, all_GO=AlvMac_allGO, all_genes=AlvMac_allgenes)


###################################################
### code chunk number 12: GOexpress-UsersGuide.Rnw:307-311
###################################################
names(deg_ab_GOenrichment) # Data slot names
head(deg_ab_GOenrichment$GO.ID[, c(1:5, 7)], n=5) # Ranked table of GO terms (subset)
head(AlvMac_results$genes[, c(1:3)], n=5) # Ranked table of genes (subset)
head(deg_ab_GOenrichment$mapping) # Gene to gene ontology mapping table (subset)


###################################################
### code chunk number 13: GOexpress-UsersGuide.Rnw:335-336 (eval = FALSE)
###################################################
AlvMac_results <- GO_analyse(eSet = AlvMac, f = "Treatment")


###################################################
### code chunk number 14: GOexpress-UsersGuide.Rnw:355-356 (eval = FALSE)
###################################################
data(microarray2dataset)


###################################################
### code chunk number 15: GOexpress-UsersGuide.Rnw:371-372 (eval = FALSE)
###################################################
results.pVal = pValue_GO(result=deg_ab_GOenrichment, N=50)
results.pVal = pValue_GO(deg_ab_GOenrichment, N = 1000, ranked.by = deg_ab_GOenrichment$rank.by, rank.by = "P", FUN.GO = deg_ab_GOenrichment$FUN.GO)

###################################################
### code chunk number 16: GOexpress-UsersGuide.Rnw:375-376
###################################################
data(AlvMac_results.pVal)


###################################################
### code chunk number 17: GOexpress-UsersGuide.Rnw:429-444
###################################################
BP.5 <- subset_scores(
  result = AlvMac_results.pVal,
  namespace = "biological_process",
  total = 5, # requires 5 or more associated genes
  p.val=0.05) 
MF.10 <- subset_scores(
  result = AlvMac_results.pVal,
  namespace = "molecular_function",
  total = 10,
  p.val=0.05)
CC.15 <- subset_scores(
  result = AlvMac_results.pVal,
  namespace = "cellular_component",
  total = 15,
  p.val=0.05)


###################################################
### code chunk number 18: GOexpress-UsersGuide.Rnw:478-482 (eval = FALSE)
###################################################
subset(
  AlvMac_results.pVal$GO,
  total_count >= 5 & p.val<0.05 & namespace_1003=='biological_process'
)


###################################################
### code chunk number 19: GOexpress-UsersGuide.Rnw:505-506
###################################################
head(BP.5$GO)


###################################################
### code chunk number 20: GOexpress-UsersGuide.Rnw:521-524
###################################################
heatmap_GO(
  go_id = "GO:0034142", result = BP.5, eSet=AlvMac, cexRow=0.4,
  cexCol=1, cex.main=1, main.Lsplit=30)


###################################################
### code chunk number 21: GOexpress-UsersGuide.Rnw:537-541
###################################################
heatmap_GO(
  go_id = "GO:0034142", result = BP.5, eSet=AlvMac, cexRow=0.4,
  cexCol=1, cex.main=1, main.Lsplit=30,
  labRow=AlvMac$Group)


###################################################
### code chunk number 22: GOexpress-UsersGuide.Rnw:549-553
###################################################
heatmap_GO(
  go_id = "GO:0034142", result = BP.5, eSet=AlvMac, cexRow=0.6,
  cexCol=1, cex.main=1, main.Lsplit=30,
  labRow='Group', subset=list(Time=c('24H','48H')))


###################################################
### code chunk number 23: GOexpress-UsersGuide.Rnw:565-569
###################################################
cluster_GO(
  go_id = "GO:0034142", result = BP.5, eSet=AlvMac,
  cex.main=1, cex=0.6, main.Lsplit=30, 
  subset=list(Time=c("24H", "48H")), f="Group")


###################################################
### code chunk number 24: GOexpress-UsersGuide.Rnw:585-586
###################################################
table_genes(go_id = "GO:0034142", result = BP.5)[,c(1:3)]


###################################################
### code chunk number 25: GOexpress-UsersGuide.Rnw:604-605
###################################################
list_genes(go_id = "GO:0034142", result = BP.5)


###################################################
### code chunk number 26: GOexpress-UsersGuide.Rnw:627-631
###################################################
expression_plot(
  gene_id = "ENSBTAG00000047107", result = BP.5, eSet=AlvMac,
  x_var = "Timepoint", title.size=1.5,
  legend.title.size=10, legend.text.size=10, legend.key.size=15)


###################################################
### code chunk number 27: GOexpress-UsersGuide.Rnw:655-659
###################################################
expression_plot(
  gene_id = "ENSBTAG00000047107", result = BP.5, eSet=AlvMac,
  x_var = "Animal", title.size=1.5, axis.text.angle=90,
  legend.title.size=10, legend.text.size=10, legend.key.size=15)


###################################################
### code chunk number 28: GOexpress-UsersGuide.Rnw:670-674 (eval = FALSE)
###################################################
expression_plot_symbol(
  gene_symbol = "BIKBA", result = BP.5, eSet=AlvMac,
  x_var = "Timepoint", title.size=1.5,
  legend.title.size=10, legend.text.size=10, legend.key.size=15)


###################################################
### code chunk number 29: GOexpress-UsersGuide.Rnw:689-693
###################################################
expression_plot_symbol(
  gene_symbol = "RPL36A", result = AlvMac_results, eSet=AlvMac,
  x_var = "Timepoint", title.size=1.5,
  legend.title.size=10, legend.text.size=10, legend.key.size=15)


###################################################
### code chunk number 30: GOexpress-UsersGuide.Rnw:718-725
###################################################
AlvMac$Animal.Treatment <- paste(AlvMac$Animal, AlvMac$Treatment, sep="_")
expression_profiles(
  gene_id = "ENSBTAG00000047107", result = AlvMac_results,
  eSet=AlvMac, x_var = "Timepoint", line.size=1,
  seriesF="Animal.Treatment", linetypeF="Animal",
  legend.title.size=10, legend.text.size=10,
  legend.key.size=15)


###################################################
### code chunk number 31: GOexpress-UsersGuide.Rnw:740-747 (eval = FALSE)
###################################################
expression_profiles(
  gene_id = "ENSBTAG00000047107", result = AlvMac_results,
  eSet=AlvMac, x_var = "Timepoint",
  lty=rep(1,10), # use line-type 1 for all 10 groups
  seriesF="Animal.Treatment", linetypeF="Animal",
  legend.title.size=10, legend.text.size=10,
  legend.key.size=15, line.size=1)


###################################################
### code chunk number 32: GOexpress-UsersGuide.Rnw:757-763 (eval = FALSE)
###################################################
expression_profiles_symbol(
  gene_symbol="TNIP3", result = AlvMac_results,
  x_var = "Timepoint", linetypeF="Animal", line.size=1,
  eSet=AlvMac, lty=rep(1,10), seriesF="Animal.Treatment", 
  title.size=1.5, legend.title.size=10, legend.text.size=10,
  legend.key.size=15)


###################################################
### code chunk number 33: GOexpress-UsersGuide.Rnw:781-785
###################################################
plot_design(
  go_id = "GO:0034134", result = BP.5, eSet=AlvMac,
  ask = FALSE, factors = c("Animal", "Treatment", "Time", "Group"),
  main.Lsplit=30)


###################################################
### code chunk number 34: GOexpress-UsersGuide.Rnw:834-862 (eval = FALSE)
###################################################
## # Load the interface to BioMart databases
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("biomartr")
library(biomaRt)
library(biomartr)
# See available resources in Ensembl release 75
getMarts()
head(getDatasets(mart = "ENSEMBL_MART_ONTOLOGY") , 20)

listMarts(host='feb2014.archive.ensembl.org')
# Connect to the Ensembl Genes annotation release 75 for Bos taurus
ensembl75 = useMart(
  host='feb2014.archive.ensembl.org',
  biomart='pride', dataset='btaurus_gene_ensembl')
## Download all the Ensembl gene annotations (no filtering)
allgenes.Ensembl = getBM(
  attributes=c('ensembl_gene_id', 'external_gene_id', 'description'),
  mart=ensembl75)
# Rename the gene identifier column to 'gene_id'
# This allows GOexpress to treat microarray and RNA-seq data identically
colnames(allgenes.Ensembl)[1] = 'gene_id'
## Download all the gene ontology annotations (no filtering)
allGO.Ensembl = getBM(
  attributes=c('go_id', 'name_1006', 'namespace_1003'),
  mart=ensembl75)
## Download all the mapping between gene and gene ontology identifiers
GOgenes.Ensembl = getBM(
  attributes=c('ensembl_gene_id', 'go_id'),
  mart=ensembl75)
# Rename the gene identifier column to 'gene_id'
colnames(GOgenes.Ensembl)[1] = 'gene_id'
# Cleanup: remove some blank fields often found in both columns
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]


###################################################
### code chunk number 35: GOexpress-UsersGuide.Rnw:879-889 (eval = FALSE)
###################################################
# save each custom annotation to a R data file
save(GOgenes.Ensembl, file='GOgenes.Ensembl75.rda')
save(allGO.Ensembl, file='allGO.Ensembl75.rda')
save(allgenes.Ensembl, file='allgenes.Ensembl75.rda')
# Run an analysis using those local annotations
GO_analyse(
  eSet=AlvMac, f='Treatment',
  GO_genes=GOgenes.Ensembl,
  all_GO=allGO.Ensembl,
  all_genes=allgenes.Ensembl)


###################################################
### code chunk number 36: GOexpress-UsersGuide.Rnw:897-900 (eval = FALSE)
###################################################
data(AlvMac_GOgenes)
data(AlvMac_allGO)
data(AlvMac_allgenes)


###################################################
### code chunk number 37: GOexpress-UsersGuide.Rnw:932-945 (eval = FALSE)
###################################################
AlvMac_results <- GO_analyse(
  eSet = AlvMac, f = "Treatment",
  subset=list(
    Time=c("6H","24H", "48H"),
    Treatment=c("CN","MB"))
)

expression_plot(
  gene_id = "ENSBTAG00000047107", result = BP.5, eSet=AlvMac,
  x_var = "Timepoint", title.size=1.5,
  legend.title.size=10, legend.text.size=10, legend.key.size=15,
  subset=list(Treatment=c("TB","MB"))
)


###################################################
### code chunk number 38: GOexpress-UsersGuide.Rnw:966-967
###################################################
overlap_GO(go_ids = head(BP.5$GO$go_id, n=5), result = BP.5, filename=NULL)


###################################################
### code chunk number 39: GOexpress-UsersGuide.Rnw:985-986
###################################################
hist_scores(result = BP.5, labels = TRUE)


###################################################
### code chunk number 40: GOexpress-UsersGuide.Rnw:993-994
###################################################
quantiles_scores(result = BP.5)


###################################################
### code chunk number 41: GOexpress-UsersGuide.Rnw:1013-1014 (eval = FALSE)
###################################################
BP.5.byScore <- rerank(result = BP.5, rank.by = "score")


###################################################
### code chunk number 42: GOexpress-UsersGuide.Rnw:1023-1024 (eval = FALSE)
###################################################
BP.5.byPval <- rerank(result = BP.5, rank.by = "p.val")


###################################################
### code chunk number 43: GOexpress-UsersGuide.Rnw:1033-1035 (eval = FALSE)
###################################################
BP.5.pVal_rank <- rerank(result = BP.5, rank.by = "rank")
BP.5.pVal_rank <- rerank(result = BP.5.pVal_rank, rank.by = "p.val")


###################################################
### code chunk number 44: GOexpress-UsersGuide.Rnw:1056-1057
###################################################
pData(AlvMac) <- droplevels(pData(AlvMac))


###################################################
### code chunk number 45: GOexpress-UsersGuide.Rnw:1073-1077 (eval = FALSE)
###################################################
subEset(
  eSet=AlvMac, subset=list(
    Time=c("2H","6H","24H"),
    Treatment=c("CN","MB")))


###################################################
### code chunk number 46: GOexpress-UsersGuide.Rnw:1233-1234
###################################################
sessionInfo()

