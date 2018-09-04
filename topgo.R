
library(topGO)
getwd()
setwd("/Users/Toyin/Desktop/GO_results/GO")

enrichGo <- function(file_name){

  deg_list <- read.delim(paste(file_name,".txt",sep=""),header=FALSE,sep="\t")

  deg_list <- deg_list[2:dim(deg_list)[1],]

  #deg_list <- deg_list[deg_list[,2] == module_num,]
  

  ##colnames(deg_list) <- c("id","baseMean","baseMeanA","baseMeanB","foldChange","log2FoldChange","pval","padj","name","description")

  geneID2GO <- readMappings(file = "MaizeGO2.txt")
  geneNames <- names(geneID2GO)

  Set1GeneList = factor(as.integer(geneNames %in% deg_list[,1]))
  names(Set1GeneList) = geneNames
  GOData_Set1BP = new("topGOdata", ontology="BP", allGenes=Set1GeneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)
  allGO = genesInTerm(GOData_Set1BP)
  Res_Set1BP_ClsFis = runTest(GOData_Set1BP, algorithm="classic", statistic = "fisher")
  Res_Set1BP_ClsFis
  str(score(Res_Set1BP_ClsFis))

  ResTable_Set1BP = GenTable(GOData_Set1BP, Fis = Res_Set1BP_ClsFis, topNodes=Res_Set1BP_ClsFis@geneData[[4]], numChar = 100)
  ResTable_Set1BP <- ResTable_Set1BP[as.numeric(ResTable_Set1BP[,6]) < 0.05,]
  ResTable_Set1BP
  
  write.table(ResTable_Set1BP, file=paste(file_name,"_GOenrichment.txt", sep = ""), quote=FALSE, sep="\t", row.names = FALSE)
}

#enrichGo("deg_ab")
#enrichGo("deg_ab")
#enrichGo("deg_ac")
#enrichGo("deg_ad")
#enrichGo("deg_ae")
#enrichGo("deg_bc")
#enrichGo("deg_bd")
#enrichGo("deg_be")
#enrichGo("deg_cd")
#enrichGo("deg_ce")
#enrichGo("deg_de")


#sg <- sigGenes(GOData_Set1BP)
#str(sg)
#numSigGenes(GOData_Set1BP)

#out <- capture.output(genesInTerm(GOData_Set1BP))
#cat("",out,file="kk.txt",sep="\n")
#write(genesInTerm(GOData_Set1BP)$"GO:0009987",file="kk2.txt")



allRes <- GenTable( deg_ab_GOenrichment , classicFisher = resFisher , topNodes = 20 )
pValue.classic <- score("deg_ab")
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl")
head(listDatasets(ensembl))
listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="aplatyrhynchos_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go_id', values = 'GO:0015979', mart = ensembl)
