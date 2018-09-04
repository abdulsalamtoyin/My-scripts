# declare a variable to hold fpkms
fpkm <- c()
# create directory names
dirr <- ("/Users/Toyin/Desktop")
#dirnames <- paste("s", 1:40, sep=""))
# loop through all directories and grab fpkm columns
for( i in (dirr) )
  fname <- paste(dirr[i], "/genes.fpkm_tracking",sep="")
  x <- read.table(file=fname, sep="\t", header=T, as.is=T)
  fpkm <- cbind(fpkm, x[,"FPKM"])

# name the columns
colnames(fpkm) <- dirr
# name the rows, they're all in the same order
rownames(fpkm) <- x[,1]

write.table(fpkm, file="fpkm.txt", sep="\t", col.names=NA)
