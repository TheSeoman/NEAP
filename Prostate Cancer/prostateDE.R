library("edgeR")

outPath <- "/home/seoman/Documents/NEAP/Prostate Cancer/R_out/"
countsPath <- "/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/"

runEdgeR <- function(dataSetName, cancerCounts, normalCounts, groupSize){
  counts <- merge(cancerCounts, normalCounts, by = 0)
  rownames(counts) <- counts$Row.names
  counts[is.na(counts)] <- 0
  counts <- counts[,-1]
  counts <- counts[rowSums(cpm(counts)>1) >= min(groupSize),]
  group <- c(rep(c("cancer"), groupSize[1]), rep(c("normal"), groupSize[2]))
  
  dge <- DGEList(counts = counts, group = group, genes = data.frame(GeneID = rownames(counts)))
  
  dge <- calcNormFactors(dge)
  dge <- estimateCommonDisp(dge)
  dge <- estimateTagwiseDisp(dge)
  et <- exactTest(dge)
  et$table$adj.PValue <- p.adjust(et$table$PValue, method = "BH")
  et$table <- et$table[order(et$table$adj.PValue), ]
  
  write.table(et$table, file = paste(outPath, dataSetName, "_DE.tsv", sep = ""), sep = "\t", quote = FALSE)
  
#  upGenes <- rownames(et$table[et$table$adj.PValue < 0.05 & et$table$logFC > 0,])
#  downGenes <- rownames(et$table[et$table$adj.PValue < 0.05 & et$table$logFC < 0,])

  upGenes <- rownames(et$table[et$table$logFC > 1,])
  downGenes <- rownames(et$table[et$table$logFC < -1,])
  
    
  write.table(upGenes, file = paste(outPath, dataSetName, "_fcup.tsv", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(downGenes, file = paste(outPath, dataSetName, "_fcdown.tsv", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  et$table
} 

GDS5072cancer <- read.table(paste(countsPath, "GDS5072_cancer", sep = ""), header = TRUE, sep = "\t", row.names = 1)
GDS5072normal <- read.table(paste(countsPath, "GDS5072_normal", sep = ""), header = TRUE, sep = "\t", row.names = 1)

GDS5072DE <- runEdgeR("GDS5072", GDS5072cancer, GDS5072normal, c(10, 1))


GDS4824cancer <- read.table(paste(countsPath, "GDS4824_cancer", sep = ""), header = TRUE, sep = "\t", row.names = 1)
GDS4824normal <- read.table(paste(countsPath, "GDS4824_normal", sep = ""), header = TRUE, sep = "\t", row.names = 1)

GDS4824DE <- runEdgeR("GDS4824", GDS4824cancer, GDS4824normal, c(13, 8))


GDS4114cancer <- apply(read.table(paste(countsPath, "GDS4114_cancer", sep = ""), header = TRUE, sep = "\t", row.names = 1), 2, exp)
GDS4114normal <- apply(read.table(paste(countsPath, "GDS4114_normal", sep = ""), header = TRUE, sep = "\t", row.names = 1), 2, exp)

GDS4114DE <- runEdgeR("GDS4114", GDS4114cancer, GDS4114normal, c(6, 6))


GDS4395cancer <- read.table(paste(countsPath, "GDS4395_cancer", sep = ""), header = TRUE, sep = "\t", row.names = 1)
GDS4395normal <- read.table(paste(countsPath, "GDS4395_normal", sep = ""), header = TRUE, sep = "\t", row.names = 1)

GDS4395DE <- runEdgeR("GDS4395", GDS4395cancer, GDS4395normal, c(10, 10))

GDS2545cancer <- read.table(paste(countsPath, "GDS2545_cancer", sep = ""), header = TRUE, sep = "\t", row.names = 1)
GDS2545normal <- read.table(paste(countsPath, "GDS2545_normal", sep = ""), header = TRUE, sep = "\t", row.names = 1)

GDS2545DE <- runEdgeR("GDS2545", GDS2545cancer, GDS2545normal, c(90, 18))

