library("clusterProfiler")
library("DOSE")
library("enrichplot")
library("org.Mm.eg.db")
library("installr")
library("geneName")

## feature 1: numeric vector
VIVO.List <- VIVO.res$log2FoldChange
## feature 2: named vector
names(VIVO.List) <- as.character(VIVO.res$Row.names)
## Remove NA
VIVO.List <- na.omit(VIVO.List)
## feature 3: decreasing order
VIVO.List <- sort(VIVO.List, decreasing = TRUE)
#Gene ontology enrichment analysis
VIVO.gse <- gseGO(geneList=VIVO.List, 
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 15, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Mm.eg.db, 
                  pAdjustMethod = "BH")
#Export tables
write.csv(VIVO.gse@result, file = "~/VIVO.GO.res.csv", row.names = FALSE)