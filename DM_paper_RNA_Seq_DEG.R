library(DESeq2)
library(readr)
library(umap)
library(Rtsne)

#Load raw STAR count data
TDM_STAR_COUNT <- as.data.frame(read_csv("~/TDM_STAR_COUNT.csv"))
TDMa_STAR_COUNT <- as.data.frame(read_csv("~/TDMa_STAR_COUNT.csv"))
VIVO_STAR_COUNT <- as.data.frame(read_csv("~/VIVO_STAR_COUNT.csv"))

#Load Ensmbl ID to gene name table
ENMUSG_to_GeneName <- as.data.frame(read_delim("~/gene_id_name.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE))

#Check for duplicates, which is not allowed for row names
ENMUSG_to_GeneName <- ENMUSG_to_GeneName[-(which(duplicated(ENMUSG_to_GeneName$`Gene stable ID`))),]

#Rearrange row.names
row.names(TDM_STAR_COUNT) <- TDM_STAR_COUNT$GENE_ID
row.names(TDMa_STAR_COUNT) <- TDMa_STAR_COUNT$GENE_ID
row.names(VIVO_STAR_COUNT) <- VIVO_STAR_COUNT$GENE_ID
row.names(ENMUSG_to_GeneName) <- ENMUSG_to_GeneName$`Gene stable ID`
TDM_STAR_COUNT <- TDM_STAR_COUNT[,-1]
TDMa_STAR_COUNT <- TDMa_STAR_COUNT[,-1]
VIVO_STAR_COUNT <- VIVO_STAR_COUNT[,-1]
temp <- as.data.frame(ENMUSG_to_GeneName[,-1])
row.names(temp) <- row.names(ENMUSG_to_GeneName)
ENMUSG_to_GeneName <- temp
colnames(ENMUSG_to_GeneName) <- c("Gene_name","NCBI_ID")

#Sort Raw expression table to remove genes with mean count of < 5
TDM_STAR_COUNT <- TDM_STAR_COUNT[which(rowSums(TDM_STAR_COUNT)>=20),]
TDMa_STAR_COUNT <- TDMa_STAR_COUNT[which(rowSums(TDMa_STAR_COUNT)>=20),]
VIVO_STAR_COUNT <- VIVO_STAR_COUNT[which(rowSums(VIVO_STAR_COUNT)>=40),]

#Check whether gene expression table is usable
TDM_STAR_COUNT[which(row.names(TDM_STAR_COUNT)=="ENSMUSG00000074151"),]
TDMa_STAR_COUNT[which(row.names(TDMa_STAR_COUNT)=="ENSMUSG00000074151"),]
VIVO_STAR_COUNT[which(row.names(VIVO_STAR_COUNT)=="ENSMUSG00000074151"),]

#Creat a dummy data.frame which for grouping, which is required by DESeq2
TDM.group <- as.data.frame(c("E","E","N","N"))
TDMa.group <- as.data.frame(c("E","E","N","N"))
VIVO.group <- as.data.frame(c("E","E","E","E","N","N","N","N"))

colnames(TDM.group)[1] <- "group"
colnames(TDMa.group)[1] <- "group"
colnames(VIVO.group)[1] <- "group"

row.names(TDM.group) <- colnames(TDM_STAR_COUNT)
row.names(TDMa.group) <- colnames(TDMa_STAR_COUNT)
row.names(VIVO.group) <- colnames(VIVO_STAR_COUNT)

#Run Deseq2 and get normalized counts
#TDM cells
TDM.dds <- DESeqDataSetFromMatrix(countData = TDM_STAR_COUNT,
                                  colData = TDM.group,
                                  design = ~ group)
TDM.dds <- DESeq(TDM.dds)
resultsNames(TDM.dds)
TDM.res <- results(TDM.dds, name="group_N_vs_E")
TDM.res <- as.data.frame(TDM.res)
TDM.count <- counts(TDM.dds, normalized=TRUE)
TDM.res[which(row.names(TDM.res)=="ENSMUSG00000074151"),]
#TDMa cells
TDMa.dds <- DESeqDataSetFromMatrix(countData = TDMa_STAR_COUNT,
                                   colData = TDMa.group,
                                   design = ~ group)
TDMa.dds <- DESeq(TDMa.dds)
resultsNames(TDMa.dds)
TDMa.res <- results(TDMa.dds, name="group_N_vs_E")
TDMa.res <- as.data.frame(TDMa.res)
TDMa.count <- counts(TDMa.dds, normalized=TRUE)
#TDMa in vivo tumors
VIVO.dds <- DESeqDataSetFromMatrix(countData = VIVO_STAR_COUNT,
                                   colData = VIVO.group,
                                   design = ~ group)
VIVO.dds <- DESeq(VIVO.dds)
resultsNames(VIVO.dds)
VIVO.res <- results(VIVO.dds, name="group_N_vs_E")
VIVO.res <- as.data.frame(VIVO.res)
VIVO.count <- counts(VIVO.dds, normalized=TRUE)

TDM.name <- data.frame(ENMUSG_to_GeneName[which(row.names(ENMUSG_to_GeneName) %in% row.names(TDM.res)),], 
                       row.names = row.names(ENMUSG_to_GeneName)[which(row.names(ENMUSG_to_GeneName) %in% row.names(TDM.res))])

colnames(TDM.name) <- "Gene_name"

#Merge DEG results with table of gene names for convenience
TDM.res <- merge(TDM.res,ENMUSG_to_GeneName,by="row.names",all=F)
TDMa.res <- merge(TDMa.res,ENMUSG_to_GeneName,by="row.names",all=F)
VIVO.res <- merge(VIVO.res,ENMUSG_to_GeneName,by="row.names",all=F)

#Generate DESeq2 normalized count table
TDM.count <- as.data.frame(log(TDM.count + 1, base = 2))
TDMa.count <- as.data.frame(log(TDMa.count + 1, base = 2))
VIVO.count <- as.data.frame(log(VIVO.count + 1, base = 2))

#Merge normalized count table  with table of gene names for convenience
TDM.count <- merge(TDM.count,ENMUSG_to_GeneName,by="row.names",all.x=T)
TDMa.count <- merge(TDMa.count,ENMUSG_to_GeneName,by="row.names",all.x=T)
VIVO.count <- merge(VIVO.count,ENMUSG_to_GeneName,by="row.names",all.x=T)

#Export tables
write.csv(TDM.res, file = "~/TDM.res.csv", row.names = FALSE)
write.csv(TDMa.res, file = "~/TDMa.res.csv", row.names = FALSE)
write.csv(VIVO.res, file = "~/VIVO.res.csv", row.names = FALSE)

write.csv(TDM.count, file = "~/TDM.count.csv", row.names = FALSE)
write.csv(TDMa.count, file = "~/TDMa.count.csv", row.names = FALSE)
write.csv(VIVO.count, file = "~/VIVO.count.csv", row.names = FALSE)

