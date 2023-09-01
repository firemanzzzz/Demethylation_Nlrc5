library(DESeq2)
library(readr)
library(umap)
library(Rtsne)

#Load raw STAR count data
TD_STAR_COUNT <- as.data.frame(read_csv("~/TDM_STAR_COUNT.csv"))
TRED_STAR_COUNT <- as.data.frame(read_csv("~/TDMa_STAR_COUNT.csv"))
VIVO_STAR_COUNT <- as.data.frame(read_csv("~/VIVO_STAR_COUNT.csv"))

#Load Ensmbl ID to gene name table
ENMUSG_to_GeneName <- as.data.frame(read_delim("~/gene_id_name.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE))

#Check for duplicates, which is not allowed for row names
ENMUSG_to_GeneName <- ENMUSG_to_GeneName[-(which(duplicated(ENMUSG_to_GeneName$`Gene stable ID`))),]

#Rearrange row.names
row.names(TD_STAR_COUNT) <- TD_STAR_COUNT$GENE_ID
row.names(TRED_STAR_COUNT) <- TRED_STAR_COUNT$GENE_ID
row.names(VIVO_STAR_COUNT) <- VIVO_STAR_COUNT$GENE_ID
row.names(ENMUSG_to_GeneName) <- ENMUSG_to_GeneName$`Gene stable ID`
TD_STAR_COUNT <- TD_STAR_COUNT[,-1]
TRED_STAR_COUNT <- TRED_STAR_COUNT[,-1]
VIVO_STAR_COUNT <- VIVO_STAR_COUNT[,-1]
temp <- as.data.frame(ENMUSG_to_GeneName[,-1])
row.names(temp) <- row.names(ENMUSG_to_GeneName)
ENMUSG_to_GeneName <- temp
colnames(ENMUSG_to_GeneName) <- c("Gene_name","NCBI_ID")

#Sort Raw expression table to remove genes with mean count of < 5
TD_STAR_COUNT <- TD_STAR_COUNT[which(rowSums(TD_STAR_COUNT)>=20),]
TRED_STAR_COUNT <- TRED_STAR_COUNT[which(rowSums(TRED_STAR_COUNT)>=20),]
VIVO_STAR_COUNT <- VIVO_STAR_COUNT[which(rowSums(VIVO_STAR_COUNT)>=40),]

#Check whether gene expression table is usable
TD_STAR_COUNT[which(row.names(TD_STAR_COUNT)=="ENSMUSG00000074151"),]
TRED_STAR_COUNT[which(row.names(TRED_STAR_COUNT)=="ENSMUSG00000074151"),]
VIVO_STAR_COUNT[which(row.names(VIVO_STAR_COUNT)=="ENSMUSG00000074151"),]

#Creat a dummy data.frame which for grouping, which is required by DESeq2
TD.group <- as.data.frame(c("E","E","N","N"))
TRED.group <- as.data.frame(c("E","E","N","N"))
VIVO.group <- as.data.frame(c("E","E","E","E","N","N","N","N"))

colnames(TD.group)[1] <- "group"
colnames(TRED.group)[1] <- "group"
colnames(VIVO.group)[1] <- "group"

row.names(TD.group) <- colnames(TD_STAR_COUNT)
row.names(TRED.group) <- colnames(TRED_STAR_COUNT)
row.names(VIVO.group) <- colnames(VIVO_STAR_COUNT)

#Run Deseq2 and get normalized counts
#TD cells
TD.dds <- DESeqDataSetFromMatrix(countData = TD_STAR_COUNT,
                                  colData = TD.group,
                                  design = ~ group)
TD.dds <- DESeq(TD.dds)
resultsNames(TD.dds)
TD.res <- results(TD.dds, name="group_N_vs_E")
TD.res <- as.data.frame(TD.res)
TD.count <- counts(TD.dds, normalized=TRUE)
TD.res[which(row.names(TD.res)=="ENSMUSG00000074151"),]
#TRED cells
TRED.dds <- DESeqDataSetFromMatrix(countData = TRED_STAR_COUNT,
                                   colData = TRED.group,
                                   design = ~ group)
TRED.dds <- DESeq(TRED.dds)
resultsNames(TRED.dds)
TRED.res <- results(TRED.dds, name="group_N_vs_E")
TRED.res <- as.data.frame(TRED.res)
TRED.count <- counts(TRED.dds, normalized=TRUE)
#TRED in vivo tumors
VIVO.dds <- DESeqDataSetFromMatrix(countData = VIVO_STAR_COUNT,
                                   colData = VIVO.group,
                                   design = ~ group)
VIVO.dds <- DESeq(VIVO.dds)
resultsNames(VIVO.dds)
VIVO.res <- results(VIVO.dds, name="group_N_vs_E")
VIVO.res <- as.data.frame(VIVO.res)
VIVO.count <- counts(VIVO.dds, normalized=TRUE)

TD.name <- data.frame(ENMUSG_to_GeneName[which(row.names(ENMUSG_to_GeneName) %in% row.names(TD.res)),], 
                       row.names = row.names(ENMUSG_to_GeneName)[which(row.names(ENMUSG_to_GeneName) %in% row.names(TD.res))])

colnames(TD.name) <- "Gene_name"

#Merge DEG results with table of gene names for convenience
TD.res <- merge(TD.res,ENMUSG_to_GeneName,by="row.names",all=F)
TRED.res <- merge(TRED.res,ENMUSG_to_GeneName,by="row.names",all=F)
VIVO.res <- merge(VIVO.res,ENMUSG_to_GeneName,by="row.names",all=F)

#Generate DESeq2 normalized count table
TD.count <- as.data.frame(log(TD.count + 1, base = 2))
TRED.count <- as.data.frame(log(TRED.count + 1, base = 2))
VIVO.count <- as.data.frame(log(VIVO.count + 1, base = 2))

#Merge normalized count table  with table of gene names for convenience
TD.count <- merge(TD.count,ENMUSG_to_GeneName,by="row.names",all.x=T)
TRED.count <- merge(TRED.count,ENMUSG_to_GeneName,by="row.names",all.x=T)
VIVO.count <- merge(VIVO.count,ENMUSG_to_GeneName,by="row.names",all.x=T)

#Export tables
write.csv(TD.res, file = "~/TD.res.csv", row.names = FALSE)
write.csv(TRED.res, file = "~/TRED.res.csv", row.names = FALSE)
write.csv(VIVO.res, file = "~/VIVO.res.csv", row.names = FALSE)

write.csv(TD.count, file = "~/TD.count.csv", row.names = FALSE)
write.csv(TRED.count, file = "~/TRED.count.csv", row.names = FALSE)
write.csv(VIVO.count, file = "~/VIVO.count.csv", row.names = FALSE)

