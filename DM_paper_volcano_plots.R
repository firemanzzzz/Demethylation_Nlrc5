library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

###TRED-I in vivo tumor differentially expressed gene analysis volcano plot (Figure 4C)###
VIVO.res.vol <- VIVO.res

#Mark genes with customized threshold
VIVO.res.vol$gene_type <- case_when(VIVO.res.vol$log2FoldChange >= 1 & VIVO.res.vol$padj <= 1e-2 ~ "up",
                                    VIVO.res.vol$log2FoldChange <= -1 & VIVO.res.vol$padj <= 1e-2 ~ "down",
                                    TRUE ~ "ns")
#Limit and adjust very small p values
VIVO.res.vol$padj[which(log10(VIVO.res.vol$padj) < -50)] <- 10^(-50)
#Mark the genes of interest
VIVO.res.NLRC5 <- rbind(VIVO.res.vol[which(VIVO.res.vol$Gene_name == "Nlrc5"),],
                        VIVO.res.vol[which(VIVO.res.vol$Gene_name == "H2-D1"),],
                        VIVO.res.vol[which(VIVO.res.vol$Gene_name == "Psmb9"),],
                        VIVO.res.vol[which(VIVO.res.vol$Gene_name == "Tap1"),],
                        VIVO.res.vol[which(VIVO.res.vol$Gene_name == "B2m"),],
                        VIVO.res.vol[which(VIVO.res.vol$Gene_name == "Gzmb"),],
                        VIVO.res.vol[which(VIVO.res.vol$Gene_name == "Stat1"),]
)
VIVO.res.ki67 <- VIVO.res.vol[which(VIVO.res.vol$Gene_name == "Mki67"),]

#Generate the volcano plot (Figure 3c)
VIVO_vol_plot <- VIVO.res.vol %>% ggplot(aes(x = log2FoldChange,
                                             y = -log10(padj),
                                             fill = gene_type,    
                                             size = gene_type,
                                             alpha = gene_type)
) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black",
             size = 3) +
  geom_hline(yintercept = 2,
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = c(0.25,0.25,0.25)) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-6, 6, 2)),       
                     limits = c(-6, 6)) +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=20, family = "sans"),
        axis.title=element_text(size=22, family = "sans"),
        plot.title = element_text(size=24,hjust = 0.5, family = "sans"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=20, family = "sans"),
        legend.text = element_text(size=20, family = "sans")
  ) +
  geom_label_repel(data = VIVO.res.NLRC5, # Add labels last to appear as the top layer  
                   aes(label = Gene_name,fontface = "italic"),
                   force = 2,
                   nudge_y = 5,
                   size=6,
                   fill = "#FFCCCC",
                   colour = "black",
                   alpha=1) +
  geom_point(data = VIVO.res.NLRC5, # New layer containing data subset il_genes       
             size = 4,
             shape = 22,
             fill = "tomato",
             colour = "black",
             alpha=2) +
  labs(
    x = "log2(fold change)",
    y = "-log10(FDR)",
    fill = "Expression \nchange")+
  geom_label_repel(data = VIVO.res.ki67, # Add labels last to appear as the top layer  
                   aes(label = Gene_name,fontface = "italic"),
                   force = 2,
                   nudge_y = 5,
                   size=6,
                   fill = "#9999FF",
                   colour = "black",
                   alpha=1) +
  geom_point(data = VIVO.res.ki67, # New layer containing data subset il_genes       
             size = 4,
             shape = 22,
             fill = "steelblue1",
             colour = "black",
             alpha=2)
#View plot
VIVO_vol_plot

###Pathway enrichment analysis in TRED-I tumors (Figure 4D)###
## feature 1: numeric vector
AMT.List <- AMT.res$log2FoldChange
## feature 2: named vector
names(AMT.List) <- as.character(AMT.res$Row.names)
## Remove NA
AMT.List <- na.omit(AMT.List)
## feature 3: decreasing order
AMT.List <- sort(AMT.List, decreasing = TRUE)
#Gene ontology enrichment analysis
AMT.gse <- gseGO(geneList=AMT.List, 
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 20, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.5, 
                  verbose = TRUE, 
                  OrgDb = org.Mm.eg.db, 
                  pAdjustMethod = "BH")

#Extract GO enrichment analysis result#
GO.res.vol <- AMT.gse@result
#Filter significantly upregulated or downregulated pathways#
GO.res.vol$gene_type <- case_when(GO.res.vol$NES >= 1.5 & GO.res.vol$pvalue <= 0.02 ~ "up",
                                  GO.res.vol$NES <= -1.5 & GO.res.vol$pvalue <= 0.02 ~ "down",
                                  TRUE ~ "ns")
#Differentiate dot sizes
do_size <- GO.res.vol$setSize
#Limit and adjust very small p values
GO.res.vol$pvalue[which(log2(GO.res.vol$pvalue) < -8)] <- 2^(-8)
#Highlight top upregulated pathways 
GO.res.NLRC5 <- GO.res.vol[which(GO.res.vol$NES > 1.636382),]
#Generate volcano plot
GO_vol_plot <- GO.res.vol %>% ggplot(aes(x = NES,
                                         y = -log2(pvalue),
                                         fill = gene_type,    
                                         alpha = gene_type)
) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black",
             size = GO.res.vol$setSize/100) +
  geom_hline(yintercept = -log10(0.02),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = do_size) + # Modify point size
  scale_alpha_manual(values = c(0.25,0.25,0.25)) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-2, 2, 1)),       
                     limits = c(-2, 2),
                     trans = squish_trans(-1, 1, 8)) +
  scale_y_continuous(breaks = c(seq(0, 8.1, 1)),       
                     limits = c(3, 8)) +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=26, family = "sans"),
        axis.title=element_text(size=26, family = "sans"),
        plot.title = element_text(size=24,hjust = 0.5, family = "sans"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=20, family = "sans"),
        legend.text = element_text(size=20, family = "sans")
  ) +
  labs(x = "NES",
       y = "-log2(p-Value)",
       fill = "Expression \nchange")+
  geom_label_repel(data = GO.res.NLRC5, # Add labels last to appear as the top layer  
                   aes(label = Description),
                   force = 2,
                   size=4,
                   fill = NA,
                   colour = "black",
                   nudge_y = -0.25,
                   alpha=1,
                   max.overlaps=8,
                   max.time=1000,
                   family="sans",
                   seed = 2) +
  geom_point(data = GO.res.NLRC5, # New layer containing data subset il_genes       
             size = GO.res.NLRC5$setSize/10,
             shape = 21,
             fill = "chartreuse",
             colour = "black",
             alpha=1)
#View volcano plot
GO_vol_plot

###TRED-I cell line differentially expressed gene analysis volcano plot (Figure S4A)###
AMT.res.vol <- AMT.res
#Filter genes with mean count of >200
AMT.res.vol <- AMT.res.vol[which(AMT.res.vol$baseMean > 128),]
#Mark genes with customized threshold
AMT.res.vol$gene_type <- case_when(AMT.res.vol$log2FoldChange >= 1 & AMT.res.vol$padj <= 1e-2 ~ "up",
                                   AMT.res.vol$log2FoldChange <= -1 & AMT.res.vol$padj <= 1e-2 ~ "down",
                                   TRUE ~ "ns")
#Limit and adjust very small p values
AMT.res.vol$padj[which(log10(AMT.res.vol$padj) < -80)] <- 10^(-80)
#Transform count to size
AMT.res.vol$size <- log(AMT.res.vol$baseMean, base = 2)

#Mark the genes of interest
AMT.res.NLRC5 <- rbind(AMT.res.vol[which(AMT.res.vol$Gene_name == "Nlrc5"),],
                       AMT.res.vol[which(AMT.res.vol$Gene_name == "H2-D1"),],
                       AMT.res.vol[which(AMT.res.vol$Gene_name == "Psmb9"),],
                       AMT.res.vol[which(AMT.res.vol$Gene_name == "Tap1"),],
                       AMT.res.vol[which(AMT.res.vol$Gene_name == "B2m"),],
                       AMT.res.vol[which(AMT.res.vol$Gene_name == "Mx1"),]
)
AMT.res.vol <- AMT.res.vol[-(which(AMT.res.vol$Gene_name %in% AMT.res.NLRC5$Gene_name)),]

#Generate the volcano plot (Figure S4A)
AMT_vol_plot <- AMT.res.vol %>% ggplot(aes(x = log2FoldChange,
                                           y = -log10(padj),
                                           fill = gene_type,    
                                           size = size,
                                           alpha = gene_type)
) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black",
             size = AMT.res.vol$size/4)+
  geom_point(data = AMT.res.NLRC5,
             shape = 21, # Specify shape and colour as fixed local parameters    
             fill = "green",
             size = 3) +
  geom_hline(yintercept = 2,
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = c(0.25,0.25,0.25)) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 5, 1)),       
                     limits = c(-4, 5)) +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=20, family = "sans"),
        axis.title=element_text(size=22, family = "sans"),
        plot.title = element_text(size=24,hjust = 0.5, family = "sans"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=20, family = "sans"),
        legend.text = element_text(size=20, family = "sans")
  ) +
  geom_label_repel(data = AMT.res.vol[c(which(AMT.res.vol$gene_type=="up"),
                                        which(AMT.res.vol$gene_type=="down")),], # Add labels last to appear as the top layer  
                   aes(label = Gene_name,fontface = "italic"),
                   label.size = NA,
                   size=3.5,
                   fill = NA,
                   colour = "black",
                   alpha=1,
                   seed = 123,
                   max.overlaps = 15)+
  geom_label_repel(data = AMT.res.NLRC5, # Add labels last to appear as the top layer  
                   aes(label = Gene_name,fontface = "bold.italic"),
                   size=3.5,
                   fill = "pink",
                   colour = "black",
                   alpha=0.7,
                   seed = 123,
                   max.overlaps = 15)+
  labs(
    x = "log2(fold change)",
    y = "-log10(p-value)",
    fill = "Expression \nchange")
#View plot
AMT_vol_plot
