library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

#TDMa in vivo tumor differentially expressed gene anlaysis volcano plot (Figure 3c)
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


#GO enrichment analysis results (Figure 3d)
#Extract GO enrichment analysis result#
GO.res.vol <- VIVO.gse@result
#Filter significantly upregulated or downregulated pathways#
GO.res.vol$gene_type <- case_when(GO.res.vol$NES >= 2 & GO.res.vol$qvalue <= 1e-3 ~ "up",
                                  GO.res.vol$NES <= -2 & GO.res.vol$qvalue <= 1e-3 ~ "down",
                                  TRUE ~ "ns")
#Differentiate dot sizes
do_size <- GO.res.vol$setSize
#Limit and adjust very small p values
GO.res.vol$padj[which(log10(GO.res.vol$padj) < -50)] <- 10^(-50)
#Highlight top upregulated pathways 
GO.res.NLRC5 <- GO.res.vol[which(GO.res.vol$NES > 2.34),]
#Generate volcano plot
GO_vol_plot <- GO.res.vol %>% ggplot(aes(x = NES,
                                         y = -log10(qvalue),
                                         fill = gene_type,    
                                         size = gene_type,
                                         alpha = gene_type)
) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black",
             size = do_size/50) +
  geom_hline(yintercept = 3,
             linetype = "dashed") + 
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = do_size) + # Modify point size
  scale_alpha_manual(values = c(0.25,0.25,0.25)) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-3, 3, 1)),       
                     limits = c(-3, 3),
                     trans = squish_trans(-2, 2, 8)) +
  scale_y_continuous(breaks = c(seq(0, 9, 1)),       
                     limits = c(0, 9),
                     trans = squish_trans(0, 6, 6)) +
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
       y = "-log10(qValue)",
       fill = "Expression \nchange")+
  geom_label_repel(data = GO.res.NLRC5, # Add labels last to appear as the top layer  
                   aes(label = Description),
                   force = 2,
                   size=6,
                   fill = "#FFCCCC",
                   colour = "black",
                   nudge_y = -1.2,
                   alpha=1,
                   max.overlaps=8,
                   max.time=1000,
                   family="sans",
                   seed = 2) +
  geom_point(data = GO.res.NLRC5, # New layer containing data subset il_genes       
             size = GO.res.NLRC5$setSize/50,
             shape = 21,
             fill = "chartreuse",
             colour = "black",
             alpha=1)
#View volcano plot
GO_vol_plot
