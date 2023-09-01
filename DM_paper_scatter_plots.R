library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

#Calculate the mean log2 normalized count of replicate samples
TD.count$Scramble_Mean <- (TD.count$TD_sgScramble_1 + TD.count$TD_sgScramble_2)/2
TD.count$Nlrc5_Mean <- (TD.count$TD_sgNlrc5_1 + TD.count$TD_sgNlrc5_2)/2
TRED.count$Scramble_Mean <- (TRED.count$TRED_sgScramble_1 + TRED.count$TRED_sgScramble_2)/2
TRED.count$Nlrc5_Mean <- (TRED.count$TRED_sgNlrc5_1 + TRED.count$TRED_sgNlrc5_2)/2

#Mark the target genes
TRED.names.up <- c("Nlrc5","H2-D1","Psmb9","Tap1","B2m")
TRED.names.down <- c("Stat1","Irf1","Mx1")

#generate the plot for TRED cells (figure 2 d)
TRED_outliner_plot <- ggplot(TRED.count, aes(Scramble_Mean, Nlrc5_Mean)) +
  geom_point(shape = 21, # Specify shape and color as fixed local parameters    
             colour = "grey",
             
             size = 1.5) +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=24, family = "sans"),
        axis.title=element_text(size=26, family = "sans"),
        plot.title = element_text(size=24,hjust = 0.5, family = "sans"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=20, family = "sans"),
        legend.text = element_text(size=20, family = "sans")) +
  geom_point(data = TRED.count[which(TRED.count$Gene_name%in%TRED.names.up),], # New layer      
             size = 5,
             shape = 22,
             fill = "tomato",
             colour = "black",
             alpha=1)+
  geom_label_repel(data = TRED.count[which(TRED.count$Gene_name%in%TRED.names.up),], # Add labels last to appear as the top layer  
                   aes(label = TRED.count[which(TRED.count$Gene_name%in%TRED.names.up),]$Gene_name,fontface = "italic"),
                   force = 2,
                   nudge_y = 3,
                   size=8,
                   fill = "#FFCCCC",
                   colour = "black",
                   alpha=0.6) +
  geom_label_repel(data = TRED.count[which(TRED.count$Gene_name%in%TRED.names.down),], # Add labels last to appear as the top layer  
                   aes(label = TRED.count[which(TRED.count$Gene_name%in%TRED.names.down),]$Gene_name,fontface = "italic"),
                   force = 2,
                   nudge_y = -4,
                   size=8,
                   fill = "#9999FF",
                   colour = "black",
                   alpha=0.6) +
  geom_point(data = TRED.count[which(TRED.count$Gene_name%in%TRED.names.down),], # New layer containing data subset il_genes       
             size = 5,
             shape = 22,
             fill = "steelblue1",
             colour = "black",
             alpha=1) +
  labs(x = "log2(Normalzied count +1) sgScramble",
       y = "log2(Normalzied count +1) sgNlrc5")

TRED_outliner_plot #View plot

#Assign gnee names
TD.names.up <- c("Nlrc5")
#generate the plot for TD cells (figure 1 f)
TD_outliner_plot <- ggplot(TD.count, aes(Scramble_Mean, Nlrc5_Mean)) +
  geom_point(shape = 21, # Specify shape and color as fixed local parameters    
             colour = "grey",
             size = 1.5) +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=24, family = "sans"),
        axis.title=element_text(size=22, family = "sans"),
        plot.title = element_text(size=24,hjust = 0.5, family = "sans"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=20, family = "sans"),
        legend.text = element_text(size=20, family = "sans")) +
  geom_point(data = TD.count[which(TD.count$Gene_name%in%TD.names.up),], # New layer containing data subset il_genes       
             size = 5,
             shape = 21,
             fill = "tomato",
             colour = "black",
             alpha=1)+
  geom_label_repel(data = TD.count[which(TD.count$Gene_name%in%TD.names.up),], # Add labels last to appear as the top layer  
                   aes(label = TD.names.up,fontface = "italic"),
                   force = 2,
                   nudge_y = 4,
                   size=8,
                   fill = "#FFCCCC",
                   colour = "black",
                   alpha=0.6)+
  labs(x = "log2(Normalzied count +1) sgScramble",
       y = "log2(Normalzied count +1) sgNlrc5")

TD_outliner_plot #View plot
