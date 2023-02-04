setwd("/Volumes/xiaohong/imac/SAR/Final/Codes/")

rm(list=ls())
graphics.off()


# load packages -----------------------------------------------------------


library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ape)
library(ggsci)
library(vegan)
library(Rmisc)
library(stringr)
library(pheatmap)


# load data ---------------------------------------------------------------

metadata <- read.table("../Files/metadata.txt", header = T, sep = "\t")
metadata$group <- gsub("C", "Con", metadata$group)
metadata$group <- gsub("S", "SAR", metadata$group)
metadata$group_id <- paste0(metadata$group_id, ".2")

rownames(metadata) <- metadata$group_id


metaphlan <- read.table("../Files/merged_metaphlan_modifed.txt", header = T, sep = "\t")
rownames(metaphlan) <- metaphlan$SampleID
metaphlan$SampleID <- NULL


# Figure 1 ----------------------------------------------------------------

# permanova
data_microbiome <- data.frame(t(metaphlan))

adonis2(data_microbiome ~ race_iden + group, data = metadata, permutations = 100)

# a diversity
# shannon index
shan_index <- diversity(t(metaphlan), "shannon")

# generate a table
a_diversity <- data.frame(colnames(metaphlan), shan_index)

colnames(a_diversity)[1] <- "group_id"
meta2com <- select(metadata, group_id, group, race_iden)

adiversity <- merge(a_diversity, meta2com, by = "group_id")
colnames(adiversity)[3:4] <-c("Group", "Race")
adiversity$Race <- gsub("han", "Han", adiversity$Race)
adiversity$Race <- gsub("zang", "Zang", adiversity$Race)
adiversity$Race <- gsub("qiang", "Qiang", adiversity$Race)

# make boxplot
p_diversity_group <- ggplot(adiversity, aes(x = as.factor(Group), y = shan_index, color = Group)) + 
  geom_boxplot() + geom_jitter(size = 0.5) +
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title = element_text(hjust = 0.5)) + theme_test(base_size = 8) + 
  scale_color_manual(values = c("#89c3eb", "#c97586")) + theme(legend.position = "none") 

pdf("../Figures/a_diversity_group.pdf", width = 4.3/2.54, height = 4.6/2.54)
p_diversity_group
graphics.off()

p_diversity_race <- ggplot(adiversity, aes(x = as.factor(Race), y = shan_index, color = Race)) + 
  geom_boxplot() + geom_jitter(size = 0.5) +
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title = element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_color_manual(values = c("#e44a35",  "#d96f0e", "#2ca25f")) + theme(legend.position = "none")

pdf("../Figures/a_diversity_race.pdf", width = 4.3/2.54, height = 4.6/2.54)
p_diversity_race
graphics.off()

# PCoA
# distance
dist <- vegdist(t(metaphlan), method = "bray")

# PCoA
pcoa <- pcoa(dist)

# choose the first and second axis
pcoa_plot <- pcoa$vectors[,1:2] %>% as.data.frame()

colnames(pcoa_plot) <- c("Axis.1", "Axis.2")
sample <- row.names(pcoa_plot)
pcoa_plot$sample <- sample
colnames(pcoa_plot)[3] <- "group_id"

# make pcoa plot
pcoa2plot <- merge(pcoa_plot, metadata, by = "group_id")

p_pcoa <- ggplot(pcoa2plot, aes(x = Axis.1, y = Axis.2, color = race_iden )) + geom_point(size = 0.5) +
  stat_ellipse(level=0.5) + theme_test(base_size = 8) +
  scale_color_manual(values = c("#e44a35",  "#d96f0e", "#2ca25f" ))  + theme(legend.position = "none")

pdf("../Figures/poca_race.pdf", width = 4 /2.54, height = 4.3/2.54)
p_pcoa
graphics.off()

# make boxplot of axis 1 and axis, respectively
poca_box1 <- ggplot(pcoa2plot, aes(x = race_iden, y = Axis.1, colour = race_iden)) + geom_boxplot() +
  scale_color_manual(values = c("#e44a35",  "#d96f0e", "#2ca25f" )) + theme_test(base_size = 8) + theme(legend.position = "none") +
  ylab("") + xlab("")

pdf("../Figures/box_pcoa_race_1.pdf", width = 2/2.54, height = 4.3/2.54) 
poca_box1 
graphics.off()


poca_box2 <- ggplot(pcoa2plot, aes(x = race_iden, y = Axis.2, colour = race_iden)) + geom_boxplot() +
  scale_color_manual(values = c("#e44a35",  "#d96f0e", "#2ca25f" )) + theme_test(base_size = 8) + theme(legend.position = "none") +
  ylab("") + xlab("")

pdf("../Figures/box_pcoa_race_2.pdf", width = 2/2.54, height = 4.3/2.54)
poca_box2
graphics.off()

# make taxonomy barplot
bar_plot <- function(micro_data, level, top_num, meta_data, sample, group){
  metadata <- meta_data[,c(which(colnames(meta_data) == sample), which(colnames(meta_data) == group))]
  colnames(metadata) <- c("Sample", "Group")
  rownames(metadata) <- metadata$Sample
  
  if(level == "phylum"){rownames(micro_data) <- str_split_fixed(rownames(micro_data), "\\|", n = 6)[,2]}
  if(level == "genus"){rownames(micro_data) <- str_split_fixed(rownames(micro_data), "\\|", n = 6)[,6]}
  micro_data = micro_data[(order(-rowSums(micro_data))), ]
  
  # combine low abundance taxonomy into Low abundance
  other = colSums(micro_data[top_num : dim(micro_data)[1], ])
  micro_data = micro_data[1:(top_num - 1), ]
  micro_data = rbind(micro_data,other)
  rownames(micro_data)[top_num] = c("Others")
  
  sample_order <- data.frame(t(micro_data))
  micro_data$Taxonomy <- rownames(micro_data)
  
  data2plot <- data.frame(melt(micro_data, id.vars = c("Taxonomy")))
  data2plot <- merge(data2plot, metadata, by.x = "variable", by.y = "Sample")
  
  
  sample_order$variable <- as.factor(rownames(sample_order))
  data4plot <- left_join(sample_order, data2plot, by = "variable")
  data4plot <- data4plot[order(data4plot[1], decreasing = T), ]
  data4plot$variable <- factor(data4plot$variable, levels = unique(data4plot$variable))
  
  p <- ggplot(data4plot, aes(x = factor(variable), y = value, fill = factor(Taxonomy, levels = rownames(micro_data)[top_num : 1]))) + 
    geom_bar(stat = "identity",position="fill", width=1) + theme_test(base_size = 8) +
    scale_y_continuous(labels = scales::percent) + 
    facet_grid( ~ factor(Group), scales = "free_x", switch = "x" ) +
    theme(strip.background = element_blank(), legend.position = "none") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    xlab(group) + ylab("Percentage")+
    labs(fill="Taxonomy") + scale_fill_manual(values = c("#BFC096", "#4dbbd5ff", "#3c54bbff",
                                                         "#f39b7fff", "#8491b4ff", "#91d1c2ff",
                                                         "#00a087ff", "#7e6148ff", "#b09c85ff",
                                                         "#e64b35ff" ))
  
  
  # average in group
  data2group <- aggregate(data4plot[1:top_num], by = data4plot[top_num+4], FUN = mean)
  data4group = as.data.frame(melt(data2group, id.vars="Group"))
  data4group$variable <- gsub("Low.abundance", "Low abundance", data4group$variable)
  
  p_phylum_group = ggplot(data4group, aes(x= Group, y = value, fill = factor(variable, levels = rownames(micro_data)[top_num : 1]))) + 
    geom_bar(stat = "identity",position="fill", width=0.7)+ 
    scale_y_continuous(labels = scales::percent) + 
    xlab(group)+ylab("Percentage") + theme_test(base_size = 8) + theme(legend.position = "right")  +
    labs(fill="Taxonomy") + scale_fill_npg()
  
  return(list(p_phylum_group, p))
} 


colnames(metadata)[2] <- "Ethnic groups"
metadata$`Ethnic groups` <- gsub("han", "Han", metadata$`Ethnic groups`)
metadata$`Ethnic groups` <- gsub("zang", "Zang", metadata$`Ethnic groups`)
metadata$`Ethnic groups` <- gsub("qiang", "Qiang", metadata$`Ethnic groups`)



p <- bar_plot(micro_data = metaphlan, level = "genus", top_num = 10, meta_data = metadata, sample = 'group_id', group = 'Ethnic groups')

pdf("../Figures/barplot_race_group.pdf", width = 8.6/2.54, height = 6/2.54)
p[1]
graphics.off()

pdf("../Figures/barplot_race_indi.pdf", width = 8.6/2.54, height = 6/2.54)
p[2]
graphics.off()





# Figure 3 ----------------------------------------------------------------


# taxonomy barplot
colnames(metadata)[4] <- "Groups"
p_groups <- bar_plot(micro_data = tax_sample, level = "genus", top_num = 10, meta_data = metadata, sample = 'group_id', group = 'Groups')


pdf("../Figures/Fig3/barplot_group_group.pdf", width = 7.1/2.54, height = 6/2.54)
p_groups[[1]]
graphics.off()

tax_sample = read.table("../Files/phylum.txt", header=T, row.names= 1, sep="\t", comment.char="") 

p_groups_phylum <- bar_plot(micro_data = tax_sample, level = "phylum", top_num = 6, meta_data = metadata, sample = 'group_id', group = 'Groups')

pdf("../Figures/Fig3/barplot_group_group_phylum.pdf", width = 6.0/2.54, height = 6/2.54)
p_groups_phylum[[1]]
graphics.off()

# Prevotella copri

p_data <- read.table("../Files/pcopri.txt",header = T, sep = "\t")

p_data_plot <- melt(p_data, id = c("sample_id", "Group"))
colnames(p_data_plot) <- c("sample_id", "Group", "P_copri", "Relative abundance")

p_plot <- ggplot(p_data_plot, aes(x = Group, y = `Relative abundance`, fill = Group)) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  geom_boxplot(outlier.shape =  NA, notch = T, outlier.size = NA, width = 0.3,
               notchwidth = 0.2 )+ ylim(0,80) + xlab("P copri") +
  theme_test(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "") 
pdf("../Figures/Fig3/barplot_pcopri.pdf", width = 4.3/2.54, height = 6/2.54)
p_plot
graphics.off()



# test
# p copri
data_p_con <- subset(p_data_plot, Group == "Con")
data_p_con_rela <- data_p_con$`Relative abundance`

data_p_sar <- subset(p_data_plot, Group == "SAR")
data_p_sar_rela <- data_p_sar$`Relative abundance`

wilcox.test(data_p_con_rela, data_p_sar_rela, alternative = "two.sided")

mean(data_p_sar_rela)
mean(data_p_con_rela)

sd(data_p_sar_rela)
sd(data_p_con_rela)

# e rectale
data_e_con <- subset(e_data_plot, Group == "Con")
data_e_con_rela <- data_e_con$`Relative abundance`

data_e_sar <- subset(e_data_plot, Group == "SAR")
data_e_sar_rela <- data_e_sar$`Relative abundance`

wilcox.test(data_e_con_rela, data_e_sar_rela, alternative = "two.sided")



# pathway
data_pathway <- read.table("../data/2022_3_16_pathway.txt", header = T, sep = "\t")
data_pathway$OR <- data_pathway$GeneRatio / data_pathway$BgRatio
data_pathway$lOR <- log2(data_pathway$OR)

plot_or <-ggplot(data_pathway, aes(OR, factor(Pathway, levels = c("Metabolism of other amino acids",
                                                                  "Amino acid metabolism",
                                                                  "Secretion system",
                                                                  "Fatty acid degradation",
                                                                  "Bacterial secretion system",
                                                                  "Valine,leucine,and isoleucine biosynthesis")))) + 
  geom_point(aes(color=-1*log10(qvalue))) + 
  scale_colour_gradient(low="blue", high="red") + 
  labs(color = expression(-log10(q)), x ="Odds ratio", y = "") +
  theme_bw(base_size = 8)
plot_or

pdf("../Figures/Fig3/pathway_scatter_or.pdf", width = 8 /2.54, height = 6.5/2.54)
plot_or
graphics.off()




# Figure 4 ----------------------------------------------------------------


# metabolite  QC
# positive mode
pos <- read.table("../Files/pos-PCA_t.txt", sep = "\t", header = T)
rownames(pos) <- pos$ID
pos$ID <- NULL
distance_pos <- dist(pos)
pcoa_pos <- pcoa(distance_pos)
pcoa_plot_pos <- pcoa_pos$vectors[,1:2] %>% as.data.frame()
pcoa_plot_pos$type <- rownames(pos)
pcoa_plot_pos$Class <- substr(pcoa_plot_pos$type, 1, 1)

pcoa_plot_pos$Class <- gsub("C", "S", pcoa_plot_pos$Class)
pcoa_plot_pos$Class <- gsub("S", "Sample", pcoa_plot_pos$Class)
pcoa_plot_pos$Class <- gsub("Q", "QC", pcoa_plot_pos$Class)

p_pcoa_1 <- ggplot(pcoa_plot_pos, aes(x = Axis.1, y = Axis.2, colour = Class)) + geom_point(size = 0.5) +
  ggtitle("") + theme_test(base_size = 8) + 
  theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values = c("#88ABDA","#C35C6A", "#9EBC19")) +
  stat_ellipse() + theme(legend.position = c(0.7,0.7))

pdf("../Figures/Fig4/pcoa_pos.pdf", width = 5.8 /2.54, height = 4.8/2.54)
p_pcoa_1 
graphics.off()


# negetive mode
neg <- read.table("../Files/neg-PCA_t.txt", sep = "\t", header = T)
rownames(neg) <- neg$ID
neg$ID <- NULL

distance_neg <- dist(neg)
pcoa_neg <- pcoa(distance_neg)

pcoa_plot_neg <- pcoa_neg$vectors[,1:2] %>% as.data.frame()
pcoa_plot_neg$type <- rownames(neg)
pcoa_plot_neg$Class <- substr(pcoa_plot_neg$type, 1, 1)

pcoa_plot_neg$Class <- gsub("C", "S", pcoa_plot_neg$Class)

pcoa_plot_neg$Class <- gsub("S", "Sample", pcoa_plot_neg$Class)
pcoa_plot_neg$Class <- gsub("Q", "QC", pcoa_plot_neg$Class)

p_pcoa_2 <- ggplot(pcoa_plot_neg, aes(x = Axis.1, y = Axis.2, colour = Class)) + geom_point(size = 0.5) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme_test(base_size = 8) + scale_colour_manual(values = c("#88ABDA","#C35C6A", "#9EBC19")) +
  stat_ellipse() + theme(legend.position = "none")

pdf("../Figures/Fig4/pcoa_neg.pdf", width = 5.8/2.54, height = 4.8/2.54)
p_pcoa_2
graphics.off()

# make volcano plot
metabolite_data <- read.table("../Files/volcanic_data.txt", sep = "\t", header = T)
colnames(metabolite_data)[3:4] <- c("FDR", "regulation")
metabolite_data$regulation <- gsub("up", "Up", metabolite_data$regulation)
metabolite_data$regulation <- gsub("down", "Down", metabolite_data$regulation)
metabolite_data$regulation <- gsub("none", "None", metabolite_data$regulation)

metabolite_data$logFDR <- -1*log10(metabolite_data$FDR)

vol <- ggplot(metabolite_data, aes(x = log2(ratio), y = -1*log10(FDR), colour = regulation)) + geom_point(size = 0.2) + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect(fill = NA)) + scale_colour_manual(values = c("#88ABDA", "#9EBC19", "#C35C6A")) +
  ylab("-1*log10(FDR)") + xlab("log2(Fold change)") + xlim(-8,8)

pdf("../Figures/Fig4/volcanic.pdf", width = 5.8 /2.54, height = 4.8/2.54)
vol
graphics.off()

# heatmap
diff_data <- read.table("../Files/pheatmap_reordered.txt", header=T, row.names= 1, sep="\t", comment.char="")
annotation <- read.table("../Files/annotation.txt", header=T, row.names= 1, sep="\t", comment.char="")
group <- read.table("../Files/annotation_group.txt", header=T, row.names= 1, sep="\t", comment.char="")

ann_colors = list(Regulated = c(Up = "#D95F02", Down = "#1B9E77"), Group = c(Control = "#89c3eb", SAR = "#c97586"))

plot <- pheatmap(ln(diff_data),cluster_cols = F ,cluster_rows = F,
                 scale = "row", cluster_col_slices = FALSE,
                 show_rownames = F, show_colnames = F, annotation_row = annotation, 
                 annotation_colors = ann_colors,
                 annotation_col = group, color=colorRampPalette(c("navy", "white","red" ))(6))


plot

# sub heatmap

sub_anno <- subset(annotation, Class == c("Lipids and lipid-like molecules", "Organic acids and derivatives"))
sub_metabolites_names <- data.frame(rownames(sub_anno))
colnames(sub_metabolites_names) <- "ID"

sub_diff_data <- read.table("../Files/pheatmap_reordered.txt", header=T, sep="\t", comment.char="")
sub_metabolites <- left_join(sub_metabolites_names, sub_diff_data, by = "ID")
rownames(sub_metabolites) <- sub_metabolites$ID
sub_metabolites$ID <- NULL

annotation$Class <- gsub("Lipids and lipid-like molecules", "LD", annotation$Class)
annotation$Class <- gsub("Organic acids and derivatives", "OD", annotation$Class)

colnames(annotation)[2] <- "Regulation"

sub_plot <- pheatmap(ln(sub_metabolites), cluster_cols = F ,cluster_rows = F,
                     scale = "row", cluster_col_slices = FALSE,
                     show_rownames = F, show_colnames = F, annotation_row = annotation, 
                     annotation_colors = ann_colors, fontsize = 8, 
                     annotation_col = group, color=colorRampPalette(c("navy", "white","red" ))(6))

pdf("../Figures/Fig3/pheatmap.pdf", width = 11.2/2.54, height = 6/2.54)
sub_plot
graphics.off()

# BCAA
group <- read.table("../Files/annotation_group.txt", header=T, row.names= 1, sep="\t", comment.char="")
diff_aa_anno <- read.table("../Files/amiacid_anno_ordered.txt", sep = "\t", header = T)

diff_aa <- subset(diff_data, rownames(diff_data) %in% c('M245T217','M229T199','M217T127','M247T92','M243T191','M205T81','M246T179','M249T188','M269T78','M246T54',
                                                        'M304T220','M358T240','M187T170'), drop = T)
diff_aa$metabolites <- rownames(diff_aa)
rownames(diff_aa) <- NULL

diff_aa_plot <- melt(diff_aa, id = "metabolites")
colnames(diff_aa_plot)[2] <- "SampleID"

group$SampleID <- rownames(group)
rownames(group) <- NULL

diff_aa_plot <- left_join(diff_aa_plot, group, by = "SampleID")

diff_aa_anno <- read.table("../Files/amiacid_anno_ordered.txt", sep = "\t", header = T)
colnames(diff_aa_plot)[1] <- "ID"
diff_aa_plot <- left_join(diff_aa_plot, diff_aa_anno, by = "ID")

e_plot <- ggplot(diff_aa_plot, aes(x = factor(metabolite, levels = c('Gly-Leu','Leu-Leu-Leu', 'Tryptophyl-Valine','Lysyl-Valine',
                                                                     'Histidinyl-Isoleucine', 'Methionyl-Valine','Glutaminylvaline',
                                                                     'Serylvaline','gamma-Glutamylisoleucine','gamma-Glutamylvaline',
                                                                     'Valyl-Valine','Isoleucylproline','Isoleucyl-Isoleucine')), 
                                   y = log(value), fill = Group)) +
  geom_boxplot(outlier.shape =  NA, notch = F, outlier.size = NA, size = 0.3)  + xlab("Branched-chain amino acids and derivatives") + ylab("log (Intensity)") +
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  coord_flip()

e_plot

pdf("../Figures/Fig3/diff_BCCA_no.pdf", width = 8.6 /2.54, height = 6/2.54)
e_plot
graphics.off()















