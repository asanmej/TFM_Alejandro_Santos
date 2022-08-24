#' 
#' title: "Proteomic_time-course_hypelipidemic_diets_experiment"
#' author: "Alejandro Santos Mejías"
#' date: "11/4/2022"
#' This scrypt have the aim to analysis de the proteomic data from dietary intervation
#' to determine HDL proteome associated to hyperlipidemia. With that in mind, 
#' the first step was to explore and preprocess the dataset, then, a differential
#' expression analysis was performed with its corresponding result visualization,
#' and lastly, an expression pattern analysis. 

library(dplyr)
library(ggplot2)

# Import the data:
#####
dietas <- read.table("./data/Log2_LFQ_intensity_ProteinGroups.txt", header = T, sep = "\t", na.strings = "NaN")
labels <- read.table("./data/Codigo_dietas.txt", sep = "\t", header = T)

#' Data exploration:
# Extract the labels and set properly the subjects' names:
colnames(dietas)
labels <- labels[,-1]
labels
colnames(dietas)[c(1:33)] <- paste("LFQ_",labels, sep = "")
colnames(dietas)
write.table(dietas[,-grep("^LFQ_C38",colnames(dietas))], file = "dietas_dataset.txt", quote = F, sep = "\t", row.names = F)

# Weird subject missing:
head(dietas)
length(grep("^C37",labels, value = T))
length(grep("^C38",labels, value = T))
length(grep("^C39",labels, value = T))
# 3 missing conditions in subject 38 ????????
grep("^C38",labels, value = T) # N diet missing


#' Cleaning the data from Potential contaminant, Only identified by sites, subject 38 
#' and N diet samples:
dietas <- dietas[-which(dietas$Potential.contaminant == "+"),]
dietas <- dietas[-which(dietas$Only.identified.by.site == "+"),]

#' Experimental design matrix, it will be useful for filtering:
groups <- gsub(pattern = "^C3[[:digit:]]", replacement = "", labels)
groups
time <- gsub("[SMNW3]+", "", groups)
time <- gsub("PP","6", time)
time <- gsub("P","3", time)
time <- gsub("B","0", time)
time
condition <- gsub("[PB]+", "", groups)
condition
subject <- gsub(pattern = "[^C3[:digit:]]", replacement = "", labels)
subject <- gsub("3$", "", subject)
subject
metadata <- data.frame(paste("LFQ_",labels, sep = ""), condition, time, subject)
colnames(metadata) <- c("Sample", "Diet", "Time", "Subject")
metadata
metadata <- metadata[order(metadata$Diet, metadata$Time),]
metadata <- metadata[metadata$Subject!="C38",]
metadata <- metadata[order(metadata$Time),]
metadata <- metadata[metadata$Diet != "N",]
write.table(metadata, file = "./data/metadata.txt", quote = F, sep = "\t", row.names = F)

# Delete subject 38, as it has many missing samples and brings noise to the dataset.
dietas <- dietas[,-grep("^LFQ_C38", colnames(dietas))]
newgroups <- dietas[,grep("LFQ_", colnames(dietas))]
newgroups <- newgroups[,metadata$Sample]
rownames(newgroups) <- dietas$Majority.protein.IDs

# After the tidying up the data, the missing values have to be dealt with. Only protein
# with at least 70% of quantified values will be taken into consideration. After that
# this filtering, the remaining NA will be imputed using the minimum value for each 
# protein, there are not many NA in our data after the filtering.
temp <- data.frame(Samples=colnames(newgroups), NA_number=colSums(is.na(newgroups)))

ggplot(temp, aes(x=NA_number, y= Samples)) + 
  geom_bar(stat = "identity", fill = ifelse (temp$NA_number >= 83 ,
                                             "firebrick","steelblue")) +
  theme_minimal() 

# In red, all the samples with more than 30% of missing values. Compare it afterwards
# the removal of those noisy proteins:
newgroups <- newgroups[which(rowMeans(!is.na(newgroups))> 0.70),]
temp <- data.frame(Samples=colnames(newgroups), NA_number=colSums(is.na(newgroups)))
ggplot(temp, aes(x=NA_number, y= Samples)) + 
  geom_bar(stat = "identity", fill = ifelse (temp$NA_number >= 83 ,
                                             "firebrick","steelblue")) +
  theme_minimal() 
temp <- dietas[match(rownames(newgroups), dietas$Majority.protein.IDs),]
write.table(temp, file = "data/dietas_preproccessed.txt", sep = "\t", quote = F, row.names = F)

# In general, there are a low quantity of missing values in each sample, 
# the maximum is 28 out of 119 proteins (23% of NA values) in LFQ_C39PS. 

# Normalization of the dataset:
# the distribution of the data should be checked, before and after the normalization.
# A multiple Shapiro Wilkinson test will be performed to each raw protein intensity in 
# order to check the normality, the p values obtained will be adjusted by an 
# stringent method (Holm method 1979
counts <- 2^newgroups
p_values <- data.frame(protein.id=rownames(counts), 
                       p.value=apply( counts,1, function(x) shapiro.test(x)$p.value))
p_values$adj.p.val <- p.adjust(p_values$p.value, method = "holm", )
table(p_values$adj.p.val>0.05)

library(NormalyzerDE)

design <- data.frame(sample = metadata$Sample, group = metadata$Time)
df_normal <- data.frame(protein.id=rownames(counts), counts)

# It is required raw data without log transformed as it was in newgroups.

write.table(x = design, file = "./data/dietas_design_matrix.tsv",quote = F, row.names = F,
            sep = "\t")
write.table(x = df_normal, file = "./data/dietas_intesity_matrix.tsv",quote = F, row.names = F,
            sep = "\t")

# This section is commented in order to save time while modifying this scrypt:
# NormalyzerDE::normalyzer(jobName = "normalizacion_newgroups",
#                          designPath = "./data/dietas_design_matrix.tsv",
#                          dataPath = "./data/dietas_intesity_matrix.tsv",
#                          outputDir = ".",
#                          requireReplicates = F,
#                          skipAnalysis = F)

# After checking the normalization report, the possible methods are Variance 
# Stabilizing Normalization (VSN) or Quantile normalization for our data. VSN was 
# the chosen one, based in bibliography, even though both are valid options. Also,
# I check the normality graphically:
library(Hmisc)

df_normal <- read.delim("normalizacion_newgroups/VSN-normalized.txt", header = T , sep = "\t")
rownames(df_normal) <- df_normal[,1]
df_normal <- df_normal[,-1]

options(repr.plot.width=100, repr.plot.height=100, res=100)
hist.data.frame(data.frame(t(df_normal)[,sample(1:119, 16)]),nclass=8)

# Then, we calculate the new p values for the normalized dataset:
p_values$norm.p.val <- apply(df_normal,1, function(x) shapiro.test(x)$p.value)
p_values$adj.norm.p.val <- p.adjust(p_values$norm.p.val, method = "holm")

table(p_values$adj.p.val<0.05)
table(p_values$adj.norm.p.val<0.05)

# For the imputation, the minimum value per protein will be taken:
message("Number of NA before imputation")
sum(is.na(df_normal))

for (i in rownames(df_normal)){
  df_normal[i,] <- impute(df_normal[i,], fun = min)
}

message("Number of NA after imputation:")
sum(is.na(df_normal))

# Save the dataset:
gene.name <- dietas$Gene.names[!is.na(match(dietas$Majority.protein.IDs,rownames(df_normal)))]
write.table(x = data.frame(protein.id= rownames(df_normal), gene.name= gene.name, df_normal), file = "./data/dietas_norm_imputed.tsv",quote = F, row.names = F, sep = "\t")

######################################
## Differential expression analysis ##
######################################

#####################
### AMICA Results ###
#####################
# Read the result from the differential analysis performed in AMICA:
amica_res <- read.csv("amica/multiple_comparisons_norm_imputed.csv")

# Volcano plot visualization:

data <- read.delim(file = "amica/amica_protein_groups.tsv")

colnames(data)

ploty <- data %>% select(Gene.names | starts_with("logFC") | starts_with("P.Value"))

ploty$Labels <- NA
ploty$Labels[69] <- ploty$Gene.names[69]
graph <- ggplot(ploty, aes(x = logFC_B__vs__PP, y = -log10(P.Value_B__vs__PP))) + 
  geom_point(stat = "identity", size = 3,
             colour = ifelse((ploty$logFC_B__vs__PP > 1 | ploty$logFC_B__vs__PP < -1) & 
                               ploty$P.Value_B__vs__PP < 0.05 , "#fc8d62", "#66c2a5")) +
  geom_text(label = ploty$Labels, nudge_x = 0.15 , size = 5, colour = "#fc8d62") +
  labs(title = " B vs PP volcano plot", x = "Log2FC (B/PP)", y = "-Log10(P value)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 50), 
        axis.title = element_text(size = 20) )

ggsave("Volcano_plot/Figure 2.jpeg", width = 30, height = 20, units = "cm",dpi = "retina" )

ploty$Labels <- NA
ploty$Labels[which((ploty$logFC_P__vs__PP > 1 | ploty$logFC_P__vs__PP < -1) & 
                     ploty$P.Value_P__vs__PP < 0.05)] <- ploty$Gene.names[which((ploty$logFC_P__vs__PP > 1 | ploty$logFC_P__vs__PP < -1) & 
                                                                                  ploty$P.Value_P__vs__PP < 0.05)]
graph <- ggplot(ploty, aes(x = logFC_P__vs__PP, y = -log10(P.Value_P__vs__PP))) + 
  geom_point(stat = "identity", size = 3,
             colour = ifelse(!is.na(ploty$Labels) , "#fc8d62", "#66c2a5")) +
  geom_text(label = ploty$Labels, nudge_x = 0.15 , size = 5, colour = "#fc8d62") +
  labs(title = " P vs PP volcano plot", x = "Log2FC (P/PP)", y = "-Log10(P value)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 50), 
        axis.title = element_text(size = 20) )

graph

ggsave("Volcano_plot/Figure 3.jpeg", width = 30, height = 20, units = "cm",dpi = "retina" )


######################
### Metaboanalyst  ###
######################
metaboanalyst <- data.frame(Protein.id=rownames(df_normal), df_normal)

B_P <- rbind(c("Label",rep(-1,6),rep(1,6)), metaboanalyst[,1:13])
B_PP <- rbind(c("Label",rep(-1,6),rep(1,6)), metaboanalyst[,c(1:7,14:19)])
P_PP <- rbind(c("Label",rep(-1,6),rep(1,6)), metaboanalyst[,c(1,8:19)])
ANOVA <- rbind(c("Label",metadata$Time), metaboanalyst)

write.table(B_P, file = "metaboanalyst/newgroups_0vs3.txt", quote = F, sep = "\t", row.names = F)
write.table(B_PP, file = "metaboanalyst/newgroups_0vs6.txt", quote = F, sep = "\t", row.names = F)
write.table(P_PP, file = "metaboanalyst/newgroups_0vs3.txt", quote = F, sep = "\t", row.names = F)
write.table(ANOVA, file = "metaboanalyst/newgroups_ANOVA.txt", quote = F, sep = "\t", row.names = F)

# No significant proteins from statistical testing

####################################
## Functional enrichment analysis ##
####################################
library(gprofiler2)

enrich <- gost(amica_res$Gene.names, organism = "hsapiens", sources = c("GO:CC", "GO:MF", "GO:BP", "KEGG", "REAC"), 
               correction_method = "fdr", user_threshold = 0.01)
gostplot(enrich)

temp <- enrich$result
plotemp <- temp %>% dplyr::select(source, term_name, intersection_size, p_value)  
go_plot <- plotemp %>% filter(source == "GO:BP" | source == "GO:CC" | source == "GO:MF") %>%
  arrange(desc(intersection_size), p_value)
cc <- go_plot %>% filter(source == "GO:CC")
mf <- go_plot %>% filter(source == "GO:MF")
bp <- go_plot %>% filter(source == "GO:BP")

pathways <- plotemp %>% filter(source == "KEGG" | source == "REAC") %>% 
  arrange(desc(intersection_size), p_value)

# Functional enrichment result visualization
# GO CC 
ggplot(cc, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(cc$intersection_size>2, "firebrick","dodgerblue"), width = 0.4) +  
  labs(x="Numbers of genes", title = "GO CC Function") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))

# GO MF
ggplot(mf, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(mf$intersection_size>2, "firebrick","dodgerblue"), width = 0.4) +  
  labs(x="Numbers of genes", title = "GO MF Function") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))

# GO BP
bp <- bp[1:20,]
ggplot(bp, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(bp$intersection_size>2, "firebrick","dodgerblue"), width = 0.4) +  
  labs(x="Numbers of genes", title = "GO BP Function") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))

# KEGG and REACTOME enriched pathways
pathways <- pathways[-14,]
ggplot(pathways, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(pathways$intersection_size>2, "firebrick","dodgerblue"), width = 0.4) +  
  labs(x="Numbers of genes", title = "Enriched Pathways") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))

###########################
## HDL PROTEOME WATCH DB ##
###########################
# 
# It was checked if significant proteins were already recorded in HDL Proteome Watch DB,
# and how many times they have been reported in HDL proteomic researches, only
# proteins with appearance in 3 or more different labs were considered HDL associated proteins
# https://homepages.uc.edu/~davidswm/HDLproteome.html

HDL <- read.delim("data/HDL_proteome_hits.tsv")

HDL$Gene.name <- gsub("_HUMAN", "", HDL$Gene.name)

accesion <- dietas$Majority.protein.IDs[match(amica_res$Gene.names, dietas$Gene.names)]
accesion<- unlist(strsplit(accesion, split = ";"))

temp <- accesion[accesion %in% HDL$Accesion]
temp
HDL[HDL$Accesion %in% temp,]

#####################################
## Time-course clustering analysis ##
#####################################
# 
# We used coseq R package to identify expression pattern in our samples depending on 
# the time.
library(coseq)
library(Biobase)
library(corrplot)

# It requires the unormalized dataset without missing values so imputation step needed:
for (i in rownames(counts)){
  counts[i,] <- impute(counts[i,], fun = min)
}

# For the clustering process recommendations from the vignettes were followed
runLogCLR <- coseq(counts, K=2:25, transformation="logclr",norm="TMM", 
                   model="kmeans",
                   nstart=5, iter.max=100, seed = 12345)

# Summary results
summary(runLogCLR)

# Visualization of the results:
# Boxplot per conditions and sample:
plot(runLogCLR, graphs = "boxplots", conds = metadata$Time)
# Profile plots: 
plot(runLogCLR, graphs="profiles")
# Summarized boxplots:
plot(runLogCLR, graphs="boxplots", conds=metadata$Time, collapse_reps = "sum")

# The following code is for saving individually the boxplot graphs:
dev.off(dev.list()["RStudioGD"]) 
library(ggforce)
p <- plot(runLogCLR, graph = "boxplots", conds = metadata$Time)$boxplots
pdf("Clustering analysis/coseq_boxplots_by_page.pdf")
for(k in seq_len(ncol(assay(runLogCLR)))){
  print(p + facet_wrap_paginate( ~labels, page = k, nrow = 1, ncol = 1))
  
}
dev.off()

# Save all the relevant clusters with the proteins in them.
t <- as.data.frame(runLogCLR@allResults$`K=8`)
cluster1 <- rownames(counts[t$Cluster_1 == 1,])
cluster3 <- rownames(counts[t$Cluster_3 == 1,])
cluster5 <- rownames(counts[t$Cluster_5 == 1,])
clusterD <- rownames(counts[t$Cluster_2 == 1 | t$Cluster_7 == 1,])

cluster1 <- dietas[match(cluster1, dietas$Majority.protein.IDs), "Gene.names"] 
cluster3 <- dietas[match(cluster3, dietas$Majority.protein.IDs), "Gene.names"]
cluster5 <- dietas[match(cluster5, dietas$Majority.protein.IDs), "Gene.names"]
clusterD <- dietas[match(clusterD, dietas$Majority.protein.IDs), "Gene.names"]

clusters <- list(cluster1, cluster3, cluster5, clusterD)
names <- c("Cluster_1", "Cluster_3", "Cluster_5", "Cluster_D")

for(i in 1:length(clusters)){
  write.table(clusters[i], file = paste("Clustering analysis/", names[i], ".txt", sep = ""), quote = F, 
              sep = "\t", row.names = F, col.names = F  )
}

# Pathways enrichment analysis of every relevant clusters:

# Cluster 1 #
enrich <- gost(cluster1, organism = "hsapiens", sources = c("KEGG", "REAC"), correction_method = "fdr", user_threshold = 0.01)

temp <- enrich$result
plotemp <- temp %>% dplyr::select(source, term_name, intersection_size, p_value)  
enrich <- plotemp %>% arrange(desc(intersection_size), p_value)
# Terms with long names are removed manually:
enrich <- enrich[-8,]
enrich$intersection_size <- enrich$intersection_size/length(cluster1)

ggplot(enrich, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(enrich$intersection_size>0.5, "firebrick","dodgerblue"), width = 0.4) +  labs(x="Numbers of genes", title = "Cluster 1 Pathways") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))
ggsave("Clustering analysis/Cluster_1_Pathways.png", width = 30, height = 20, units = "cm",dpi = "retina")

# Cluster 3 # 
enrich <- gost(cluster3, organism = "hsapiens", sources = c("KEGG", "REAC"), correction_method = "fdr", user_threshold = 0.01)

temp <- enrich$result
plotemp <- temp %>% dplyr::select(source, term_name, intersection_size, p_value)  
enrich <- plotemp %>% arrange(desc(intersection_size), p_value)
enrich <- enrich[-7,]
enrich$intersection_size <- enrich$intersection_size/length(cluster3)


ggplot(enrich, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(enrich$intersection_size>0.5, "firebrick","dodgerblue"), width = 0.4) +  labs(x="Numbers of genes", title = "Cluster 3 Pathways") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))
ggsave("Clustering analysis/Cluster_3_Pathways.png", width = 30, height = 20, units = "cm",dpi = "retina")

# Cluster 5 #
enrich <- gost(cluster5, organism = "hsapiens", sources = c("KEGG", "REAC"), correction_method = "fdr", user_threshold = 0.01)

temp <- enrich$result
plotemp <- temp %>% dplyr::select(source, term_name, intersection_size, p_value)  
enrich <- plotemp %>% arrange(desc(intersection_size), p_value)
enrich <- enrich[1:22,]
enrich <- enrich[-c(10,19),]
enrich$intersection_size <- enrich$intersection_size/length(cluster5)

ggplot(enrich, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(enrich$intersection_size>0.5, "firebrick","dodgerblue"), width = 0.4) +  labs(x="Numbers of genes", title = "Cluster 5 Pathways") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))
ggsave("Clustering analysis/Cluster_5_Pathways.png", width = 30, height = 20, units = "cm",dpi = "retina")

# Cluster 2 and 7
enrich <- gost(clusterD, organism = "hsapiens", sources = c("KEGG", "REAC"), correction_method = "fdr", user_threshold = 0.01)

temp <- enrich$result
plotemp <- temp %>% dplyr::select(source, term_name, intersection_size, p_value)  
enrich <- plotemp %>% arrange(desc(intersection_size), p_value)
enrich <- enrich[-c(4,21),]
enrich$intersection_size <- enrich$intersection_size/length(clusterD)

ggplot(enrich, aes(x = intersection_size, y = (reorder(term_name, +intersection_size)))) + 
  geom_bar(stat="identity", fill = ifelse(enrich$intersection_size>0.5, "firebrick","dodgerblue"), width = 0.4) +  labs(x="Numbers of genes", title = "Cluster D Pathways") +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(), 
        plot.title = element_text(face = "bold", size = 50, hjust = 0.5))
ggsave("Clustering analysis/Cluster_D_Pathways.png", width = 30, height = 20, units = "cm",dpi = "retina")

