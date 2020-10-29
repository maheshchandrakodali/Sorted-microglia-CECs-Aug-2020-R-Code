
library(tximport)
library(tximportData)
library(apeglm)
library(AnnotationHub)
library(ensembldb)
library(stringr)
library(dplyr)
library(DESeq2)
library(tibble)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(pheatmap)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(ggupset)
library(org.Mm.eg.db)
library(tidyverse)
library(DOSE)
library(pathview)
library(SPIA)


install.packages("ggupset")
update.packages()

# For creating Tx2gene file

ah <- AnnotationHub()
mus_ens <- query(ah, c("Mus musculus", "EnsDb"))
mus_ens
mus_ens <- mus_ens[["AH78811"]]
genes(mus_ens, return.type = "data.frame") %>% View()

# Create a transcript dataframe

txdbb <- transcripts(mus_ens, return.type = "data.frame") %>% dplyr::select(tx_id, gene_id)

txdbb <- txdbb[grep("ENSMUST", txdbb$tx_id),]

# Create a gene-level dataframe
genedb <- genes(mus_ens, return.type = "data.frame")  %>%  dplyr::select(gene_id, entrezid, symbol)

# Merge the two dataframes together
annotations <- inner_join(txdbb, genedb)
annotations %>% View()


# List the file location of the data  

samples <- read.table(file.path("quants", "samples.txt"), header = TRUE)
samples

files <- file.path("C:/Users/Mahesh/Documents/Sorted Cells 6 18 2020 Full/quants", samples$sample, "quant.sf")

files
names(files) <- paste0("N", 1:24)
file.exists(files)
all(file.exists((files)))

#---------------MG only
samplesMG <- read.table(file.path("quantsMG", "samplesMG.txt"), header = TRUE)
samplesMG
filesMG <- file.path("C:/Users/Mahesh/Documents/Sorted Cells 6 18 2020 Full/quantsMG", samplesMG$sample, "quant.sf")

filesMG
names(filesMG) <- paste0("N", 1:12)
file.exists(filesMG)
all(file.exists((filesMG)))

#---------------EC only
samplesEC <- read.table(file.path("quantsEC", "samplesEC.txt"), header = TRUE)
samplesEC
filesEC <- file.path("C:/Users/Mahesh/Documents/Sorted Cells 6 18 2020 Full/quantsEC", samplesEC$sample, "quant.sf")

filesEC
names(filesEC) <- paste0("N", 13:24)
file.exists(filesEC)
all(file.exists((filesEC)))


# Run tximport
txi <- tximport(files, type="salmon", tx2gene=annotations[,c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)

# Run tximport ------------MG only
txiMG <- tximport(filesMG, type="salmon", tx2gene=annotations[,c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)

# Run tximport---------------EC only
txiEC <- tximport(filesEC, type="salmon", tx2gene=annotations[,c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)

attributes(txi)
colnames(txi$counts)

metaa <- data.frame(Cell_Type = factor(c(rep("Microglia", 12), rep("CECs", 12))), Treatment = factor(c(rep("PBS", 3), rep("LPS_30min", 3), rep ("LPS_1hr", 3), rep("LPS_2hrs", 3), rep("PBS", 3), rep("LPS_30min", 3), rep ("LPS_1hr", 3), rep("LPS_2hrs", 3))))

metaa$Cell_Type <- factor(metaa$Cell_Type,levels = c("Microglia", "CECs"))

metaa$Treatment <- factor(metaa$Treatment, levels = c("PBS", "LPS_30min", "LPS_1hr", "LPS_2hrs"))

row.names(metaa) <- colnames(txi$counts)

metaa

all(colnames(txi$counts) %in% rownames(metaa))
all(colnames(txi$counts) == rownames(metaa))

#--------------------MG metadata
metaaMG <- data.frame(Treatment = factor(c(rep("PBS", 3), rep("LPS_30min", 3), rep ("LPS_1hr", 3), rep("LPS_2hrs", 3))))

metaaMG$Treatment <- factor(metaaMG$Treatment, levels = c("PBS", "LPS_30min", "LPS_1hr", "LPS_2hrs"))

row.names(metaaMG) <- colnames(txiMG$counts)

metaaMG

all(colnames(txiMG$counts) %in% rownames(metaaMG))
all(colnames(txiMG$counts) == rownames(metaaMG))

#--------------------EC metadata
metaaEC <- data.frame(Treatment = factor(c(rep("PBS", 3), rep("LPS_30min", 3), rep ("LPS_1hr", 3), rep("LPS_2hrs", 3))))

metaaEC$Treatment <- factor(metaaEC$Treatment, levels = c("PBS", "LPS_30min", "LPS_1hr", "LPS_2hrs"))


row.names(metaaEC) <- colnames(txiEC$counts)

metaaEC

all(colnames(txiEC$counts) %in% rownames(metaaEC))
all(colnames(txiEC$counts) == rownames(metaaEC))


#-------------------------DeSeq for combined
head(txi)

dds <- DESeqDataSetFromTximport(txi, colData = metaa, design = ~ Cell_Type + Treatment + Cell_Type:Treatment)

vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = c("Cell_Type", "Treatment"))

pcaData <- plotPCA(vsd, intgroup = c("Cell_Type", "Treatment"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x=PC1, y=PC2, color = Treatment, shape = Cell_Type)) + 
  geom_point(size = 2) + geom_point(alpha = 1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_bw() + scale_colour_manual(values = c("#FF8C00", "#00ff00", "#FF00FF", "#CCCC00")) + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + ggtitle("")

#-------LRT

dds_LRT <- DESeq(dds, test = "LRT", reduced = ~ Cell_Type + Treatment)

metaa_relevel_for_EC <- metaa

metaa_relevel_for_EC$Cell_Type <- factor(metaa_relevel_for_EC$Cell_Type, levels = c("CECs", "Microglia"))

metaa_relevel_for_EC$Cell_Type

dds_relevel_for_EC<- DESeqDataSetFromTximport(txi, colData = metaa_relevel_for_EC, design = ~ Cell_Type + Treatment + Cell_Type:Treatment)

dds_relevel_for_EC_LRT <- DESeq(dds_relevel_for_EC, test = "LRT", reduced = ~ Cell_Type + Treatment)

#---Normalised counts -- to be used for plotting heatmaps

normalized_counts <- counts(dds_LRT, normalized=T) %>%  data.frame() %>%  rownames_to_column(var="gene") 

normalized_counts <- normalized_counts[,c(1:25)] 

#--------annotations

grch99annot <- annotations %>% 
  dplyr::select(gene_id, symbol) %>% 
  dplyr::distinct()

normalized_counts <- merge(normalized_counts, grch99annot, by.x="gene", by.y="gene_id")

normalized_counts <- normalized_counts %>%
  as_tibble()

#--- Below results give only the genes different in expression from microglia at a given time point, if they are expressed similar to Microglia, then they arent included**

library(DESeq2)

res_dds_LRT_30min <- results(dds_LRT, name = "Cell_TypeCECs.TreatmentLPS_30min", test = "Wald", alpha = 0.05)

res_dds_LRT_30min_DF <- as.data.frame(res_dds_LRT_30min)

res_dds_LRT_30min_DF$significant <- ifelse(res_dds_LRT_30min_DF$padj < .05 & abs(res_dds_LRT_30min_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_30min_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change")  + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (90)", "Not significant")) + theme_bw() + ggtitle("LPS 30min CECs vs LPS 30min Microglia") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_1hr <- results(dds_LRT, name = "Cell_TypeCECs.TreatmentLPS_1hr", test = "Wald", alpha = 0.05)

res_dds_LRT_1hr_DF <- as.data.frame(res_dds_LRT_1hr)

res_dds_LRT_1hr_DF$significant <- ifelse(res_dds_LRT_1hr_DF$padj < .05 & abs(res_dds_LRT_1hr_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_1hr_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (202)", "Not significant")) + theme_bw() + ggtitle("LPS 1hr CECs vs LPS 1hr Microglia") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_2hrs <- results(dds_LRT, name = "Cell_TypeCECs.TreatmentLPS_2hrs", test = "Wald", alpha = 0.05)


res_dds_LRT_2hrs_DF <- as.data.frame(res_dds_LRT_2hrs)

res_dds_LRT_2hrs_DF$significant <- ifelse(res_dds_LRT_2hrs_DF$padj < .05 & abs(res_dds_LRT_2hrs_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_2hrs_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (2078)", "Not significant")) + theme_bw() + ggtitle("LPS 2hrs CECs vs LPS 2hrs Microglia") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_30min_tb <- res_dds_LRT_30min %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_dds_LRT_1hr_tb <- res_dds_LRT_1hr %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_dds_LRT_2hrs_tb <- res_dds_LRT_2hrs %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

padj.cutoff <- 0.05

sig_res_dds_LRT_30min <- res_dds_LRT_30min_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_1hr <- res_dds_LRT_1hr_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_2hrs <- res_dds_LRT_2hrs_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

nrow(sig_res_dds_LRT_30min)
nrow(sig_res_dds_LRT_1hr)
nrow(sig_res_dds_LRT_2hrs)


#---extracting the sig genes from the LRT test, for the Treatments  compared to control i.e. Microglia (was the reference Factor Level)

res_dds_LRT_30min_MGonly <- results(dds_LRT, name = "Treatment_LPS_30min_vs_PBS", test = "Wald", alpha = 0.05)

res_dds_LRT_1hr_MGonly <- results(dds_LRT, name = "Treatment_LPS_1hr_vs_PBS", test = "Wald", alpha=0.05)

res_dds_LRT_2hrs_MGonly <- results(dds_LRT, name = "Treatment_LPS_2hrs_vs_PBS", test = "Wald", alpha=0.05)


res_dds_LRT_30min_MGonly_DF <- as.data.frame(res_dds_LRT_30min_MGonly)

res_dds_LRT_30min_MGonly_DF$significant <- ifelse(res_dds_LRT_30min_MGonly_DF$padj < .05 & abs(res_dds_LRT_30min_MGonly_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_30min_MGonly_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (17)", "Not significant")) + theme_bw() + ggtitle("LPS 30min vs PBS - Microglia") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_1hr_MGonly_DF <- as.data.frame(res_dds_LRT_1hr_MGonly)

res_dds_LRT_1hr_MGonly_DF$significant <- ifelse(res_dds_LRT_1hr_MGonly_DF$padj < .05 & abs(res_dds_LRT_1hr_MGonly_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_1hr_MGonly_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (58)", "Not significant")) + theme_bw() + ggtitle("LPS 1hr vs PBS - Microglia") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")

res_dds_LRT_2hrs_MGonly_DF <- as.data.frame(res_dds_LRT_2hrs_MGonly)

res_dds_LRT_2hrs_MGonly_DF$significant <- ifelse(res_dds_LRT_2hrs_MGonly_DF$padj < .05 & abs(res_dds_LRT_2hrs_MGonly_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_2hrs_MGonly_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (601)", "Not significant")) + theme_bw() + ggtitle("LPS 2hrs vs PBS - Microglia") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_30min_MGonly_tb <- res_dds_LRT_30min_MGonly %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_dds_LRT_1hr_MGonly_tb <- res_dds_LRT_1hr_MGonly %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_dds_LRT_2hrs_MGonly_tb <- res_dds_LRT_2hrs_MGonly %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_res_dds_LRT_30min_MGonly <- res_dds_LRT_30min_MGonly_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_30min_MGonly_Txt <- merge(sig_res_dds_LRT_30min_MGonly, grch99annot, by.x="gene", by.y="gene_id")

sig_res_dds_LRT_1hr_MGonly <- res_dds_LRT_1hr_MGonly_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_1hr_MGonly_Txt <- merge(sig_res_dds_LRT_1hr_MGonly, grch99annot, by.x="gene", by.y="gene_id")

sig_res_dds_LRT_2hrs_MGonly <- res_dds_LRT_2hrs_MGonly_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_2hrs_MGonly_Txt <- merge(sig_res_dds_LRT_2hrs_MGonly, grch99annot, by.x="gene", by.y="gene_id")

nrow(sig_res_dds_LRT_30min_MGonly)
nrow(sig_res_dds_LRT_1hr_MGonly)
nrow(sig_res_dds_LRT_2hrs_MGonly)



#---extracting the sig genes from the LRT test, for the Treatments  compared to control, for  Endothelial Cell (EC is relevelled as the reference Factor Level, to get the sig genes for treatment)


res_relevel_for_EC_LRT <- results(dds_relevel_for_EC_LRT)


resultsNames(dds_relevel_for_EC_LRT)


res_dds_LRT_30min_EConly <- results(dds_relevel_for_EC_LRT, name = "Treatment_LPS_30min_vs_PBS", test = "Wald", alpha = 0.05)


res_dds_LRT_30min_EConly_DF <- as.data.frame(res_dds_LRT_30min_EConly)

res_dds_LRT_30min_EConly_DF$significant <- ifelse(res_dds_LRT_30min_EConly_DF$padj < .05 & abs(res_dds_LRT_30min_EConly_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_30min_EConly_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (318)", "Not significant")) + theme_bw() + ggtitle("LPS 30min vs PBS - CECs") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_1hr_EConly <- results(dds_relevel_for_EC_LRT, name = "Treatment_LPS_1hr_vs_PBS", test = "Wald", alpha = 0.05)

res_dds_LRT_1hr_EConly_DF <- as.data.frame(res_dds_LRT_1hr_EConly)

res_dds_LRT_1hr_EConly_DF$significant <- ifelse(res_dds_LRT_1hr_EConly_DF$padj < .05 & abs(res_dds_LRT_1hr_EConly_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_1hr_EConly_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (663)", "Not significant")) + theme_bw() + ggtitle("LPS 1hr vs PBS - CECs") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_dds_LRT_2hrs_EConly <- results(dds_relevel_for_EC_LRT, name = "Treatment_LPS_2hrs_vs_PBS", test = "Wald", alpha = 0.05)

res_dds_LRT_2hrs_EConly_DF <- as.data.frame(res_dds_LRT_2hrs_EConly)

res_dds_LRT_2hrs_EConly_DF$significant <- ifelse(res_dds_LRT_2hrs_EConly_DF$padj < .05 & abs(res_dds_LRT_2hrs_EConly_DF$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_dds_LRT_2hrs_EConly_DF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (4244)", "Not significant")) + theme_bw() + ggtitle("LPS 2hrs vs PBS - CECs") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")



res_dds_LRT_30min_EConly_tb <- res_dds_LRT_30min_EConly %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_dds_LRT_1hr_EConly_tb <- res_dds_LRT_1hr_EConly %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_dds_LRT_2hrs_EConly_tb <- res_dds_LRT_2hrs_EConly %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


sig_res_dds_LRT_30min_EConly <- res_dds_LRT_30min_EConly_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_30min_EConly_Txt <- merge(sig_res_dds_LRT_30min_EConly, grch99annot, by.x="gene", by.y="gene_id")

sig_res_dds_LRT_1hr_EConly <- res_dds_LRT_1hr_EConly_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_1hr_EConly_Txt <- merge(sig_res_dds_LRT_1hr_EConly, grch99annot, by.x="gene", by.y="gene_id")

sig_res_dds_LRT_2hrs_EConly <- res_dds_LRT_2hrs_EConly_tb %>%  filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig_res_dds_LRT_2hrs_EConly_Txt <- merge(sig_res_dds_LRT_2hrs_EConly, grch99annot, by.x="gene", by.y="gene_id")

nrow(sig_res_dds_LRT_30min_EConly)
nrow(sig_res_dds_LRT_1hr_EConly)
nrow(sig_res_dds_LRT_2hrs_EConly)


#-------------------------------

#---Normalised counts -- to be used for plotting heatmaps

normalized_counts <- counts(dds_LRT, normalized=T) %>%  data.frame() %>%  rownames_to_column(var="gene") 

normalized_counts <- normalized_counts[,c(1:25)] 

#--------annotations

grch99annot <- annotations %>% 
  dplyr::select(gene_id, symbol) %>% 
  dplyr::distinct()

normalized_counts <- merge(normalized_counts, grch99annot, by.x="gene", by.y="gene_id")

normalized_counts <- normalized_counts %>%
  as_tibble()

#--LRT ego

res_ids_LRT_30min <- inner_join(res_dds_LRT_30min_tb, annotations, by=c("gene"="gene_id"))

all_LRT_30min_genes <- as.character(res_ids_LRT_30min$gene)

LRT_30min_genes <- as.character(sig_res_dds_LRT_30min$gene)

ego_LRT_30min <- enrichGO(gene = LRT_30min_genes, universe = all_LRT_30min_genes,
                          keyType = "ENSEMBL",
                          OrgDb = org.Mm.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = TRUE)

cluster_summary_ego_LRT_30min <- data.frame(ego_LRT_30min)
write.csv(cluster_summary_ego_LRT_30min, "data/cluster_summary_ego_LRT_30min.csv")


dotplot(ego_LRT_30min, showCategory=20)

#---1hr LRT

res_ids_LRT_1hr <- inner_join(res_dds_LRT_1hr_tb, annotations, by=c("gene"="gene_id"))

all_LRT_1hr_genes <- as.character(res_ids_LRT_1hr$gene)

LRT_1hr_genes <- as.character(sig_res_dds_LRT_1hr$gene)

ego_LRT_1hr <- enrichGO(gene = LRT_1hr_genes, universe = all_LRT_1hr_genes,
                        keyType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)

cluster_summary_ego_LRT_1hr <- data.frame(ego_LRT_1hr)
write.csv(cluster_summary_ego_LRT_1hr, "data/cluster_summary_ego_LRT_1hr.csv")


dotplot(ego_LRT_1hr, showCategory=50)

#----2hrs LRT

res_ids_LRT_2hrs <- inner_join(res_dds_LRT_2hrs_tb, annotations, by=c("gene"="gene_id"))

all_LRT_2hrs_genes <- as.character(res_ids_LRT_2hrs$gene)

LRT_2hrs_genes <- as.character(sig_res_dds_LRT_2hrs$gene)

ego_LRT_2hrs <- enrichGO(gene = LRT_2hrs_genes, universe = all_LRT_2hrs_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)

cluster_summary_ego_LRT_2hrs <- data.frame(ego_LRT_2hrs)
write.csv(cluster_summary_ego_LRT_2hrs, "data/cluster_summary_ego_LRT_2hrs.csv")


dotplot(ego_LRT_2hrs, showCategory=50)

#---ego for MG only from the LRT test, for the Treatments  compared to control i.e. Microglia (was the reference Factor Level)

#--30min MG only LRT

res_ids_LRT_MGonly30min <- inner_join(res_dds_LRT_30min_MGonly_tb, annotations, by=c("gene"="gene_id"))

all_LRT_MGonly30min_genes <- as.character(res_ids_LRT_MGonly30min$gene)


dds_LRT_30min_MGonly_genes <- as.character(sig_res_dds_LRT_30min_MGonly$gene)

nrow(sig_res_dds_LRT_30min_MGonly)

sig_res_dds_LRT_30min_MGonly

ego_LRT_30min_MGonly <- enrichGO(gene = dds_LRT_30min_MGonly_genes, universe = all_LRT_MGonly30min_genes,
                                 keyType = "ENSEMBL",
                                 OrgDb = org.Mm.eg.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)

cluster_summary_ego_LRT_30min_MGonly <- data.frame(ego_LRT_30min_MGonly)

write.csv(cluster_summary_ego_LRT_30min_MGonly, "data/cluster_summary_ego_LRT_30min_MGonly.csv")

cluster_summary_ego_LRT_30min_MGonly <- data.frame(ego_LRT_30min_MGonly)

dotplot(ego_LRT_30min_MGonly, showCategory=50)


#1hr MG only LRT

res_ids_LRT_MGonly1hr <- inner_join(res_dds_LRT_1hr_MGonly_tb, annotations, by=c("gene"="gene_id"))

all_LRT_MGonly1hr_genes <- as.character(res_ids_LRT_MGonly1hr$gene)

dds_LRT_1hr_MGonly_genes <- as.character(sig_res_dds_LRT_1hr_MGonly$gene)

ego_LRT_1hr_MGonly <- enrichGO(gene = dds_LRT_1hr_MGonly_genes, universe = all_LRT_MGonly1hr_genes,
                               keyType = "ENSEMBL",
                               OrgDb = org.Mm.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = TRUE)

cluster_summary_ego_LRT_1hr_MGonly <- data.frame(ego_LRT_1hr_MGonly)
write.csv(cluster_summary_ego_LRT_1hr_MGonly, "data/cluster_summary_ego_LRT_1hr_MGonly.csv")


dotplot(ego_LRT_1hr_MGonly, showCategory=50)

#--2hr MG only LRT 

res_ids_LRT_MGonly2hrs <- inner_join(res_dds_LRT_2hrs_MGonly_tb, annotations, by=c("gene"="gene_id"))

all_LRT_MGonly2hrs_genes <- as.character(res_ids_LRT_MGonly2hrs$gene)

dds_LRT_2hrs_MGonly_genes <- as.character(sig_res_dds_LRT_2hrs_MGonly$gene)

ego_LRT_2hrs_MGonly <- enrichGO(gene = dds_LRT_2hrs_MGonly_genes, universe = all_LRT_MGonly2hrs_genes,
                                keyType = "ENSEMBL",
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = TRUE)

cluster_summary_ego_LRT_2hrs_MGonly <- data.frame(ego_LRT_2hrs_MGonly)
write.csv(cluster_summary_ego_LRT_2hrs_MGonly, "data/cluster_summary_ego_LRT_2hrs_MGonly.csv")


dotplot(ego_LRT_2hrs_MGonly, showCategory=50)


#---ego for EC only from the LRT test after relevelling the Cell_Type Factor, for the Treatments  compared to control i.e. (**Endothelial_Cells is the reference Factor Level)

#--30min EC only LRT

res_ids_LRT_EConly30min <- inner_join(res_dds_LRT_30min_EConly_tb, annotations, by=c("gene"="gene_id"))

all_LRT_EConly30min_genes <- as.character(res_ids_LRT_EConly30min$gene)

dds_LRT_30min_EConly_genes <- as.character(sig_res_dds_LRT_30min_EConly$gene)

ego_LRT_30min_EConly <- enrichGO(gene = dds_LRT_30min_EConly_genes, universe = all_LRT_EConly30min_genes,
                                 keyType = "ENSEMBL",
                                 OrgDb = org.Mm.eg.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)

cluster_summary_ego_LRT_30min_EConly <- data.frame(ego_LRT_30min_EConly)
write.csv(cluster_summary_ego_LRT_30min_EConly, "data/cluster_summary_ego_LRT_30min_EConly.csv")

dotplot(ego_LRT_30min_EConly, showCategory=12)

#1hr EC only LRT

res_ids_LRT_EConly1hr <- inner_join(res_dds_LRT_1hr_EConly_tb, annotations, by=c("gene"="gene_id"))

all_LRT_EConly1hr_genes <- as.character(res_ids_LRT_EConly1hr$gene)


dds_LRT_1hr_EConly_genes <- as.character(sig_res_dds_LRT_1hr_EConly$gene)

ego_LRT_1hr_EConly <- enrichGO(gene = dds_LRT_1hr_EConly_genes, universe = all_LRT_EConly1hr_genes,
                               keyType = "ENSEMBL",
                               OrgDb = org.Mm.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = TRUE)

cluster_summary_ego_LRT_1hr_EConly <- data.frame(ego_LRT_1hr_EConly)
write.csv(cluster_summary_ego_LRT_1hr_EConly, "data/cluster_summary_ego_LRT_1hr_EConly.csv")


dotplot(ego_LRT_1hr_EConly, showCategory=50)

#--2hr EC only LRT

res_ids_LRT_EConly2hrs <- inner_join(res_dds_LRT_2hrs_EConly_tb, annotations, by=c("gene"="gene_id"))

all_LRT_EConly2hrs_genes <- as.character(res_ids_LRT_EConly2hrs$gene)

dds_LRT_2hrs_EConly_genes <- as.character(sig_res_dds_LRT_2hrs_EConly$gene)

ego_LRT_2hrs_EConly <- enrichGO(gene = dds_LRT_2hrs_EConly_genes, universe = all_LRT_EConly2hrs_genes,
                                keyType = "ENSEMBL",
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = TRUE)

cluster_summary_ego_LRT_2hrs_EConly <- data.frame(ego_LRT_2hrs_EConly)
write.csv(cluster_summary_ego_LRT_2hrs_EConly, "data/cluster_summary_ego_LRT_2hrs_EConly.csv")


dotplot(ego_LRT_2hrs_EConly, showCategory=50)

#---Plotting top inflammation categories

ego_top_EC30min <- ego_LRT_30min_EConly 

ego_top_EC30min <- dplyr::filter(ego_top_EC30min, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_EC30min, showCategory = 10)

ego_top_EC1hr <- ego_LRT_1hr_EConly 

ego_top_EC1hr <- dplyr::filter(ego_top_EC1hr, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_EC1hr, showCategory = 10)


ego_top_EC2hrs <- ego_LRT_2hrs_EConly 

ego_top_EC2hrs <- dplyr::filter(ego_top_EC2hrs, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_EC2hrs, showCategory = 10)

ego_vessel_related_2hrs_EConly <- dplyr::filter(ego_LRT_2hrs_EConly, ID %in% c("GO:2001233", "GO:0034330", "GO:0061028", "GO:0007044", "GO:0034332", "GO:0007045", "GO:0090109", "GO:2000351", "GO:0043297", "GO:0120192", "GO:1903392", "GO:0016264", "GO:2000353", "GO:0045216", "GO:0034329", "GO:0070830"))

dotplot(ego_vessel_related_2hrs_EConly, showCategory=15)

ego_top_MG30min <- ego_LRT_30min_MGonly 

ego_top_MG30min <- dplyr::filter(ego_top_MG30min, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_MG30min, showCategory = 10)

ego_top_MG1hr <- ego_LRT_1hr_MGonly 

ego_top_MG1hr <- dplyr::filter(ego_top_MG1hr, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_MG1hr, showCategory = 10)

ego_top_MG2hrs <- ego_LRT_2hrs_MGonly 

ego_top_MG2hrs <- dplyr::filter(ego_top_MG2hrs, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_MG2hrs, showCategory = 10)

#------ego plot for LRT

ego_top_30min <- ego_LRT_30min 

ego_top_30min <- dplyr::filter(ego_top_30min, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_30min, showCategory = 10)


ego_top_1hr <- ego_LRT_1hr 

ego_top_1hr <- dplyr::filter(ego_top_1hr, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_1hr, showCategory = 10)

ego_top_2hrs <- ego_LRT_2hrs 

ego_top_2hrs <- dplyr::filter(ego_top_2hrs, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_2hrs, showCategory = 10)



#-------For heatmaps

#-----Heatmaps of the selected pathway genes from each GO term, genes significant at each timepoint in ECminus and MG saparately---------------


#------------------------IEGs

IEGs <- dplyr::filter(normalized_counts, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

resp_lps_EC30min <- dplyr::filter(normalized_counts, symbol %in% c("Tnfaip3", "Nfkbia", "Ly86", "Litaf", "Tnf", "Mrc1", "Cxcl1", "Arid5a", "Gch1", "Zfp36", "Cd14", "Cx3cr1", "Jun", "Junb", "Ticam2", "Cxcl2", "Irak2", "Jund"))

nfkb_signalling_EC30min <- dplyr::filter(normalized_counts, symbol %in%c("Irf1", "Tnfaip3", "Nfkbia", "Ajuba", "Litaf", "Tnf", "Rhoh", "Nfkbid", "Icam1", "Cd14", "Cx3cr1", "Ticam2", "Irak2"))


pos_reg_cyt_prdn_EC30min <- dplyr::filter(normalized_counts, symbol %in% c("Brca1", "Irf1", "Irf4", "Aif1", "Tnf", "Csf1r", "Ptprc", "Tyrobp", "Arid5a", "Egr1", "Atf4", "Cd14", "Ticam2", "Fcer1g", "Fcgr3", "Lilra5"))


tnf_prdn_EC30min <- dplyr::filter(normalized_counts, symbol %in% c("Tnfaip3", "Ptprc", "Tyrobp", "Arid5a", "Zfp36", "Cd14", "Cx3cr1", "Fcer1g", "Fcgr3", "Lilra5"))

go_all_four_EC30min <- dplyr::filter(normalized_counts, symbol %in% c("Tnfaip3", "Nfkbia", "Ly86", "Litaf", "Tnf", "Mrc1", "Cxcl1", "Arid5a", "Gch1", "Zfp36", "Cd14", "Cx3cr1", "Jun", "Junb", "Ticam2", "Cxcl2", "Irak2", "Jund", "Irf1", "Ajuba", "Rhoh", "Nfkbid", "Icam1", "Brca1", "Irf4", "Aif1", "Csf1r", "Ptprc", "Tyrobp", "Egr1", "Atf4", "Fcer1g", "Fcgr3", "Lilra5"))

# after deleting no change in heatmap

go_all_four_EC30min <- dplyr::filter(normalized_counts, symbol %in% c("Tnfaip3", "Nfkbia", "Litaf", "Tnf", "Cxcl1", "Arid5a", "Gch1", "Zfp36", "Jun", "Junb", "Cxcl2", "Irak2", "Jund", "Irf1", "Ajuba", "Nfkbid", "Icam1", "Egr1", "Atf4"))

go_all_four_EC30min <- as.data.frame(go_all_four_EC30min)

row.names(go_all_four_EC30min) <- go_all_four_EC30min$symbol



#----------example for extracting the selected results (fold changes, from deseq2 results) into text file
EC30min_Txt <- dplyr::filter(sig_res_dds_LRT_30min_EConly_Txt, symbol %in% c("Tnfaip3", "Nfkbia", "Litaf", "Tnf", "Prdx3", "Bcl10", "Cxcl1", "Arid5a", "Gch1", "Zfp36", "Jun", "Junb", "Cxcl2", "Irak2"))


write.csv(EC30min_Txt, "data/EC30min_Txt.csv")


EC30_vessel_dia <- dplyr::filter(normalized_counts, symbol %in% c("Hpn", "Pmp22", "Ctss", "Irf1", "Tnfaip3", "Nf1", "Dusp1", "Ptprc", "March7", "Pde3b", "Plxnb3", "Nfkbid", "Rnd1", "C1qa", "C1qb", "C1qc", "Icam1", "Gch1", "Dusp5", "Tnf"))


EC30_vessel_dia_Txt <- dplyr::filter(sig_res_dds_LRT_30min_EConly_Txt, symbol %in% c("Hpn", "Pmp22", "Ctss", "Irf1", "Tnfaip3", "Nf1", "Dusp1", "Ptprc", "March7", "Pde3b", "Plxnb3", "Nfkbid", "Rnd1", "C1qa", "C1qb", "C1qc", "Icam1", "Gch1", "Dusp5", "Tnf"))

write.csv(EC30_vessel_dia_Txt, "data/EC30_vessel_dia_Txt.csv")



EC30_vessel_dia <- as.data.frame(EC30_vessel_dia)
row.names(EC30_vessel_dia) <- EC30_vessel_dia$symbol


annoEC30 <- read.csv(file.path("annorow.csv"), header = TRUE, sep = ",", row.names = 1)

annoEC30<-as.data.frame(annoEC30)
annoEC30$GO.0032496 <- factor(annoEC30$GO.0032496, exclude = "")
annoEC30$GO.0051092 <- factor(annoEC30$GO.0051092, exclude = "")
annoEC30$GO.0001819 <- factor(annoEC30$GO.0001819, exclude = "")
annoEC30$GO.0032640 <- factor(annoEC30$GO.0032640, exclude = "")

ann_colorsEC30 = list(
  Treatment = c(PBS = "#FF8C00", LPS_30min = "#00ff00"), Cell_Type = c(Microglia = "#000000", CECs = "#2083FB"), GO.0032496 = c("Response to lipopolysaccharide" = "#008080"), GO.0051092 = c("Positive regulation of NF-kappaB transcription factor activity" = "#F5DEB3"), GO.0001819 = c("Positive regulation of cytokine production" = "#A52A2A"), GO.0032640 = c("Tumor necrosis factor production" = "#FF7F50"))

pheatmap(go_all_four_EC30min[c(2:7, 14:19)],          color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = 6,
         show_rownames = T,
         show_colnames = F,
         treeheight_row = 10,
         annotation = (metaa),
         annotation_row = annoEC30,
         annotation_names_row = F,
         annotation_colors = ann_colorsEC30,
         legend_breaks = c(-2,0,2),
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20,
         cellwidth = 25) 


ann_colorsEC30dia = list(Treatment = c(PBS = "#FF8C00", LPS_30min = "#00ff00"), Cell_Type = c(CECs = "#D95F02"))


annorowEC30dia <- read.csv(file.path("annorowEC30morphology.csv"), header = TRUE, sep = ",", row.names = 1)

annorowEC30dia<-as.data.frame(annorowEC30dia)
annorowEC30dia$GO.0071711 <- factor(annorowEC30dia$GO.0071711, exclude = "")
annorowEC30dia$GO.0007162 <- factor(annorowEC30dia$GO.0007162, exclude = "")
annorowEC30dia$GO.0150146 <- factor(annorowEC30dia$GO.0150146, exclude = "")
annorowEC30dia$GO.0097746 <- factor(annorowEC30dia$GO.0097746, exclude = "")

ann_colorsves30 = list(Treatment = c(PBS = "#FF8C00", LPS_30min = "#00ff00"), Cell_Type = c(CECs = "#2083FB"))

pheatmap(EC30_vessel_dia[c(14:19)],       
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         treeheight_row = 10,
         annotation_row = annorowEC30dia,
         annotation_names_row = F,
         annotation = (metaa),
         legend_breaks = c(-1.5,0,1.5),
         annotation_colors = ann_colorsves30,
         annotation_names_col = F,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20,
         cellwidth = 40) 


#EC 1hr

resp_lps_EC1hr <- dplyr::filter(normalized_counts, symbol %in% c("Tgfb1", "Bcr", "Cxcl16", "Tnfaip3", "Nfkbia", "Litaf", "Cd86", "Noct", "Tnf", "Rela", "Ncf2", "Mrc1", "Bcl10", "Cxcl1", "Cx3cl1", "Ptgs2", "Cxcl10", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Irf8", "Zc3h12a", "Zfp36", "Cd14", "Jun", "Junb", "Nod2", "Cebpb", "Cxcl2", "Cxcl11", "Irak2", "Jund", "Cd80", "Selenos"))

nfkb_signalling_EC1hr <- dplyr::filter(normalized_counts, symbol %in%c("Itgb2", "Tgfb1", "Relb", "Irf1", "Tnfaip3", "Rel", "Nfkbia", "Edn1", "Clu", "Litaf", "Tnf", "Rela", "Nfkb2", "Alpk1", "Bcl10", "Prkcz", "Cx3cl1", "Birc3", "Nfkbid", "Icam1", "Ripk2", "Zc3h12a", "Tifa", "Cd14", "Bcl3", "Nod2", "Irak2"))


pos_reg_cyt_prdn_EC1hr <- dplyr::filter(normalized_counts, symbol %in% c("Icosl", "Ccl3", "Tgfb1", "Hyal2", "Cybb", "Irf1", "Rel", "Clu", "Tnf", "Rela", "Bcl10", "Nr4a3", "Prkcz", "Il4ra", "Cx3cl1", "Ptgs2", "P2ry2", "Tnfsf9", "Serpine1", "Arid5a", "Egr1", "Akap12", "Ripk2", "Irf8", "Atf4", "Ccr2", "Cd14", "Bcl3", "Nod2", "Cebpg", "Cebpb", "Fcer1g", "Lilra5"))

tnf_prdn_EC1hr <- dplyr::filter(normalized_counts, symbol %in% c("Ccl3", "Ptpn6", "Cybb", "Tnfaip3", "Igf1", "Clu", "Errfi1", "Cx3cl1", "Arid5a", "Akap12", "Ripk2", "Zc3h12a", "Zfp36", "Ccr2", "Trex1", "Cd14", "Bcl3", "Nod2", "Fcer1g", "Lilra5", "Selenos"))

go_all_four_EC1hr <- dplyr::filter(normalized_counts, symbol %in% c("Tgfb1", "Bcr", "Cxcl16", "Tnfaip3", "Nfkbia", "Litaf", "Cd86", "Noct", "Tnf", "Rela", "Ncf2", "Mrc1", "Bcl10", "Cxcl1", "Cx3cl1", "Ptgs2", "Cxcl10", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Irf8", "Zc3h12a", "Zfp36", "Cd14", "Jun", "Junb", "Nod2", "Cebpb", "Cxcl2", "Cxcl11", "Irak2", "Jund", "Cd80", "Selenos", "Itgb2", "Relb", "Irf1", "Rel", "Edn1", "Clu",  "Nfkb2", "Alpk1", "Prkcz", "Birc3", "Nfkbid", "Icam1", "Tifa", "Bcl3", "Icosl", "Ccl3", "Hyal2", "Cybb", "Irf1", "Rel", "Clu", "Nr4a3", "Prkcz", "Il4ra", "P2ry2", "Tnfsf9",  "Egr1", "Akap12", "Atf4", "Ccr2", "Bcl3", "Cebpg", "Fcer1g", "Lilra5", "Ptpn6", "Igf1", "Errfi1", "Trex1"))

nrow(go_all_four_EC1hr)

# after deleting no change in heatmap

go_all_four_EC1hr <- dplyr::filter(normalized_counts, symbol %in% c("Tgfb1", "Bcr", "Tnfaip3", "Nfkbia", "Litaf", "Noct", "Tnf", "Rela", "Bcl10", "Cxcl1", "Cx3cl1", "Ptgs2", "Cxcl10", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Zc3h12a", "Zfp36", "Jun", "Junb", "Nod2", "Cxcl2", "Cxcl11", "Irak2", "Jund", "Relb", "Irf1", "Rel", "Edn1", "Clu",  "Nfkb2", "Alpk1", "Prkcz", "Birc3", "Nfkbid", "Icam1", "Tifa", "Bcl3", "Hyal2", "Irf1", "Rel", "Clu", "Nr4a3", "P2ry2", "Egr1", "Akap12", "Atf4", "Bcl3", "Errfi1", "Trex1"))

EC1hr_Txt <- dplyr::filter(sig_res_dds_LRT_1hr_EConly_Txt, symbol %in% c("Tgfb1", "Bcr", "Tnfaip3", "Nfkbia", "Litaf", "Noct", "Tnf", "Rela", "Bcl10", "Cxcl1", "Cx3cl1", "Ptgs2", "Cxcl10", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Zc3h12a", "Zfp36", "Jun", "Junb", "Nod2", "Cxcl2", "Cxcl11", "Irak2", "Jund", "Relb", "Irf1", "Rel", "Edn1", "Clu",  "Nfkb2", "Alpk1", "Prkcz", "Birc3", "Nfkbid", "Icam1", "Tifa", "Bcl3", "Hyal2", "Irf1", "Rel", "Clu", "Nr4a3", "P2ry2", "Egr1", "Akap12", "Atf4", "Bcl3", "Errfi1", "Trex1"))


write.csv(EC1hr_Txt, "data/EC1hr_Txt.csv")


go_all_four_EC1hr <- as.data.frame(go_all_four_EC1hr)

row.names(go_all_four_EC1hr) <- go_all_four_EC1hr$symbol


annorowEC1hr <- read.csv(file.path("annorowEC1hr.csv"), header = TRUE, sep = ",", row.names = 1)

annorowEC1hr<-as.data.frame(annorowEC1hr)
annorowEC1hr$GO.0032496 <- factor(annorowEC1hr$GO.0032496, exclude = "")
annorowEC1hr$GO.0051092 <- factor(annorowEC1hr$GO.0051092, exclude = "")
annorowEC1hr$GO.0001819 <- factor(annorowEC1hr$GO.0001819, exclude = "")
annorowEC1hr$GO.0032640 <- factor(annorowEC1hr$GO.0032640, exclude = "")


ann_colors1hr = list(
  Treatment = c(PBS = "#FF8C00", LPS_30min = "#00ff00", LPS_1hr = "#FF00FF"), Cell_Type = c(Microglia = "#000000", CECs = "#2083FB"), GO.0032496 = c("Response to lipopolysaccharide" = "#008080"), GO.0051092 = c("Positive regulation of NF-kappaB transcription factor activity" = "#F5DEB3"), GO.0001819 = c("Positive regulation of cytokine production" = "#A52A2A"), GO.0032640 = c("Tumor necrosis factor production" = "#FF7F50"))


metaanno <- data.frame(Cell_Type = factor(c(rep("Microglia", 12), rep("CECs", 12))), Treatment = factor(c(rep("PBS", 3), rep("LPS_30min", 3), rep ("LPS_1hr", 3), rep("LPS_2hrs", 3), rep("PBS", 3), rep("LPS_30min", 3), rep ("LPS_1hr", 3), rep("LPS_2hrs", 3))))

metaanno$Cell_Type <- factor(metaanno$Cell_Type,levels = c("Microglia", "CECs"))

metaanno$Treatment <- factor(metaanno$Treatment, levels = c("PBS", "LPS_30min", "LPS_1hr", "LPS_2hrs"))

metaanno$Cell_Type <- as.character(metaanno$Cell_Type)
metaanno$Treatment <- as.character(metaanno$Treatment)

row.names(metaanno) <- colnames(txi$counts)

metaanno

pheatmap(go_all_four_EC1hr[c(2:10, 14:22)],          color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = 9,
         treeheight_row = 10,
         show_rownames = T,
         show_colnames = F,
         annotation = (metaanno),
         annotation_row = annorowEC1hr,
         annotation_names_row = F,
         annotation_colors = ann_colors1hr,
         legend_breaks = c(-2.5,0,2.5),
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20,
         cellwidth = 30)



#------------MG 1hr (has no 30min GO)

resp_lps_MG1hr <- dplyr::filter(normalized_counts, symbol %in% c("Nfkbia", "Tnf", "Il1b", "Cxcl1", "Ptgs2", "Cxcl10", "Cxcl2"))

#----for MG1hr only resp to lps, no nfkb signalling, pos reg cytokine prdn and pos reg TNF------****

go_all_four_MG1hr <- dplyr::filter(normalized_counts, symbol %in% c("Nfkbia", "Tnf", "Il1b", "Cxcl1", "Ptgs2", "Cxcl10", "Cxcl2"))


MG1hr_Txt <- dplyr::filter(sig_res_dds_LRT_1hr_MGonly_Txt, symbol %in% c("Nfkbia", "Tnf", "Il1b", "Cxcl1", "Ptgs2", "Cxcl10", "Cxcl2"))

write.csv(MG1hr_Txt, "data/MG1hr_Txt.csv")

go_all_four_MG1hr <- as.data.frame(go_all_four_MG1hr)

row.names(go_all_four_MG1hr) <- go_all_four_MG1hr$symbol


annorowMG1hr <- read.csv(file.path("annorowMG1hr.csv"), header = TRUE, sep = ",", row.names = 1)

annorowMG1hr<-as.data.frame(annorowMG1hr)
annorowMG1hr$GO.0032496 <- factor(annorowMG1hr$GO.0032496, exclude = "")

annorowMG1hr$GO.0032496

ann_colorsMG1hr = list(Treatment = c(PBS = "#FF8C00", LPS_1hr = "#FF00FF"), Cell_Type = c(Microglia = "#000000"))


#--go all four has only one GO i.e. resp to lps---

pheatmap(go_all_four_MG1hr[c(2:4, 8:10)],          color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = 0,
         treeheight_row = 10,
         show_rownames = T,
         show_colnames = F,
         annotation = (metaanno),
         annotation_names_row = F,
         annotation_colors = ann_colorsMG1hr,
         legend_breaks = c(-1.5,0,1.5),
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cellwidth = 25) 



#----2hrs EC est of endothelial barrier

est_endothelial_barrier_EC2hrs <- dplyr::filter(normalized_counts, symbol %in% c("Rock2", "Marveld2", "Rapgef3", "Rock1", "Tnf", "Tjp2", "Il1b", "Tnfrsf1a", "Msn", "Ikbkb", "Myd88", "Icam1", "F11r", "Rapgef1", "Cldn5", "Sox18", "Rap1b", "Afdn", "Myadm","Mapk7", "Angptl4", "Hipk1", "Gata2", "Nfe2l2", "H2-M3", "Cd40", "Tnfaip3", "Xbp1", "Rgcc", "Abl2", "Pak4", "Ccl12", "Serpine1", "Thbs1", "Ndnf", "Gper1", "Tnip2", "Kdr", "Bmpr2"))

EC_est_endothelial_barrier_EC2hrs_Txt <- dplyr::filter(sig_res_dds_LRT_2hrs_EConly_Txt, symbol %in% c("Rock2", "Marveld2", "Rapgef3", "Rock1", "Tnf", "Tjp2", "Il1b", "Tnfrsf1a", "Msn", "Ikbkb", "Myd88", "Icam1", "F11r", "Rapgef1", "Cldn5", "Sox18", "Rap1b", "Afdn", "Myadm","Mapk7", "Angptl4", "Hipk1", "Gata2", "Nfe2l2", "H2-M3", "Cd40", "Tnfaip3", "Xbp1", "Rgcc", "Abl2", "Pak4", "Ccl12", "Serpine1", "Thbs1", "Ndnf", "Gper1", "Tnip2", "Kdr", "Bmpr2"))

write.csv(EC_est_endothelial_barrier_EC2hrs_Txt, "data/EC_est_endothelial_barrier_EC2hrs_Txt.csv")

est_endothelial_barrier_EC2hrs <- as.data.frame(est_endothelial_barrier_EC2hrs)

row.names(est_endothelial_barrier_EC2hrs) <- est_endothelial_barrier_EC2hrs$symbol

ann_colorsEC2hrs = list(Treatment = c(PBS = "#FF8C00", LPS_2hrs = "#CCCC00"), Cell_Type = c(Microglia = "#000000", CECs = "#2083FB"))


annorowEC2hrs <- read.csv(file.path("annorowEC2hrs_est_endo_barrier_apoptosis.csv"), header = TRUE, sep = ",", row.names = 1)

annorowEC2hrs<-as.data.frame(annorowEC2hrs)
annorowEC2hrs$GO.0061028 <- factor(annorowEC2hrs$GO.0061028, exclude = "")
annorowEC2hrs$GO.0072577 <- factor(annorowEC2hrs$GO.0072577, exclude = "")

metaanno <- as.character(metaanno$Treatment)


pheatmap(est_endothelial_barrier_EC2hrs[c(2:4,11:13,14:16, 23:25)],          color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = 6,
         treeheight_row = 10,
         show_rownames = T,
         show_colnames = F,
         annotation = (metaanno),
         border_color = NA,
         legend_breaks = c(-2,0,2),
         annotation_colors = ann_colorsEC2hrs,
         annotation_row = annorowEC2hrs,
         annotation_names_row = F,
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cellwidth=30)


IEGs <- dplyr::filter(normalized_counts, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

EC30minIEG_Txt <- dplyr::filter(sig_res_dds_LRT_30min_EConly_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

EC1hrIEG_Txt <- dplyr::filter(sig_res_dds_LRT_1hr_EConly_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

EC2hrsIEG_Txt <- dplyr::filter(sig_res_dds_LRT_2hrs_EConly_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))


write.csv(EC30minIEG_Txt, "data/EC30minIEG_Txt.csv")
write.csv(EC1hrIEG_Txt, "data/EC1hrIEG_Txt.csv")
write.csv(EC2hrsIEG_Txt, "data/EC2hrsIEG_Txt.csv")

MG30minIEG_Txt <- dplyr::filter(sig_res_dds_LRT_30min_MGonly_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

MG1hrIEG_Txt <- dplyr::filter(sig_res_dds_LRT_1hr_MGonly_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

MG2hrsIEG_Txt <- dplyr::filter(sig_res_dds_LRT_2hrs_MGonly_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))


write.csv(MG30minIEG_Txt, "data/MG30minIEG_Txt.csv")
write.csv(MG1hrIEG_Txt, "data/MG1hrIEG_Txt.csv")
write.csv(MG2hrsIEG_Txt, "data/MG2hrsIEG_Txt.csv")


IEGs <- as.data.frame(IEGs)

row.names(IEGs) <- IEGs$symbol

ann_colorsIEGs = list(Treatment = c(PBS = "#FF8C00", LPS_30min = "#00ff00", LPS_1hr = "#FF00FF", LPS_2hrs = "#CCCC00"), Cell_Type = c(Microglia = "#000000", CECs = "#2083FB"))


pheatmap(IEGs[c(2:25)],          
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = 12,
         show_rownames = T,
         show_colnames = F,
         annotation = (metaanno),
         border_color = NA,
         legend_breaks = c(-3,0,3),
         annotation_colors = ann_colorsIEGs,
         treeheight_row = 10,
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cellwidth = 15)

library(RColorBrewer)


cell_specific <- dplyr::filter(normalized_counts, symbol %in% c("Tmem119", "Aif1", "P2ry12", "Siglech", "Csf1r", "Cdh5", "Ocln", "Slc2a1", "Cldn5", "Tjp1"))

cell_specific <- as.data.frame(cell_specific)

row.names(cell_specific) <- cell_specific$symbol

ann_colors_cellspec = list(Cell_Type = c(Microglia = "#000000", CECs = "#2083FB"))


annocelltype <- read.csv(file.path("annocelltype.csv"), header = TRUE, sep = ",", row.names = 1)

annocelltype<-as.data.frame(annocelltype)
annocelltype
annocelltype$Cell_Type_specific_genes <- factor(annocelltype$Cell_Type_specific_genes, levels = c("Enriched in Microglia", "Enriched in CECs"))

annocelltype 

meta_celltype <- metaanno

meta_celltype$Treatment <- NULL

pheatmap(cell_specific[c(2:4, 14:16)],          
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = 3,
         treeheight_row = 10,
         show_rownames = T,
         show_colnames = F,
         annotation = (meta_celltype),
         annotation_row = annocelltype,
         annotation_names_row = F,
         annotation_colors = ann_colors_cellspec,
         border_color = NA,
         legend_breaks = c(-1,0,1),
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

#------------KEGG GSEA and SPIA

#---For MG 30min

res_dds_LRT_30min_MGonly_tb_entrez <- merge(res_dds_LRT_30min_MGonly_tb, annotations[, c("gene_id", "entrezid")], by.x="gene", by.y="gene_id") %>% as_tibble()


res_dds_LRT_30min_MGonly_tb_entrez  <- filter(res_dds_LRT_30min_MGonly_tb_entrez, entrezid != "NA")

res_dds_LRT_30min_MGonly_tb_entrez <- res_dds_LRT_30min_MGonly_tb_entrez[which(duplicated(res_dds_LRT_30min_MGonly_tb_entrez$entrezid) == F), ]

foldchanges_res_dds_LRT_30min_MGonly <- res_dds_LRT_30min_MGonly_tb_entrez$log2FoldChange

names(foldchanges_res_dds_LRT_30min_MGonly) <- res_dds_LRT_30min_MGonly_tb_entrez$entrezid

foldchanges_res_dds_LRT_30min_MGonly <- sort(foldchanges_res_dds_LRT_30min_MGonly, decreasing = TRUE)


sig_30min_MGonly_entrezid <- dplyr::filter(res_dds_LRT_30min_MGonly_tb_entrez, padj < 0.05 & abs(log2FoldChange) > 0.58)

sig_30min_MGonly_entrezid <<-dplyr::filter(sig_30min_MGonly_entrezid, entrezid != "NA")

entrezid_sig_30min_MGonly_entrezid <- as.character(sig_30min_MGonly_entrezid$entrezid)

kk_30min_MGonly <- enrichKEGG(gene = entrezid_sig_30min_MGonly_entrezid, universe=names(foldchanges_res_dds_LRT_30min_MGonly), organism = 'mmu', pvalueCutoff = 0.05, qvalueCutoff = 0.01, keyType =  "ncbi-geneid")


enrichKEGG_kk_30min_MGonly <- data.frame(kk_30min_MGonly)

write.csv(enrichKEGG_kk_30min_MGonly, "data/enrichKEGG_kk_30min_MGonly.csv")

dotplot(kk_30min_MGonly, 
        showCategory = 10, 
        title = "Enriched Pathways in Microglia at LPS 30min",
        font.size = 8)


sig_foldchanges_30min_MGonly_entrezid <- sig_30min_MGonly_entrezid$log2FoldChange

names(sig_foldchanges_30min_MGonly_entrezid) <- sig_30min_MGonly_entrezid$entrezid

spia_30min_MGonly <- spia(de=sig_foldchanges_30min_MGonly_entrezid, all = names(foldchanges_res_dds_LRT_30min_MGonly), organism = "mmu", nB=3000, plots = TRUE)

view(spia_30min_MGonly)



#---For MG 1hr

res_dds_LRT_1hr_MGonly_tb_entrez <- merge(res_dds_LRT_1hr_MGonly_tb, annotations[, c("gene_id", "entrezid")], by.x="gene", by.y="gene_id") %>% as.tibble()


res_dds_LRT_1hr_MGonly_tb_entrez  <- filter(res_dds_LRT_1hr_MGonly_tb_entrez, entrezid != "NA")

res_dds_LRT_1hr_MGonly_tb_entrez <- res_dds_LRT_1hr_MGonly_tb_entrez[which(duplicated(res_dds_LRT_1hr_MGonly_tb_entrez$entrezid) == F), ]

foldchanges_res_dds_LRT_1hr_MGonly <- res_dds_LRT_1hr_MGonly_tb_entrez$log2FoldChange

names(foldchanges_res_dds_LRT_1hr_MGonly) <- res_dds_LRT_1hr_MGonly_tb_entrez$entrezid

foldchanges_res_dds_LRT_1hr_MGonly <- sort(foldchanges_res_dds_LRT_1hr_MGonly, decreasing = TRUE)


sig_1hr_MGonly_entrezid <- dplyr::filter(res_dds_LRT_1hr_MGonly_tb_entrez, padj < 0.05 & abs(log2FoldChange) > 0.58)

sig_1hr_MGonly_entrezid <<-dplyr::filter(sig_1hr_MGonly_entrezid, entrezid != "NA")

entrezid_sig_1hr_MGonly_entrezid <- as.character(sig_1hr_MGonly_entrezid$entrezid)

kk_1hr_MGonly <- enrichKEGG(gene = entrezid_sig_1hr_MGonly_entrezid, universe=names(foldchanges_res_dds_LRT_1hr_MGonly), organism = 'mmu', pvalueCutoff = 0.05, qvalueCutoff = 0.01, keyType =  "ncbi-geneid")

enrichKEGG_kk_1hr_MGonly <- data.frame(kk_1hr_MGonly)

write.csv(enrichKEGG_kk_1hr_MGonly, "data/enrichKEGG_kk_1hr_MGonly.csv")

dotplot(kk_1hr_MGonly, 
        showCategory = 10, 
        title = "Enriched Pathways in Microglia at LPS 1hr",
        font.size = 8)


sig_foldchanges_1hr_MGonly_entrezid <- sig_1hr_MGonly_entrezid$log2FoldChange

names(sig_foldchanges_1hr_MGonly_entrezid) <- sig_1hr_MGonly_entrezid$entrezid

spia_1hr_MGonly <- spia(de=sig_foldchanges_1hr_MGonly_entrezid, all = names(foldchanges_res_dds_LRT_1hr_MGonly), organism = "mmu", nB=3000, plots = TRUE)

view(spia_1hr_MGonly)



#---For MG 2hrs

res_dds_LRT_2hrs_MGonly_tb_entrez <- merge(res_dds_LRT_2hrs_MGonly_tb, annotations[, c("gene_id", "entrezid")], by.x="gene", by.y="gene_id") %>% as.tibble()


res_dds_LRT_2hrs_MGonly_tb_entrez  <- filter(res_dds_LRT_2hrs_MGonly_tb_entrez, entrezid != "NA")

res_dds_LRT_2hrs_MGonly_tb_entrez <- res_dds_LRT_2hrs_MGonly_tb_entrez[which(duplicated(res_dds_LRT_2hrs_MGonly_tb_entrez$entrezid) == F), ]

foldchanges_res_dds_LRT_2hrs_MGonly <- res_dds_LRT_2hrs_MGonly_tb_entrez$log2FoldChange

names(foldchanges_res_dds_LRT_2hrs_MGonly) <- res_dds_LRT_2hrs_MGonly_tb_entrez$entrezid

foldchanges_res_dds_LRT_2hrs_MGonly <- sort(foldchanges_res_dds_LRT_2hrs_MGonly, decreasing = TRUE)


sig_2hrs_MGonly_entrezid <- dplyr::filter(res_dds_LRT_2hrs_MGonly_tb_entrez, padj < 0.05 & abs(log2FoldChange) > 0.58)

sig_2hrs_MGonly_entrezid <<-dplyr::filter(sig_2hrs_MGonly_entrezid, entrezid != "NA")

entrezid_sig_2hrs_MGonly_entrezid <- as.character(sig_2hrs_MGonly_entrezid$entrezid)

kk_2hrs_MGonly <- enrichKEGG(gene = entrezid_sig_2hrs_MGonly_entrezid, universe=names(foldchanges_res_dds_LRT_2hrs_MGonly), organism = 'mmu', pvalueCutoff = 0.05, qvalueCutoff = 0.01, keyType =  "ncbi-geneid")

enrichKEGG_kk_2hrs_MGonly <- data.frame(kk_2hrs_MGonly)

write.csv(enrichKEGG_kk_2hrs_MGonly, "data/enrichKEGG_kk_2hrs_MGonly.csv")

dotplot(kk_2hrs_MGonly, 
        showCategory = 10, 
        title = "Enriched Pathways in Microglia at LPS 2hrs",
        font.size = 8)


sig_foldchanges_2hrs_MGonly_entrezid <- sig_2hrs_MGonly_entrezid$log2FoldChange

names(sig_foldchanges_2hrs_MGonly_entrezid) <- sig_2hrs_MGonly_entrezid$entrezid

spia_2hrs_MGonly <- spia(de=sig_foldchanges_2hrs_MGonly_entrezid, all = names(foldchanges_res_dds_LRT_2hrs_MGonly), organism = "mmu", nB=3000, plots = TRUE)

view(spia_2hrs_MGonly)

subset(spia_2hrs_MGonly, ID == "04064")


#---For EC 30min

res_dds_LRT_30min_EConly_tb_entrez <- merge(res_dds_LRT_30min_EConly_tb, annotations[, c("gene_id", "entrezid")], by.x="gene", by.y="gene_id") %>% as.tibble()


res_dds_LRT_30min_EConly_tb_entrez  <- filter(res_dds_LRT_30min_EConly_tb_entrez, entrezid != "NA")

res_dds_LRT_30min_EConly_tb_entrez <- res_dds_LRT_30min_EConly_tb_entrez[which(duplicated(res_dds_LRT_30min_EConly_tb_entrez$entrezid) == F), ]

foldchanges_res_dds_LRT_30min_EConly <- res_dds_LRT_30min_EConly_tb_entrez$log2FoldChange

names(foldchanges_res_dds_LRT_30min_EConly) <- res_dds_LRT_30min_EConly_tb_entrez$entrezid

foldchanges_res_dds_LRT_30min_EConly <- sort(foldchanges_res_dds_LRT_30min_EConly, decreasing = TRUE)


sig_30min_EConly_entrezid <- dplyr::filter(res_dds_LRT_30min_EConly_tb_entrez, padj < 0.05 & abs(log2FoldChange) > 0.58)

sig_30min_EConly_entrezid <<-dplyr::filter(sig_30min_EConly_entrezid, entrezid != "NA")

entrezid_sig_30min_EConly_entrezid <- as.character(sig_30min_EConly_entrezid$entrezid)

kk_30min_EConly <- enrichKEGG(gene = entrezid_sig_30min_EConly_entrezid, universe=names(foldchanges_res_dds_LRT_30min_EConly), organism = 'mmu', pvalueCutoff = 0.05, qvalueCutoff = 0.01, keyType =  "ncbi-geneid")

enrichKEGG_kk_30min_EConly <- data.frame(kk_30min_EConly)

write.csv(enrichKEGG_kk_30min_EConly, "data/enrichKEGG_kk_30min_EConly.csv")

dotplot(kk_30min_EConly, 
        showCategory = 10, 
        title = "Enriched Pathways in CECs at LPS 30min",
        font.size = 8)



sig_foldchanges_30min_EConly_entrezid <- sig_30min_EConly_entrezid$log2FoldChange

names(sig_foldchanges_30min_EConly_entrezid) <- sig_30min_EConly_entrezid$entrezid

spia_30min_EConly <- spia(de=sig_foldchanges_30min_EConly_entrezid, all = names(foldchanges_res_dds_LRT_30min_EConly), organism = "mmu", nB=3000, plots = TRUE)

view(spia_30min_EConly)

pathview(gene.data = foldchanges_res_dds_LRT_30min_EConly, pathway.id = "mmu04064", species = "mmu", limit = list(gene = 1.58))


#---For CEC 1hr

res_dds_LRT_1hr_EConly_tb_entrez <- merge(res_dds_LRT_1hr_EConly_tb, annotations[, c("gene_id", "entrezid")], by.x="gene", by.y="gene_id") %>% as.tibble()


res_dds_LRT_1hr_EConly_tb_entrez  <- filter(res_dds_LRT_1hr_EConly_tb_entrez, entrezid != "NA")

res_dds_LRT_1hr_EConly_tb_entrez <- res_dds_LRT_1hr_EConly_tb_entrez[which(duplicated(res_dds_LRT_1hr_EConly_tb_entrez$entrezid) == F), ]

foldchanges_res_dds_LRT_1hr_EConly <- res_dds_LRT_1hr_EConly_tb_entrez$log2FoldChange

names(foldchanges_res_dds_LRT_1hr_EConly) <- res_dds_LRT_1hr_EConly_tb_entrez$entrezid

foldchanges_res_dds_LRT_1hr_EConly <- sort(foldchanges_res_dds_LRT_1hr_EConly, decreasing = TRUE)


sig_1hr_EConly_entrezid <- dplyr::filter(res_dds_LRT_1hr_EConly_tb_entrez, padj < 0.05 & abs(log2FoldChange) > 0.58)

sig_1hr_EConly_entrezid <<-dplyr::filter(sig_1hr_EConly_entrezid, entrezid != "NA")

entrezid_sig_1hr_EConly_entrezid <- as.character(sig_1hr_EConly_entrezid$entrezid)

kk_1hr_EConly <- enrichKEGG(gene = entrezid_sig_1hr_EConly_entrezid, universe=names(foldchanges_res_dds_LRT_1hr_EConly), organism = 'mmu', pvalueCutoff = 0.05, qvalueCutoff = 0.01, keyType =  "ncbi-geneid")

enrichKEGG_kk_1hr_EConly <- data.frame(kk_1hr_EConly)

write.csv(enrichKEGG_kk_1hr_EConly, "data/enrichKEGG_kk_1hr_EConly.csv")

dotplot(kk_1hr_EConly, 
        showCategory = 10, 
        title = "Enriched Pathways in CECs at LPS 1hr",
        font.size = 8)

sig_foldchanges_1hr_EConly_entrezid <- sig_1hr_EConly_entrezid$log2FoldChange

names(sig_foldchanges_1hr_EConly_entrezid) <- sig_1hr_EConly_entrezid$entrezid

spia_1hr_EConly <- spia(de=sig_foldchanges_1hr_EConly_entrezid, all = names(foldchanges_res_dds_LRT_1hr_EConly), organism = "mmu", nB=3000, plots = TRUE)

view(spia_1hr_EConly)



#---For MG 2hrs

res_dds_LRT_2hrs_EConly_tb_entrez <- merge(res_dds_LRT_2hrs_EConly_tb, annotations[, c("gene_id", "entrezid")], by.x="gene", by.y="gene_id") %>% as.tibble()


res_dds_LRT_2hrs_EConly_tb_entrez  <- filter(res_dds_LRT_2hrs_EConly_tb_entrez, entrezid != "NA")

res_dds_LRT_2hrs_EConly_tb_entrez <- res_dds_LRT_2hrs_EConly_tb_entrez[which(duplicated(res_dds_LRT_2hrs_EConly_tb_entrez$entrezid) == F), ]

foldchanges_res_dds_LRT_2hrs_EConly <- res_dds_LRT_2hrs_EConly_tb_entrez$log2FoldChange

names(foldchanges_res_dds_LRT_2hrs_EConly) <- res_dds_LRT_2hrs_EConly_tb_entrez$entrezid

foldchanges_res_dds_LRT_2hrs_EConly <- sort(foldchanges_res_dds_LRT_2hrs_EConly, decreasing = TRUE)


sig_2hrs_EConly_entrezid <- dplyr::filter(res_dds_LRT_2hrs_EConly_tb_entrez, padj < 0.05 & abs(log2FoldChange) > 0.58)

sig_2hrs_EConly_entrezid <<-dplyr::filter(sig_2hrs_EConly_entrezid, entrezid != "NA")

entrezid_sig_2hrs_EConly_entrezid <- as.character(sig_2hrs_EConly_entrezid$entrezid)

kk_2hrs_EConly <- enrichKEGG(gene = entrezid_sig_2hrs_EConly_entrezid, universe=names(foldchanges_res_dds_LRT_2hrs_EConly), organism = 'mmu', pvalueCutoff = 0.05, qvalueCutoff = 0.01, keyType =  "ncbi-geneid")

enrichKEGG_kk_2hrs_EConly <- data.frame(kk_2hrs_EConly)

write.csv(enrichKEGG_kk_2hrs_EConly, "data/enrichKEGG_kk_2hrs_EConly.csv")

dotplot(kk_2hrs_EConly, 
        showCategory = 10, 
        title = "Enriched Pathways in CECs at LPS 2hrs",
        font.size = 8)


sig_foldchanges_2hrs_EConly_entrezid <- sig_2hrs_EConly_entrezid$log2FoldChange

names(sig_foldchanges_2hrs_EConly_entrezid) <- sig_2hrs_EConly_entrezid$entrezid

spia_2hrs_EConly <- spia(de=sig_foldchanges_2hrs_EConly_entrezid, all = names(foldchanges_res_dds_LRT_2hrs_EConly), organism = "mmu", nB=3000, plots = TRUE)

view(spia_2hrs_EConly)


