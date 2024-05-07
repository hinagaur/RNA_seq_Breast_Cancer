#Importing the required packages
library("tximportData")
library("tximport")
library("readr")
library("DESeq2")
library("dplyr")
library("GenomicFeatures")
library("ggplot2")
library("apeglm")
library("writexl")
library("tibble")
library(pheatmap)
library(ggrepel)
library(cowplot)
library(heatmap)
#set the maximum overlap between the points and labels
ggrepel.max.overlaps = 10000

#Creating specific reference files using RefSeq gff3 files and GenomicFeatures package
txdb <- makeTxDbFromGFF(file="/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/refTranscriptome/GRCh38_latest_genomic.gff.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene = tx2gene %>% filter(!is.na(TXNAME))

#Creating a table in .tsv format for the reference transcriptome to map transcript ID to Gene ID
write.table(tx2gene, "/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/refTranscriptome/tx2gene.RefSeq.All.tsv", 
            sep = ",", 
            row.names = FALSE)

#Reading the TSV file created
tx2gene <- read_csv("/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/refTranscriptome/tx2gene.RefSeq.All.tsv") %>% as.data.frame()
head(tx2gene)

#set the current working directory
setwd("~/biol6150/ProjectSubmissions/Group18-AloHA/Project56/transcriptQuant/")

# Read the samples map file containing sample ID and corresponding quant.sf path and condition
sample_info <- read.table("/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/SampleMetadata.tsv", header=TRUE) 
#get Sample IDs from sample_info
samples = sample_info %>% pull(SampleID)
print(samples)

#get quant.sf file paths from sample_info
file_paths = sample_info %>% pull(SampleQuantPath)
print(file_paths)
#check if the paths exist
all(file.exists(file_paths))

# Create a TXI object. Read quant.sf files created by Salmon
txi.salmon <- tximport(file_paths, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
print(head(txi.salmon$counts))

# creating map file in format as expected by DESeq2.
sampleTable <- data.frame(condition = factor(c(rep("breast_tumor",30), rep("normal_breast_tissue",30))
))

#get the column headers of the TXI object
rownames(sampleTable) <- colnames(txi.salmon$counts)
sampleTable

#Create the dds object to prepare data to run diffrential expression analysis
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)

#Run the DESeq pipeline
dds <- DESeq(dds)

#Get differential expression results
res <- results(dds)

#export res object to CSV
write.csv(as.data.frame(res), file = "/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/res.csv")

#filter res for p-adjusted value<0.05 and export to CSV
res0.05 <- res %>% subset(padj<0.05)
head(res0.05)
write.csv(as.data.frame(res0.05), file = "/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/res0.05.csv")

#filter res for p-adjusted value<0.01 and export to CSV
res0.01 <- res %>% subset(padj<0.01)
head(res0.01)
write.csv(as.data.frame(res0.01), file = "/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/res0.01.csv")

#create df from res object which we exported as CSV
res_csv <- read_csv("/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/res.csv")

#Renaming the first column with transcript IDs as "ID"
colnames(res_csv)[1] = "ID"

#Annotating the gene IDs with the transcript IDs
gene_meta <- read.csv("/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/refTranscriptome/tx2gene.RefSeq.All.tsv", header = TRUE)

#merging our data with metadata containing gene names for corresponding transcript IDs
ann_df <- merge(gene_meta, res_csv, by.x = 'TXNAME' , by.y = 'ID')

#Visualizing the results

#creating a new column stating which genes are UP and DOWN regulated
ann_df$diffexpressed <- "NO"
ann_df$diffexpressed[ann_df$log2FoldChange > 1 & ann_df$padj < 0.05] <- "UP"
ann_df$diffexpressed[ann_df$log2FoldChange <  -1 & ann_df$padj< 0.05] <- "DOWN"

#create delabel column having top10 significantly expressed gene names
ann_df$delabel <- ifelse(ann_df$TXNAME %in% head(ann_df[order(ann_df$padj), "TXNAME"], 10), ann_df$TXNAME, NA)

# Setting colours for UP, DOWN, AND INSIGNIFICANT GENES
mycolors <- c('blue', 'red', 'gray')
names(mycolors) <- c("DOWN", "UP", "NO")

#volcano plot
ggplot(data = ann_df, aes(x= log2FoldChange, y = -log(pvalue), col= diffexpressed , label = delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel(aes(label = ann_df$delabel, col = 'black'), max.overlaps = Inf)+
  labs(color = "Diffexpressed",
       x= expression("log"[2]*"FC"),y= expression("-log"[10]*"(p-value)"))+
  scale_color_manual(values = mycolors)+
  theme(
    axis.text.x = element_text(size = 12),  # Adjust size of x-axis labels
    axis.text.y = element_text(size = 12),  # Adjust size of y-axis labels
    axis.title.x = element_text(size = 14),  # Adjust size of x-axis title
    axis.title.y = element_text(size = 14),  # Adjust size of y-axis title
    legend.text = element_text(size = 12),  # Adjust size of legend text
    legend.title = element_text(size = 14)  # Adjust size of legend text
  )+
  scale_x_continuous(breaks = c(seq(-5, 5, 2.5)), # Modify x-axis tick intervals  
                     limits = c(-5, 5))+
  scale_y_continuous(limits = c(0, 60))


#Getting the names of top two most significantly expressed genes 
filt_df <- subset(ann_df, diffexpressed == 'UP' | diffexpressed == 'DOWN') %>% arrange(padj)

#Printing the top 2 and 10 most significantly expressed up- or down- regulated genes
print(filt_df[1:2,])
print(filt_df[1:10,])

#creating the count matrix from TXI object
count <- txi.salmon$abundance
colnames(count) <- sampleTable$condition

#export count data to CSV
write.csv(as.data.frame(count), file = "/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/count.csv")

#loading count matrix from CSV
count_csv <- read.csv("/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/count.csv", header = TRUE, check.names = TRUE)

# making transcript IDs as rownames in count df
row.names(count_csv) <- count_csv$X
count_csv[1] <- NULL

# Transposing the count matrix
count_t <- t(count_csv)

#converting 
count_t <- tibble::rownames_to_column(as.data.frame(count_t), "Sample") # converting rownames to first column

count_t$Sample <- sub('\\..*', '', count_t$Sample)

#Box plot
gg_p1 <- ggplot(count_t, aes(x = Sample, y = NM_001854.4, fill = Sample)) +
  geom_boxplot() +
  labs(title = "NM_001854.4 (COL11A1)", y = "Counts across all the samples") +
  scale_y_continuous(breaks = c(seq(0, 1, .2)), # Modify x-axis tick intervals  
                     limits = c(0, 1))+
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold", family = "Helvetica"),
    legend.title = element_text(size = 15),
    title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )
gg_p1

#box plot
gg_p2 <- ggplot(count_t, aes(x = Sample, y = XM_047440587.1, fill = Sample)) +
  geom_boxplot() +
  labs(title = "XM_047440587.1 (B4GALT5)", y = "Counts across all the samples") +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold", family = "Helvetica"),
    legend.title = element_text(size = 15),
    title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )
gg_p2


#printing combined plot

p_all <- cowplot::plot_grid(gg_p1, gg_p2, nrow = 1)
p_all

##heatmap 
#Omiting NAs to make a dataset of significantly expressed transcripts
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05,]

#Writing the data to a CSV File
write.csv(sigs, "/storage/ice-shared/biol6150/ProjectSubmissions/Group18-AloHA/Project56/sigs.csv")

#Converting sigs.df to a dataframe 
sigs.df <- as.data.frame(sigs)

#Applying filters to visualise the heatmap
sigs.df <- sigs.df[(sigs.df$baseMean > 150) & (abs(sigs.df$log2FoldChange) > 2),]

#extracting counts of significant transcripts
mat <- counts(dds, normalized = TRUE)[rownames(sigs.df),]

#getting the z scores
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- sampleTable$condition
  
#Creating the heatmap
pheatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z),
        name = "Z-score")



##PCA

#Creating a PCA which shows the samples in the 2D plane spanned by their first two principal components. 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup=c("condition"))

#checking the numbers of up and downregulated genes

#Filtering of up and downregulated genes
up_reg <- subset(ann_df, diffexpressed == "UP")
down_reg <- subset(ann_df, diffexpressed == 'DOWN')

