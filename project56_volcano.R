#installing and loading packages
install.packages('pacman')
pacman::p_load(pacman, dplyr, tidyr,readxl, ggplot2, ggrepel, readr, cowplot)

ggrepel.max.overlaps = Inf
#Importing data 

mydata <- read_xlsx("/Users/hinagaur/Downloads/res1.xlsx")
# mydata <- mydata %>% mutate_all(~ifelse(is.na(.), 0, .))


#Displaying column names

colnames(mydata)

#creating a new column stating which genes are UP and DOWN regulated
mydata$diffexpressed <- "NO"
mydata$diffexpressed[mydata$log2FoldChange > 1 & mydata$padj < 0.05] <- "UP"
mydata$diffexpressed[mydata$log2FoldChange <  -1 & mydata$padj< 0.05] <- "DOWN"

mydata$delabel <- NA

head(mydata)
View(mydata)

#Setting colours for UP, DOWN, AND INSIGNIFICANT GENES
mycolors <- c('blue', 'red', 'gray')
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(data = mydata, aes(x= log2FoldChange, y = -log(pvalue), col= diffexpressed , label = delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  labs(color = "Diffexpressed",
       x= expression("log"*"FC"),y= expression("-log"[10]*"(p-value)"))+
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

#Filtering of up and downregulated genes
up_reg <- subset(mydata, diffexpressed == "UP")
down_reg <- subset(mydata, diffexpressed == 'DOWN')

#Annotating the gene IDs with the transcript IDs

gene_meta <- read.csv("/Users/hinagaur/Downloads/tx2gene.RefSeq.All.tsv", header = TRUE)

ann_df <- merge(gene_meta, mydata, by.x = 'TXNAME' , by.y = 'ID')

##New volcano plot with the top 10 up and down regulated genes 

#creating a new column stating which genes are UP and DOWN regulated
ann_df$diffexpressed <- "NO"
ann_df$diffexpressed[ann_df$log2FoldChange > 1 & ann_df$padj < 0.05] <- "UP"
ann_df$diffexpressed[ann_df$log2FoldChange <  -1 & ann_df$padj< 0.05] <- "DOWN"

ann_df$delabel <- ifelse(ann_df$TXNAME %in% head(ann_df[order(ann_df$padj), "TXNAME"], 10), ann_df$TXNAME, NA)


#Setting colours for UP, DOWN, AND INSIGNIFICANT GENES
mycolors <- c('blue', 'red', 'gray')
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(data = ann_df, aes(x= log2FoldChange, y = -log(pvalue), col= diffexpressed , label = delabel))+
  geom_point()+
  theme_minimal()+
  labs(color = "Diffexpressed",
       x= expression("log"[2]*"FC"),y= expression("-log"[10]*"(p-value)"))+
  scale_color_manual(values = mycolors)+
  geom_text_repel(aes(colour= NULL), max.overlaps = Inf)+
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

df <- subset(ann_df, diffexpressed == 'UP' | diffexpressed == 'DOWN')
new_df <- arrange(df, padj)

#Printing the first three rows, displaying the top two most significantly expressed genes
print(new_df[1:2,])

#loading count matrix
count <- read.csv("/Users/hinagaur/Downloads/count.csv", header = TRUE, check.names = TRUE) 
row.names(count) <- count$X #specifying 
count[1] <- NULL

count_t <- t(count)
# rownames(count_t) <- colnames(count)
count_t <- tibble::rownames_to_column(as.data.frame(count_t), "Sample") # converting rownames to first column

count_t$Sample <- sub('\\..*', '', count_t$Sample)

# write.csv(count_t, "~/desktop/count_t.csv")

#Box plot
gg_p1 <- ggplot(count_t, aes(x = Sample, y = NM_001854.4, fill = Sample)) +
  geom_boxplot() +
  labs(title = "Boxplot for NM_001854.4", y = "Counts across all the samples") +
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
  labs(title = "Boxplot for XM_047440587.1", y = "Counts across all the samples") +
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



filt_df <- subset(ann_df, !is.na(delabel))

rownames(count_t)
