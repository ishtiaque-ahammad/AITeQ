#Read the CSV file
data <- read.csv("all_merged_new.csv", header = TRUE, sep = ",", row.names = 1)
# Select only the columns with gene expression values (excluding the first column "Geneid")
#expression_data <- data[, -1]
#pdata <- expression_data+1
# Apply log transformation
#log_transformed_data <- log(expdata)
# Add the gene IDs back to the log-transformed data
#log_transformed_data <- cbind(data$Geneid, log_transformed_data)
# Convert the log-transformed data back to a data frame
#log_transformed_data <- as.data.frame(log_transformed_data)
# Save the log-transformed data to a new CSV file
#write.csv(log_transformed_data, "log_transformed_counts.csv", row.names = FALSE)

#batch correction using ComBatseq
library(nlme)
library(mgcv)
library(genefilter)
library(BiocParallel)
library(sva)

allsampleinfo1 <- read.csv("ML_sampleinfo.csv", header = TRUE, row.names = 1)
#allmerged  <- read.csv("log_transformed_counts.csv", row.names = 1, header = TRUE)
count_matrix <- as.matrix(data)
batch <- allsampleinfo1$batch
adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)

#DESeq2
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(usethis)
library(devtools)
#install_github('jokergoo/ComplexHeatmap')
library(ComplexHeatmap)
#making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix) %in% rownames(allsampleinfo1))
#are they in the same order?
all(colnames(count_matrix) == rownames(allsampleinfo1))

#Step 2: construct a DESeqDataSet object --------
#count_matrix <- round(count_matrix)
dds <- DESeqDataSetFromMatrix(countData=count_matrix,
                              colData = allsampleinfo1,
                              design = ~batch+condition)
dds
#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total
filtered <- rowSums(counts(dds)) >= 50

#if the rows fulfill the condition then the value will be True, if not it will be False
#To use the logical (T/F) values to subset the DESeq2 data set object
dds_new <- dds[filtered,]
dds_new

#Set reference factor for doing differential gene expression analysis
dds_new$condition <- as.factor(dds_new$condition) #converting in a factor
factorlevel <- relevel(dds_new$condition, ref = "control")
factorlevel
#Now Control is the reference level. 
#If we don't explicitly mention which level we want to use as a reference level, it will automatically set the level alphabatically.

#NOTE: Remove collapse technical replicates, if present 
#Step 3: Run DESeq2
dds_new_deseq <- DESeq(dds_new)
res <- results(dds_new_deseq)
res
#Here the log2fold change is calculated for in the design factor Condition between Alzheimers Disease and Control.
#So, whatever values we see here for the fold changes are all in the Alzheimers disease because these are compared with the control
#The statistical test is Wald test 
#basemean is the average of the normalized counts taken over all the samples
#The log2foldchange is the fold change of the genes in the disease condition compared with the control #positive value = upregulated #negative value = downregulated
#lfcSE = standard error that estimates for the log2fold change
#stat = Wald test values of the genes
#pvalue = test statistic for these genes
#padj = the corrected pvalues for multiple testing #NOTE: we need to correct pvalues for multiple testing 
#because whenever we perform a statistical test, we use a pvalue of 0.05, so 5% of our DEGs are not really differentially expressed but they are only due to random chance and the genes have no real effect on patients
#that means in our dataset we have around 42 thousand genes, so, 5% of the 42000 genes is 2100 genes. so in our list, 2100 of those genes are false positive. that is why padj values method is being used to avoid false positive genes
head(results(dds_new_deseq, tidy=TRUE))
#Explore results
summary(res)
res0.05 <- results(dds_new_deseq, alpha = 0.05)
summary(res0.05)
#contrasts
resultsNames(dds_new_deseq)
#MA plot
plotMA(res, ylim =c(-5,5))
#PCA plot
vsdB <- varianceStabilizingTransformation(dds_new_deseq)
plotPCA(vsdB, intgroup=c("condition", "batch"))

#variance
vsdB_table <- as.data.frame(assay(vsdB))
write.csv(vsdB_table, "vsdB_table.csv")
vsdB_table_rowsum <- transform(vsdB_table, sum = rowSums(vsdB_table))
colnames(vsdB_table_rowsum)
selected <- order(vsdB_table_rowsum$sum, decreasing = TRUE) [1:100]
vsdB_table[selected,]
#Volcano plot
with(res0.05, plot(log2FoldChange, -log10(pvalue), pch=20, main = "Volcano plot", xlim = c(-5,5), ylim = c(0,8)))

with(res0.05, plot(log2FoldChange, -log10(pvalue), pch=20, main = "Differentally Expressed Genes (DEGs)", xlim = c(-5,5), ylim = c(0,8),
                   col = ifelse(pvalue < 0.05, "orange", "red")))




#gene naming
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05,]
sigs
write.csv(sigs, file = "deseq_results.csv")

sigs.df <- as.data.frame(sigs)
#sigs.df <- sigs.df[(abs(sigs.df$log2FoldChange) > 1),]
sigs.df <- sigs.df[(sigs.df$baseMean > 100) & (abs(sigs.df$log2FoldChange) > 1.45),]
sigs.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", 
                         column = "SYMBOL")
sigs.df
naomit_sigsdf <- na.omit(sigs.df)
naomit_sigsdf
write.csv(naomit_sigsdf, file = "ControlvsAD_gene_name.csv")
#sigs.df1 <- read.csv("ControlvsAD_gene_name.csv", header = TRUE, row.names = 1)
# Apply the filtering condition
#filtered_df <- sigs.df[(sigs.df$baseMean > 3.47) & (abs(sigs.df$log2FoldChange) > 0.095), ]
#naomit_filtered_df <- na.omit(filtered_df)
#write.csv(naomit_filtered_df, file = "ControlvsAD_gene_name.csv")

#heatmap
mat <- counts(dds_new_deseq, normalized = T)[rownames(naomit_sigsdf),]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(colData)
mat.z
h <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels =colnames(mat.z),
             name = "Z-score", row_labels = sigs.df[rownames(mat.z),]$symbol)
png('simple_heatmap.png', res = 250, width = 1000, height = 2000)
print(h)
dev.off()