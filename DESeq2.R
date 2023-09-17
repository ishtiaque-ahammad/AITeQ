#Read the CSV file
data <- read.csv("merged_counts.csv", header = TRUE, sep = ",", row.names = 1)

#batch correction using ComBatseq
library(nlme)
library(mgcv)
library(genefilter)
library(BiocParallel)
library(sva)

#Please provide sample metadata in a file called sampleinfo.csv. 
#The following columns are needed: sample_id, condition, batch (in case of meta-analysis involving multiple projects). 
sampleinfo <- read.csv("sampleinfo.csv", header = TRUE, row.names = 1)

#For batch correction, run the following 3 lines of code.
count_matrix <- as.matrix(data)
batch <- sampleinfo$batch
adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)

#DESeq2
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(usethis)
library(devtools)

#Make sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix) %in% rownames(sampleinfo))

#Make sure they in the same order
all(colnames(count_matrix) == rownames(sampleinfo))

#Step 2: construct a DESeqDataSet object --------
dds <- DESeqDataSetFromMatrix(countData=count_matrix,
                              colData = sampleinfo,
                              design = ~batch+condition)
dds

#Filter out rows with low gene counts
filtered <- rowSums(counts(dds)) >= 50
dds_new <- dds[filtered,]
dds_new

#Set the reference factor for differential gene expression analysis
dds_new$condition <- as.factor(dds_new$condition) #converting each condition to a factor
factorlevel <- relevel(dds_new$condition, ref = "control")
factorlevel

#Now Control is the reference level. 
#If we don't explicitly mention which level we want to use as a reference level, it will automatically set the level alphabatically.

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

#Contrasts
resultsNames(dds_new_deseq)

#Variance Stabilizing Transformation
VST_tranform <- varianceStabilizingTransformation(dds_new_deseq)
VST <- as.data.frame(assay(VST_tranform))
write.csv(VST, "VST.csv")