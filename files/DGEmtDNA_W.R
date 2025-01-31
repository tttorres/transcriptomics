# Analysis of differential gene expression using DESeq2----

# We will use several R, that need to be installed once, and loaded every time we
# need them.

# First time:
# setRepositories() 
# install.packages("DESeq2")
# install.packages("tidyverse")
# install.packages("RColorBrewer")
# install.packages("pheatmap")
# install.packages("DEGreport")
# install.packages("tximport")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("beeper")

# Setup----
## Loading Bioconductor and CRAN libraries
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport) 
library(ggplot2)
library(ggrepel)
library(beepr)

# Step 1: Load data ----
## List all directories containing data  
samples <- list.files(path = "/Users/Informatica/Downloads/06-DGEmtDNA", full.names = T, pattern="salmon$")

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Give unique name for each element
names(files) <- str_replace(samples, "/Users/Informatica/Downloads/06-DGEmtDNA/", "") %>% 
  str_replace("_mtDNA-salmon", "") 

## Load the annotation table for the mtDNA genes
tx2gene <- read.delim("/Users/Informatica/Downloads/06-DGEmtDNA/mtDNA_annot.txt")

## Check the annotation table 
tx2gene %>% View()

## Run tximport
txi <- tximport(files, type="salmon", 
                tx2gene=tx2gene[,c("transcript", "gene")], 
                countsFromAbundance="lengthScaledTPM")
beep(sound = 2)

## Check the object atributes
class(txi)
names(txi)
attributes(txi)

## View at the counts
txi$counts %>% View()

## Create an object with the counts
data <- txi$counts %>% 
  round() %>% 
  data.frame()

## Create metadata
sampleNames <- colnames(data) # Extract sample names

sampletype <- substr(sampleNames, 1, 5)
species <- substr(sampleNames, 1, 4) # Species = First four characters
stage <- substr(sampleNames, 5, 5)   # Stage = Fifth characters (F or L)

stage <- ifelse(stage == "F", "female", "larva") # Replacing F with "female"              
                                                 # and L with "larva"

meta <- data.frame(sampletype, species, stage, row.names = sampleNames)
meta #show dataframe


print("Step 1 complete!")

# Step 2: Data exploration----
## Distribution of RNA-Seq counts
ggplot(data) +
  geom_histogram(aes(x = CalbF_1), stat = "bin", bins = 10) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

## Assessing mean and variance of the data
mean_counts <- apply(data[,1:32], 1, mean) #vector of means       
variance_counts <- apply(data[,1:32], 1, var) #a vector of variances
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e12)) +
  scale_x_log10(limits = c(1,1e12)) +
  geom_abline(intercept = 0, slope = 1, color="red")

# Step 3: Count normalization----
## Check if metadata and counts (txi) match
all(colnames(txi$counts) %in% rownames(meta)) #checks if all the column names of 
                                              #the matrix txi$counts are present 
                                              #in the rows of the object meta
all(colnames(txi$counts) == rownames(meta)) # checks if the column names of txi$counts 
                                            # are identical to the row names of the 
                                            # object meta and in the exact same order.

## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
View(counts(dds))

## Normalize counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalizationFactors(dds) # View normalization factors 

# Step 4: Exploratory analysis----
## Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
beep(sound = 2)

## Plot PCA
plotPCA(rld, intgroup="sampletype")
plotPCA(rld, intgroup="species")
plotPCA(rld, intgroup="stage")

## Extract the log matrix from the object
rld_mat <- assay(rld)    

## Compute pairwise correlation values
rld_cor <- cor(rld_mat)    #cor() is a base R function
rld_cor #checking the output of cor()

## Plot the heatmap
pheatmap(rld_cor, annotation = meta)

# Step 5: Testing Gene Expression Differences----
## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
beep(sound = 2)

## Run DESeq analysis
dds <- DESeq(dds)
beep(sound = 2)

# Step 6: Exploring results----
## Define contrasts 
## contrast <- c("condition", "level_to_compare", "base_level")
contrast_oe <- c("sampletype", "ChomF", "CmacF")

## Extract results for ChomF vs CmacF
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.05)

## Check what type of object is returned
class(res_tableOE)

## View information stored in results
res_tableOE %>% data.frame() %>% View()

## Plot log2 fold changes (on the y-axis) versus the mean of normalized counts 
## (on the x-axis).
plotMA(res_tableOE, ylim=c(-3,3))
       
## Summarize results
summary(res_tableOE, alpha = 0.05)

## Extracting significant differentially expressed genes
## Set thresholds
padj.cutoff <- 0.05

## Create a tibble (an enhanced version of a data frame) of results
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_tableOE_tb

## Filter the tibble to keep only significant genes
sigOE <- res_tableOE_tb %>%
  dplyr::filter(padj < padj.cutoff)
sigOE

## Plot the expression of a single gene       
plotCounts(dds, gene="NAD5 ", intgroup="sampletype") 
plotCounts(dds, gene="NAD6 ", intgroup="sampletype") 
plotCounts(dds, gene="ATP8 ", intgroup="sampletype") 

## Using ggplot2 to plot expression of a single gene
## Save plotcounts to a data frame object
d <- plotCounts(dds, gene="NAD5 ", intgroup="sampletype", returnData=TRUE)
d %>% View() # View the output of plotCounts()

## Plot
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  theme_bw() +
  ggtitle("NAD5 ") +
  theme(plot.title = element_text(hjust = 0.5))

## Volcano Plot
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold 
## change > 1.5 

res_tableOE_tb <- res_tableOE_tb %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("mtDNA expression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))         