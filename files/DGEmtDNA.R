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
samples <- list.files(path = "~/rnaseq/06-DGEmtDNA", full.names = T, pattern="salmon$")

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Give unique name for each element
names(files) <- str_replace(samples, "/Users/tatiana/rnaseq/06-DGEmtDNA/", "") %>% 
  str_replace("_mtDNA-salmon", "") 

## Load the annotation table for the mtDNA genes
tx2gene <- read.delim("~/rnaseq/00-Databases/mtDNA_annot.txt")

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
       
       