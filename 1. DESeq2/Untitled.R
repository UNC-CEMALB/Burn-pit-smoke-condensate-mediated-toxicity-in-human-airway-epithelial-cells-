setwd('/Users/alexis/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/CEMALB_DataAnalysisPM/Projects/P1011. Emission Mixtures/P1011.3. Analyses/P1011.3.1. DESeq2/')
Output = ('/Users/alexis/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/CEMALB_DataAnalysisPM/Projects/P1011. Emission Mixtures/P1011.3. Analyses/P1011.3.1. DESeq2/Output')
cur_date = "090423"

library(readxl)
library(openxlsx)
library(tidyverse)
library(reshape2)
library(data.table)
library(factoextra)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(RUVSeq)
library(DESeq2)

# reading in file
cytokine_df = data.frame(read_excel("Input/Cytokine_Data_050423.xlsx", sheet = 2)) %>%
  dplyr::rename(Dose = Condensate_Conc) %>%
  # only interested in doses 1, 25, and control
  filter(Dose %in% c(1,25,NA))

# first creating a `coldata` object that contains all the metadata for each sample
coldata = unique(cytokine_df[,1:5]) %>%
  filter(Dose %in% c(1,NA)) %>%
  select(-Dose) %>%
  # creating sample ids
  unite("SampleID", c(colnames(cytokine_df))[c(1,3,4)], sep = "_", remove = FALSE)

head(coldata)


# making a `countdata` obj that contains cytokines as rows and sample names as cols
countdata = cytokine_df %>%
  filter(Dose %in% c(1,NA)) %>%
  select(-Dose) %>%
  # creating sample ids
  unite("SampleID", c(colnames(cytokine_df)[c(1,3:4)]), sep = "_") %>%
  select(-c("Cytokine_Conc_pslog2", "Subject_ID")) %>%
  pivot_wider(names_from = "SampleID", values_from = "Cytokine_Conc") %>%
  column_to_rownames(var = "Cytokine")

head(countdata)

# Set the rownames of coldata and column names of countdata to be in the same order 
countdata <- setcolorder(countdata, as.character(coldata$SampleID))

# replacing the sample ids in the countdata file with the ids
#colnames(countdata) <- coldata$ID

head(countdata)


# First count the total number of samples and save it as a value in the global environment
nsamp <- ncol(countdata)

# Then, calculate the median expression level across all genes and all samples and save it as a value
total_median <- median(as.matrix(countdata))

# We need to temporarily add back in the cytokine column to the countdata 
# so we can filter for genes that pass the background filter
countdata <- countdata %>% 
  rownames_to_column("Cytokine")


# filtering for cytokines that have an expression greater than the total median in at least 20% of the samples
cytokines_above_background <- countdata %>% 
  pivot_longer(cols =! Cytokine, names_to = "sampleID", values_to = "expression") %>% 
  # indicates whether the expression of a cytokine for the corresponding exposure condition is above (1) or not 
  # above (0) the median of all count data
  mutate(above_median = ifelse(expression > total_median,1,0)) %>% 
  group_by(Cytokine) %>% 
  # For each cytokine, count the number of exposure conditions where the expression was greater than the median 
  # of all count data
  summarize(total_above_median = sum(above_median)) %>% 
  # Filter for cytokines that have expression above the median in at least 20% of the samples
  filter(total_above_median >= 0.2*nsamp) %>% 
  select(Cytokine) 

# Then filter the original 'countdata' dataframe for only the cytokines above background. 
countdata <- left_join(cytokines_above_background, countdata, by = "Cytokine")

dim(countdata)

countdata_T <- countdata %>% 
  pivot_longer(cols =! Cytokine, names_to = "sampleID",values_to = "expression") %>% 
  pivot_wider(names_from = Cytokine, values_from = expression)

# Then add in a column to the transposed countdata dataframe that sums expression across all cytokines for each 
# exposure condition
countdata_T$rowsum <- rowSums(countdata_T[2:ncol(countdata_T)])

# Remove samples that have no expression
countdata_T <- countdata_T %>% 
  filter(rowsum != 0)

# Take the count data filtered for correct samples, remove the 'rowsums' column
countdata_T <- countdata_T %>% 
  select(!rowsum) 

# Then, transpose it back to the correct format for analysis
countdata <- countdata_T %>%
  pivot_longer(cols =! sampleID, names_to = "Cytokine", values_to = "expression") %>% 
  pivot_wider(names_from = sampleID, values_from = "expression") %>%
  column_to_rownames(var = "Cytokine")

dim(countdata)

# IMPUTATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IF NECESSARY

# Calculate principal components using transposed count data
pca <- prcomp(t(countdata))  # WHY WAS THIS DATA NOT SCALED AND CENTERED???

options(repr.plot.width = 10, repr.plot.height = 5) #changing size
fviz_eig(pca, addlabels = TRUE)


# Make dataframe for PCA plot generation using first two components and the sample name
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID = colnames(countdata))

# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)

# Organize dataframe so we can color our points by burn condition
pca_df <- pca_df %>% 
  separate(SampleID, into = c("SubjectNo", "Condensate", "Burn_Condition"), sep = "_", remove = FALSE)

head(pca_df)

options(repr.plot.width=12, repr.plot.height=7) #changing size

# color by burn condition and condensate
ggplot(pca_df, aes(PC1, PC2, color = Burn_Condition, shape = Condensate)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = SampleID), size = 4) +
  labs(x = paste0("PC1 (",pca_percent[1],"%)"), y = paste0("PC2 (",pca_percent[2],"%)"))

countdata_for_clustering <- t(countdata)

pheatmap(scale(countdata_for_clustering), 
         cluster_rows = TRUE, cluster_cols = FALSE, fontsize_col = 7, treeheight_row = 60, show_colnames = FALSE)


# First we store the treatment IDs and exposure conditions as a separate vector
ID <- coldata$SampleID

# And differentiate our treatments and control conditions, first by grabbing groups associated with each sample
groups <- as.factor(coldata$Burn_Condition)
groups

# Then setting a control label
ctrl <- "PBS"

# And extracting a vector of just our treatment groups
trt_groups <- setdiff(groups, ctrl)

# Let's view this vector
trt_groups

exprSet <- newSeqExpressionSet(as.matrix(countdata), phenoData = 
                                 data.frame(groups, row.names = colnames(countdata)))
colors <- brewer.pal(4, "Set2")
plotRLE(exprSet, outline = FALSE, ylim = c(-4, 4), col = colors[groups])

plotPCA(exprSet, col = colors[groups], cex = 1.2)

# Construct a matrix specifying the replicates (samples of the same exposure condition) 
# for running RUV
differences <- makeGroups(groups)

# Viewing this new matrix
head(differences)

# Now capture 1 factor (k=1) of unwanted variation
ruv_set <- RUVs(exprSet, rownames(countdata), k = 1, differences) 

pData(ruv_set)

# Viewing the head of the normalized count data, accounting for unwanted variation
head(normCounts(ruv_set))

plotRLE(ruv_set, outline=FALSE, ylim=c(-4, 4), col=colors[groups]) # GAVE ME A WARNING CAUSE IT KNOWS IT'S NOT COUNT DATA

# Set up our experiment using our RUV adjusted count and phenotype data.
# Our design indicates that our count data is dependent on the exposure condition 
# (groups variable) and our factor of unwanted variation, and we have specified 
# that there not be an intercept term through the use of '~0'

## HAD TO ROUND...NOT SURE IF THAT'S OK
dds <- DESeqDataSetFromMatrix(countData = round(counts(ruv_set)), # Grabbing count data from the 'ruv_set' object
                              colData = pData(ruv_set), # Grabbing the phenotype data and corresponding factor of unwanted variation from the 'ruv_set' object
                              design = ~0+groups+W_1) # Setting up the statistical formula (see below)

# Estimate size factors from the dds object that was just created as the experiment above
dds <- estimateSizeFactors(dds)
sizeFactors(dds)  # viewing the size factors


# Calculating normalized count data
normcounts <- as.data.frame(counts(dds, normalized=TRUE))

# Transforming normalized counts through variance stabilization
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_matrix <- as.matrix(assay(vsd))
