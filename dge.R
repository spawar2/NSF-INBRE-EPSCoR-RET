# load package ------------------------------------------------------------
library(edgeR)
library(hgu133plus2.db)
library(tidyverse)

# add annotation information ----------------------------------------------
# specify processed cancer data
processed_data <- gastric_cel_processed

# view annotation of processed data
pData(processed_data)

# add metadata
cancer <- rep("cancer", 10)
normal <- rep("normal", 10)

# add metadata to processed data (phenoData)
pData(processed_data)[, 1] <- c(cancer, normal)
pData(processed_data) <- rownames_to_column(pData(processed_data))

# rename columns
colnames(pData(processed_data)) <- c("samples", "condition")

# view annotation of processed data
pData(processed_data)

# store phenoData
phenodata <- pData(processed_data)

# set condition as a factor
phenodata$condition <- factor(phenodata$condition)

# create study design -----------------------------------------------------
# design
design <- model.matrix(~0+phenodata$condition)

# rename columns of design
colnames(design) <- c("cancer", "normal")

# fit linear model --------------------------------------------------------
fit <- lmFit(processed_data, design)

# make contrasts ----------------------------------------------------------
contrast <- makeContrasts(cancer - normal, levels=design)

# fit linear models with contrast -----------------------------------------
fit2 <- contrasts.fit(fit,contrasts=contrast)

# moderated t-test with ebayes --------------------------------------------
fit2 <- eBayes(fit2)

# summary of results
summary(decideTests(fit2,lfc=1))


# filter expression -------------------------------------------------------
# all genes that are differentially expressed with a adjusted p-value of 
# less than 0.05, with at fold change of at least two (log fold change at least 1)
fit2_adj <- topTable(fit2,number=Inf,p.value = 0.05,lfc=1)


# save log fold change results ----------------------------------------------------
write.csv(fit2_adj, "gastric_cancer_fold_change.csv")

