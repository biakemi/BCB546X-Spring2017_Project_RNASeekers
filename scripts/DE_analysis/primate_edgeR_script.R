#load the edgeR package
library(edgeR)
#Read in the targets file
targets <- readTargets("primates_target.txt", row.names="sample.id")
targets_species <- readTargets("primates_target_species.txt", row.names="sample.id")

#Read in the count data file
x <- read.delim("primate_readcounts.txt", row.names=1, stringsAsFactors=FALSE)

#Assign column names to the count data using the targets file
colnames(x) = targets$sample.id

#Generate a DGEList of data for furture use with edgeR functions one for sex 
#comparisons and another for by species comparisons
y <- DGEList(counts=x, group=targets$trt)
colnames(y) <- targets$sample.id

y_species <- DGEList(counts=x, group=targets_species$trt)
colnames(y_species) <- targets_species$sample.id

#Estimate the normalization factor for each sample using TMM normalization
y <-calcNormFactors(y)

y_species <-calcNormFactors(y_species)

#Estimate the dispersion parameter for negative binomial model for 
#each gene assuming equal dispersion parameter across all four treatment groups
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)
length(y$tagwise.dispersion)

y_species <- estimateCommonDisp(y_species, verbose = TRUE)
y_species <- estimateTagwiseDisp(y_species)

#Write the normalized count data to a .csv file
write.csv(y$pseudo.counts, "primate_norm_read_count.csv", row.names=T)

######DE Gene Comparisons################################################
#To determine how many genes are differentially expressed between pairs 
#at 5% level of FDR. 1st using an exact tests for differences between two 
#groups of negative-binomial counts

#Human male to chimpanzee male
etHMCM <- exactTest(y, pair = c("HM", "CM"))
sum(etHMCM$table[, 3] <= 0.05)

#Human female to chimpanzee female
etHFCF <- exactTest(y, pair = c("HF", "CF"))
sum(etHFCF$table[, 3] <= 0.05)

#Human male to rhesus male
etHMRM <- exactTest(y, pair = c("HM", "RM"))
sum(etHMRM$table[, 3] <= 0.05)

#Human female to rhesus female
etHFRF <- exactTest(y, pair = c("HF", "RF"))
sum(etHFRF$table[, 3] <= 0.05)

#Chimpanzee male to rhesus male
etCMRM <- exactTest(y, pair = c("CM", "RM"))
sum(etCMRM$table[, 3] <= 0.05)

#Chimpanzee female to rhesus female
etCFRF <- exactTest(y, pair = c("CF", "RF"))
sum(etCFRF$table[, 3] <= 0.05)

#All human to chimpanzee
etHC <- exactTest(y_species, pair = c("H", "C"))
sum(etHC$table[, 3] <= 0.05)

#All human to rhesus
etHR <- exactTest(y_species, pair = c("H", "R"))
sum(etHR$table[, 3] <= 0.05)

#All chimpanzee to rhesus
etCR <- exactTest(y_species, pair = c("C", "R"))
sum(etCR$table[, 3] <= 0.05)




















