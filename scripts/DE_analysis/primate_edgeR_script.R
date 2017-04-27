#load the edgeR package
library(edgeR)
library(qvalue)
#Read in the targets file
targets <- readTargets("primates_target.txt", row.names="sample.id")
targets_species <- readTargets("primates_target_species.txt", row.names="sample.id")

#Read in the count data file with technical replicates pooled together 
x <- read.delim("primate_readcounts_pooled.txt", row.names=1, stringsAsFactors=FALSE)

#Assign column names to the count data using the targets file
colnames(x) = targets$sample.id

# data filtering to remove rows with counts below 1
keep <- rowSums(cpm(x) > 1) >= 2
x <- x[keep, ]

#Generate a DGEList of data for furture use with edgeR functions one for sex 
#comparisons and another for by species comparisons

group <- factor(targets$trt)
y <- DGEList(counts=x, group=group)
colnames(y) <- targets$sample.id

group <- factor(targets_species$trt)
y_species <- DGEList(counts=x, group=group)
colnames(y_species) <- targets_species$sample.id

#Estimate the normalization factor for each sample using TMM normalization
y <-calcNormFactors(y, method = "TMM")
y_species <-calcNormFactors(y_species, method = "TMM")

#Estimate the dispersion parameter for negative binomial model for 
#each gene assuming equal dispersion parameter across all four treatment groups
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)

y_species <- estimateCommonDisp(y_species, verbose = TRUE)
y_species <- estimateTagwiseDisp(y_species)

#Write the normalized count data to a .csv file
write.csv(y$pseudo.counts, "primate_norm_read_count_pooled.csv", row.names=T)

######DE Gene Comparisons################################################
#To determine how many genes are differentially expressed between pairs 
#at 5% level of FDR. To correct for multiple testing, the approach of Storey and 
#Tibshirani is used. 1st using an exact tests for differences between two 
#groups of negative-binomial counts

#Human male to chimpanzee male
etHMCM <- exactTest(y, pair = c("HM", "CM"))
qobjHMCM <- qvalue(etHMCM$table$PValue)
sum(qobjHMCM$qvalues<0.05)


#Human female to chimpanzee female
etHFCF <- exactTest(y, pair = c("HF", "CF"))
qobjHFCF <- qvalue(etHFCF$table$PValue)
sum(qobjHFCF$qvalues<0.05)

#Human male to rhesus male
etHMRM <- exactTest(y, pair = c("HM", "RM"))
qobHMRM <- qvalue(etHMRM$table$PValue)
sum(qobHMRM$qvalues<0.05)

#Human female to rhesus female
etHFRF <- exactTest(y, pair = c("HF", "RF"))
qobHFRF <- qvalue(etHFRF$table$PValue)
sum(qobHFRF$qvalues<0.05)

#Chimpanzee male to rhesus male
etCMRM <- exactTest(y, pair = c("CM", "RM"))
qobCMRM <- qvalue(etCMRM$table$PValue)
sum(qobCMRM$qvalues<0.05)

#Chimpanzee female to rhesus female
etCFRF <- exactTest(y, pair = c("CF", "RF"))
qobCFRF  <- qvalue(etCFRF$table$PValue)
sum(qobCFRF$qvalues<0.05)

#All human to chimpanzee
etHC <- exactTest(y_species, pair = c("H", "C"))
qobHC  <- qvalue(etHC$table$PValue)
sum(qobHC$qvalues<0.05)

#All human to rhesus
etHR <- exactTest(y_species, pair = c("H", "R"))
qobHR  <- qvalue(etHR$table$PValue)
sum(qobHR$qvalues<0.05)

#All chimpanzee to rhesus
etCR <- exactTest(y_species, pair = c("C", "R"))
qobCR  <- qvalue(etCR$table$PValue)
sum(qobCR$qvalues<0.05)




















