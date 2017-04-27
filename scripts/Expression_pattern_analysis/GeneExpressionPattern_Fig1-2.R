# Install and load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
# Version information: R version 3.3.3 (2017-03-06) ggplot2_2.2.1 tidyr_0.6.1   dplyr_0.5.0   readr_1.1.0

########################################################
# Analysis using the normalization method in the paper #
########################################################

# Read in raw data of absolute counts
readcounts_raw <- read_delim("primate_readcounts.txt","\t", escape_double = FALSE,trim_ws = TRUE)
names(readcounts_raw)[1] <-c('EnsemblGeneID')
# Normalize readcounts using the method in the paper, which is the square root of the row proportions 
# use prop.table to get proportion
readcounts_prop <- cbind(readcounts_raw[1],prop.table(as.matrix(readcounts_raw[-1]),2))
readcounts_prop_sqrt <- cbind(readcounts_raw[1],sqrt(prop.table(as.matrix(readcounts_raw[-1]),2)))

# Analysis of the gene PEX14 expression
# Extract genes with ensemble gene ID
PEX14 <- readcounts_prop_sqrt[readcounts_prop_sqrt$EnsemblGeneID=="ENSG00000142655",]
PEX14_summary <- 
  PEX14 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c("H","C","R","H","C","R","R","H","C","R","H","C","R","H","C","R",
                "H","C","R","H","H","R","C","C","R","H","R","H","C","R","R","C"
                ,"C","C","H","H")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)

PEX14_summary$type <- factor(PEX14_summary$type, levels = c("H","C","R"))
# Plot the expression data
# color by sex
cols <- c("H" = "red", "C" = "blue", "R" = "black")

ggplot(PEX14_summary,aes(type,mean))+geom_point(aes(colour= type),size=1)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= type), width=.2,position=position_dodge(0.1))+
  ggtitle("PEX14")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
  axis.text.y = element_text(colour="grey20",size=8,angle=90,hjust=1,vjust=0,face="plain"))+
  scale_colour_manual(values = cols)+
  scale_y_continuous(limits = c(0, 0.01))
ggsave("norm1_PEX4_expression.jpg",width = 3, height = 5, units = "in") 


# Analysis of the gene GPC1 expression
# Extract genes with ensemble gene ID
GPC1 <- readcounts_prop_sqrt[readcounts_prop_sqrt$EnsemblGeneID=="ENSG00000063660",]
GPC1_summary <- 
  GPC1 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c("H","C","R","H","C","R","R","H","C","R","H","C","R","H","C","R",
                "H","C","R","H","H","R","C","C","R","H","R","H","C","R","R","C"
                ,"C","C","H","H")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)

GPC1_summary$type <- factor(GPC1_summary$type, levels = c("H","C","R"))
# Plot the expression data
# color by sex
cols <- c("H" = "red", "C" = "blue", "R" = "black")

ggplot(GPC1_summary,aes(type,mean))+geom_point(aes(colour= type),size=1)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= type), width=.2,position=position_dodge(0.1))+
  ggtitle("GPC1")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=8,angle=90,hjust=1,vjust=0,face="plain"))+
  scale_colour_manual(values = cols)+
  scale_y_continuous(limits = c(0, 0.01))
ggsave("norm1_GPC1_expression.jpg",width = 3, height = 5, units = "in") 


# Analysis of the gene RBM4 expression
# Extract genes with ensemble gene ID
RBM4 <- readcounts_prop_sqrt[readcounts_prop_sqrt$EnsemblGeneID=="ENSG00000173933",]
# Format the data to get summary statistics
RBM4_summary <- 
  RBM4 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","HF","CM","RF",
                "HM","CF","RM","HM","HF","RM","CF","CM","RF","HM","RF","HM","CF","RM","RF","CM"
                ,"CM","CF","HF","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
RBM4_summary$type <- factor(RBM4_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
RBM4_summary$sex <- factor(RBM4_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")
                  
ggplot(RBM4_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("RBM4")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm1_RBM4_expression.jpg",width = 5, height = 5, units = "in") 

# Analysis of the gene PAQR3 expression
# Extract genes with ensemble gene ID
PAQR3<- readcounts_prop_sqrt[readcounts_prop_sqrt$EnsemblGeneID=="ENSG00000163291",]
# Format the data to get summary statistics
PAQR3_summary <- 
  PAQR3 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","HF","CM","RF",
                "HM","CF","RM","HM","HF","RM","CF","CM","RF","HM","RF","HM","CF","RM","RF","CM"
                ,"CM","CF","HF","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
PAQR3_summary$type <- factor(PAQR3_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
PAQR3_summary$sex <- factor(PAQR3_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(PAQR3_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("PAQR3")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm1_PAQR3_expression.jpg",width = 5, height = 5, units = "in") 


# Analysis of the gene PGM2 expression
# Extract genes with ensemble gene ID
PGM2<- readcounts_prop_sqrt[readcounts_prop_sqrt$EnsemblGeneID=="ENSG00000169299",]
# Format the data to get summary statistics
PGM2_summary <- 
  PGM2 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","HF","CM","RF",
                "HM","CF","RM","HM","HF","RM","CF","CM","RF","HM","RF","HM","CF","RM","RF","CM"
                ,"CM","CF","HF","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
PGM2_summary$type <- factor(PGM2_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
PGM2_summary$sex <- factor(PGM2_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(PGM2_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("PGM2")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm1_PGM2_expression.jpg",width = 5, height = 5, units = "in") 


# Analysis of the gene KDM5C expression
# Extract genes with ensemble gene ID
KDM3C<- readcounts_prop_sqrt[readcounts_prop_sqrt$EnsemblGeneID=="ENSG00000126012",]
# Format the data to get summary statistics
KDM3C_summary <- 
  KDM3C %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","HF","CM","RF",
                "HM","CF","RM","HM","HF","RM","CF","CM","RF","HM","RF","HM","CF","RM","RF","CM"
                ,"CM","CF","HF","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
KDM3C_summary$type <- factor(KDM3C_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
KDM3C_summary$sex <- factor(KDM3C_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(KDM3C_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("KDM3C")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm1_KDM3C_expression.jpg",width = 5, height = 5, units = "in") 

#################################################
# Analysis using the normalized data from edgeR #
#################################################

# This is the normalized read count data after the technical replicates have been pooled together(average), rows with zero reads are removed, and normalized using TMM. 

# Read in normalized data generated by edgeR, which uses a trimmed mean of Mvalues (TMM)
readcounts <- read_csv("primate_norm_read_count_pooled.csv")
names(readcounts)[1] <-c('EnsemblGeneID')

# Analysis of the gene RBM4 expression
# Extract genes with ensemble gene ID
RBM4 <- readcounts[readcounts$EnsemblGeneID=="ENSG00000173933",]
# Format the data to get summary statistics
RBM4_summary <- 
  RBM4 %>% gather(sample, measurement, 2:19) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","RF","HM","CF","CM","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
RBM4_summary$type <- factor(RBM4_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
RBM4_summary$sex <- factor(RBM4_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(RBM4_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("RBM4")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm2_RBM4_expression.jpg",width = 5, height = 5, units = "in") 


# Analysis of the gene PAQR3 expression
# Extract genes with ensemble gene ID
PAQR3 <- readcounts[readcounts$EnsemblGeneID=="ENSG00000163291",]
# Format the data to get summary statistics
PAQR3_summary <- 
  PAQR3 %>% gather(sample, measurement, 2:19) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","RF","HM","CF","CM","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
PAQR3_summary$type <- factor(PAQR3_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
PAQR3_summary$sex <- factor(PAQR3_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(PAQR3_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("PAQR3")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm2_PAQR3_expression.jpg",width = 5, height = 5, units = "in") 


# Analysis of the gene PGM2 expression
# Extract genes with ensemble gene ID
PGM2<- readcounts[readcounts$EnsemblGeneID=="ENSG00000169299",]
# Format the data to get summary statistics
PGM2_summary <- 
  PGM2 %>% gather(sample, measurement, 2:19) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","RF","HM","CF","CM","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
PGM2_summary$type <- factor(PGM2_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
PGM2_summary$sex <- factor(PGM2_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(PGM2_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("PGM2")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm2_PGM2_expression.jpg",width = 5, height = 5, units = "in") 


# Analysis of the gene KDM5C expression
# Extract genes with ensemble gene ID
KDM5C <- readcounts[readcounts$EnsemblGeneID=="ENSG00000126012",]
# Format the data to get summary statistics
KDM5C_summary <- 
  KDM5C %>% gather(sample, measurement, 2:19) %>%
  mutate(type=c("HM","CF","RM","HF","CM","RF","RF","HM","CF","RM","HF","CM","RM","RF","HM","CF","CM","HF")) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
KDM5C_summary$type <- factor(KDM5C_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
# color by sex
KDM5C_summary$sex <- factor(KDM5C_summary$sex)
cols <- c("F" = "springgreen1", "M" = "dodgerblue3")

ggplot(KDM5C_summary,aes(type,mean))+geom_point(aes(colour= sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour= sex), width=.2,position=position_dodge(0.1))+
  ggtitle("KDM5C")+
  ylab("Normalized Expression") + xlab("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"))+
  scale_colour_manual(values = cols)
ggsave("norm2_KDM5C_expression.jpg",width = 5, height = 5, units = "in") 



