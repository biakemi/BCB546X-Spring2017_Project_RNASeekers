library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
# Read in the data
readcounts <- read_csv("primate_norm_read_count.csv")
names(readcounts)[1] <-c('EnsemblGeneID')

# Analysis of the gene RBM4 expression

# Extract genes with ensemble gene ID
RBM4 <- readcounts[readcounts$EnsemblGeneID=="ENSG00000173933",]
# Format the data to get summary statistics
RBM4_summary <- 
  RBM4 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c(rep("CF",6), rep("CM",6), rep("HF",6), rep("HM",6), rep("RF",6), rep("RM",6))) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
RBM4_summary$type <- factor(RBM4_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
ggplot(RBM4_summary,aes(type,mean))+geom_point(aes(color=sex),size=3)+
geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.1))+
  ggtitle("RBM4")+
  theme_bw()
ggsave("RBM4_expression.jpg")

# Analysis of the gene PAQR3 expression
# Extract genes with ensemble gene ID
PAQR3 <- readcounts[readcounts$EnsemblGeneID=="ENSG00000163291",]
# Format the data to get summary statistics
PAQR3_summary <- 
  PAQR3 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c(rep("CF",6), rep("CM",6), rep("HF",6), rep("HM",6), rep("RF",6), rep("RM",6))) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
PAQR3_summary$type <- factor(PAQR3_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
ggplot(PAQR3_summary,aes(type,mean))+geom_point(aes(color=sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.1))+
  ggtitle("PAQR3")+
  theme_bw()
ggsave("PAQR3_expression.jpg")

# Analysis of the gene PGM2 expression
# Extract genes with ensemble gene ID
PGM2 <- readcounts[readcounts$EnsemblGeneID=="ENSG00000169299",]
# Format the data to get summary statistics
PGM2_summary <- 
  PGM2 %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c(rep("CF",6), rep("CM",6), rep("HF",6), rep("HM",6), rep("RF",6), rep("RM",6))) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
PGM2_summary$type <- factor(PAQR3_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
ggplot(PGM2_summary,aes(type,mean))+geom_point(aes(color=sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.1))+
  ggtitle("PGM2")+
  theme_bw()
ggsave("PGM2_expression.jpg")


# Analysis of the gene KDM5C expression
# Extract genes with ensemble gene ID
KDM5C <- readcounts[readcounts$EnsemblGeneID=="ENSG00000126012",]
# Format the data to get summary statistics
KDM5C_summary <- 
  KDM5C %>% gather(sample, measurement, 2:37) %>%
  mutate(type=c(rep("CF",6), rep("CM",6), rep("HF",6), rep("HM",6), rep("RF",6), rep("RM",6))) %>%
  group_by(type) %>%
  summarise_each(funs(mean,se=sd(.)/sqrt(n())), measurement)%>%
  mutate(sex=c("F","M","F","M","F","M"))
# Change the order of the levels of factors to make the figure be consistent with the paper
KDM5C_summary$type <- factor(KDM5C_summary$type, levels = c("HM","HF","CM","CF","RM","RF"))
# Plot the expression data
ggplot(KDM5C_summary,aes(type,mean))+geom_point(aes(color=sex),size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.1))+
  ggtitle("KDM5C")+
  theme_bw()
ggsave("KDM5C_expression.jpg")



