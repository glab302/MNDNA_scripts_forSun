#!/usr/bin/env Rscript
source('utils.r')
## set the right library paths
library(optparse)
option_list = list(
    make_option(c("-p", "--path"), type = "character", default = "",  help = "a file of weeks data(row1:a group of mus(diff weeks); row2:another group of mus(diff weeks))"),
    make_option(c("-i", "--files"), type = "character", default = "NormalizedData/Mus1.8m.nodup.q30.bam.10kb_copyNumbersCorrect.t.NorDataSetSamplesInfo.log",  help = "a file of weeks data(row1:a group of mus(diff weeks); row2:another group of mus(diff weeks))"),
    make_option(c("-t", "--samples_type1"), type = "character", default = "A_Ctrl",  help = "a file of weeks data(row1:a group of mus(diff weeks); row2:another group of mus(diff weeks))"),
    make_option(c("-q", "--samples_type2"), type = "character", default = "A_Apc",  help = "a file of weeks data(row1:a group of mus(diff weeks); row2:another group of mus(diff weeks))"),
    make_option(c("-f", "--foldchange_or_not"), type = "character", default = "yes",  help = "Option is yes/ no"),
    make_option(c("-s", "--test.use"), type = "character", default = "anova",  help = "Option is anova/ wilcoxon/ null"),
    make_option(c("-a", "--pvalue.threshold"), type = "double", default = 0.01,  help = "Pvalue cutoff for significant regions"),
    make_option(c("-d", "--lg2fc.threshold"), type = "integer", default = 1,  help = "Option is 1"),
    make_option(c("-o", "--output"), type = "character", default = "2.ApcWT",  help = "Prefix of all files name")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)
path <- as.character(opt$path)
###########################################################################
###########################################################################
files <- opt$files
files.temp <- read.table(str_c(path,"/",as.character(files)),head=T,sep="\t", row.names=1)
files.temp1 <- as.data.frame(t(files.temp))
files.temp <- files.temp1
###########################################################################
########################################################################### 
# samples_type <- opt$samples_type
# files.temp <- files.temp[files.temp['task_label']==samples_type,]
# table(files.temp['sample_label'])
# samples_type1 <- files.temp[files.temp['sample_label']==as.character(opt$samples_type1), 1:(dim(files.temp)[2]-2)]
if(length(unlist(strsplit(opt$samples_type1, "\\|")))>1){
      samples_type1_tmp <- subset(files.temp, (sample_label==unlist(strsplit(opt$samples_type1, "\\|"))[1]) | (sample_label==unlist(strsplit(opt$samples_type1, "\\|"))[2]), select=c(-task_label, -sample_label))
}else{
      samples_type1_tmp <- subset(files.temp, sample_label==as.character(opt$samples_type1), select=c(-task_label, -sample_label))
}
samples_type1 <- samples_type1_tmp
# samples_type2 <- files.temp[files.temp['sample_label']==as.character(opt$samples_type2), 1:(dim(files.temp)[2]-2)]#files[,grep(samples_type_type2,colnames(files))]
if(length(unlist(strsplit(opt$samples_type2, "\\|")))>1){
      samples_type2_tmp <- subset(files.temp, (sample_label==unlist(strsplit(opt$samples_type2, "\\|"))[1]) | (sample_label==unlist(strsplit(opt$samples_type2, "\\|"))[2]), select=c(-task_label, -sample_label))
}else{
      samples_type2_tmp <- subset(files.temp, sample_label==as.character(opt$samples_type2), select=c(-task_label, -sample_label))
}
samples_type2 <- samples_type2_tmp
###########################################################################
#Option:#
#############################################intersect()##############################
foldchange_or_not <- as.character(opt$foldchange_or_not)
test.use <- as.character(opt$test.use)
pvalue.threshold <- as.numeric(opt$pvalue.threshold)
lg2fc.threshold <- as.numeric(opt$lg2fc.threshold)
Pseudocount=1/10000
###########################################################################
###########################################################################
setwd(str_c(path,"/UnsupervisedClustering/"))
title <- str_c(as.character(opt$output))
dir.create(title)
setwd(title)
###########################################################################
# Select regions and Default output:table and pvalueFC plot#
###########################################################################
D_samples_type1 <- as.data.frame(t(samples_type1))
D_samples_type2 <- as.data.frame(t(samples_type2))
D_samples_type2 <- D_samples_type2[, setdiff(colnames(D_samples_type2), c("A_yApc1", "A_yApc2", "A_yApc4","A_yApc5","A_yApc6"))]

D_samples_type1_1 <- c()
for(i in colnames(D_samples_type1)){
      D_samples_type1_1 <- cbind(D_samples_type1_1, as.numeric(as.character(as.matrix(D_samples_type1[, i]))))
}
D_samples_type1_1 <- as.data.frame(D_samples_type1_1)
colnames(D_samples_type1_1) <- colnames(D_samples_type1)
rownames(D_samples_type1_1) <- rownames(D_samples_type1)

D_samples_type2_1 <- c()
for(i in colnames(D_samples_type2)){
      D_samples_type2_1 <- cbind(D_samples_type2_1, as.numeric(as.character(as.matrix(D_samples_type2[, i]))))
}
D_samples_type2_1 <- as.data.frame(D_samples_type2_1)
colnames(D_samples_type2_1) <- colnames(D_samples_type2)
rownames(D_samples_type2_1) <- rownames(D_samples_type2)

grp2vgrp1 <- FindFeatures(D_samples_type1_1,D_samples_type2_1,test.use=test.use,lg2fc.threshold,pvalue.threshold,Pseudocount,title)
