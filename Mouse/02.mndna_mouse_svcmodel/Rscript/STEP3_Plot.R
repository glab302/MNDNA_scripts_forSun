#!/usr/bin/env Rscript
source('utils.r')
## set the right library paths
library(optparse)

option_list = list(
    make_option(c("-p", "--path"), type = "character", default = "",  help = "a file of weeks data(row1:a group of mus(diff weeks); row2:another group of mus(diff weeks))"),
    make_option(c("-i", "--files"), type = "character", default = "NormalizedData/CovInBins_1000079.t.NorDataSetSamplesInfo.log",  help = "a file of weeks data(row1:a group of mus(diff weeks); row2:another group of mus(diff weeks))"),
    make_option(c("-s", "--test.use"), type = "character", default = "anova",  help = "Option is anova/ wilcoxon/ null"),
    make_option(c("-a", "--pvalue.threshold"), type = "double", default = 0.01,  help = "Pvalue cutoff for significant regions"),
    make_option(c("-d", "--lg2fc.threshold"), type = "integer", default = 1,  help = "Option is 1"),
    make_option(c("-o", "--output"), type = "character", default = "2.Apc",  help = "Prefix of all files name")
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
transformat <- function(df){
  D_samples <- as.data.frame(t(df))
  D_samples_1 <- c()
  for(i in colnames(D_samples)){
        D_samples_1 <- cbind(D_samples_1, as.numeric(as.character(as.matrix(D_samples[, i]))))
  }
  D_samples_1 <- as.data.frame(D_samples_1)
  colnames(D_samples_1) <- colnames(D_samples); rownames(D_samples_1) <- rownames(D_samples)
  return(D_samples_1)
}

# 2.Apc
# features = read.table(str_c(path, '/UnsupervisedClustering/2.Apc/2.ApcDataMatrix_SignedFindMarkers_40regions.log'), sep='\t')
features = read.table(str_c(path, '/UnsupervisedClustering/2.Apc/2.ApcDataMatrix_SignedFindMarkers_40regions.log'), sep='\t')
samples_type1 <- files.temp[files.temp['sample_label']=='A_Ctrl', rownames(features)]#files[,grep(samples_type_type1,colnames(files))]
samples_type2 <- files.temp[files.temp['sample_label']=='A_Apc', rownames(features)]#files[,grep(samples_type_type1,colnames(files))]
print(str_c("GroupA-sample:",dim(samples_type1)[1],";  GroupB-sample:",dim(samples_type2)[1]))
groups_label = 'Con_APC'
pic_name = str_c(str_c(path, '/UnsupervisedClustering/2.Apc/pheatmap_Con_APC_name.pdf'))
pheatmap_f_g2(groups_label, transformat(samples_type1), transformat(samples_type2), pic_name)