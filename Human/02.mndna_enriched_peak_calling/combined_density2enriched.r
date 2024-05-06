###
# rm(list=ls())
library(stringr)
library(Rmisc)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(scales)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(rtracklayer)
library(tidyverse)
library(GenomicRanges)
load('data/ref.RData') # genome_size, feature_1000kb_anno, feature_100kb_anno
load('data/HD_deep_10sample.RData')

addChromatinInfo <- function(annotated_bed, peak_feature, rbcDNA_med, gDNA_med){
    ### genome size
    genome_size <- genome_size[(genome_size$V1!='M')&(genome_size$V1!='X')&(genome_size$V1!='Y'), ]
    genome_size$V1 <- str_c('chr', genome_size$V1)
    genome_size$Genome_scaled_size <- scale(genome_size$V2)
    ### annotation from step4
    GenomicAnno = annotated_bed ###from step4 annotation
    GenomicAnno <- merge(GenomicAnno, genome_size, by.x='chr', by.y='V1')
    GenomicAnno$Genome_scaled_size_y <- -1
    GenomicAnno$genedensity_y <- -2
    GenomicAnno$genedensity_scale <- scale(GenomicAnno$overlappedGeneNumDensity)/10
    GenomicAnno$CFS_y <- -3
    GenomicAnno$CFS_region<-0; GenomicAnno[(!is.na(GenomicAnno$CFS))|(!is.na(GenomicAnno$ERFS)),'CFS_region'] <- 0.5
    # GenomicAnno$ERFS_y <- -4
    # GenomicAnno$ERFS_region<-0; GenomicAnno[!is.na(GenomicAnno$ERFS),'ERFS_region'] <- 1
    armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p","12q",
                "13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q",
                "19p", "19q","20p","20q","21p","21q","22p","22q")
   
    ### add peaks results (添加feature区域)
    peak_feature = peak_feature[order(peak_feature$count,decreasing=T),]
    peak_feature$region = str_c('chr',peak_feature$chromosome,':',peak_feature$start,'-',peak_feature$end)

    ### add fold change to label rbcDNA significant regions
    rbcDNA_med$fc = rbcDNA_med$median / gDNA_med$median
    rbcDNA_med$fc_trans = gDNA_med$median / rbcDNA_med$median
    rbcDNA_fea = rbcDNA_med[rbcDNA_med$fc>1.2,]
    rbcDNA_fea_region = rbcDNA_fea[, c('feature','chromosome','start','end')]
    colnames(rbcDNA_fea_region) = c('feature','chromosome','start','end')
    rbcDNA_fea_region.hg38 = GRanges(rbcDNA_fea_region) ### filtered 1Mb region
    top_feature_region = peak_feature[,c('region','chromosome','start','end')]
    colnames(top_feature_region) = c('feature','chromosome','start','end')
    top_feature_region.hg38 = GRanges(top_feature_region)
    top_feature_region.hg38 = top_feature_region.hg38[subjectHits(findOverlaps(rbcDNA_fea_region.hg38, top_feature_region.hg38))]
    top_feature_region.hg38 ### 12个样本最终剩下6653regions ### 10个样本最终剩下6567regions

    ### find overlaps to selected regions
    CN_smooth_r1_region = rbcDNA_med[, c('feature','chromosome','start','end')]
    CN_smooth_r1_region.hg38 = GRanges(CN_smooth_r1_region)
    CN_smooth_r1_region$broadPeak = 0
    CN_smooth_r1_region[queryHits(findOverlaps(CN_smooth_r1_region.hg38, top_feature_region.hg38)), 'broadPeak'] = 1
    CN_smooth_r1 = merge(GenomicAnno, CN_smooth_r1_region[, c('feature','broadPeak')], by.x = 'region', by.y = 'feature')
    CN_smooth_r1$broadPeak_y = -4
    CN_smooth_r1$arm <- factor(CN_smooth_r1$arm, levels=armlevels)
    CN_smooth_r1 <- CN_smooth_r1[order(CN_smooth_r1$chr, CN_smooth_r1$start), ]
    return(CN_smooth_r1)

}

feature = read.table(str_c('result/', 'result_d_0_merge_t_rbcDNA_c_gDNA.10samples.broadPeak.addsort.bed'), sep='\t', head=T)
rbcDNA_med = HD12_60m_1000k_MN
gDNA_med = HD12_60m_1000k_gDNA
CN_smooth_r1_1000k <- addChromatinInfo(feature_1000kb_anno, feature, rbcDNA_med, gDNA_med)

rbcDNA_med = HD12_60m_100k_MN
gDNA_med = HD12_60m_100k_gDNA
CN_smooth_r1_100k <- addChromatinInfo(feature_100kb_anno, feature, rbcDNA_med, gDNA_med)

save(CN_smooth_r1_1000k, CN_smooth_r1_100k, file = 'result/HD_deep_10sample_anno.RData')

