# rm(list=ls())#
library(stringr)#
library(BSgenome.Hsapiens.UCSC.hg38)#
library(rtracklayer)#
library(tidyverse)#
library(httr)#
source('addGenomicInfo.R')#
load('feature/features.RData')

### for human:#
genome <- "hg38"#
mySession <- browserSession()#
genome(mySession) <- genome#
arms <- deal_Chr_arm(mySession)#
cyto.hg38 <- deal_Chr_band(mySession)#
cfs.hg38 <- deal_fragile()#
erfs.hg38 <- deal_ERFS()#
knownGene.hg38 <- deal_GeneAnno(mySession)

#########################################################################################################################################
### 10samples mnDNA+ 位点（6656个特征）#
filename='MNDNApos6566.bed'#
feature_bed = MNDNApos6566; feature_bed = feature_bed[, c('chr','start','end')]#
colnames(feature_bed) <- c('chromosome', 'start', 'end')#
add_GenomicAnnotation_userbed(feature_bed, 'result/', 'MNDNApos6566')#
#
### CRC features#
filename='CRC0106.overlap.bed'#
feature_bed = CRC0106; feature_bed = feature_bed[, c('chr','start','end')]#
colnames(feature_bed) <- c('chromosome', 'start', 'end')#
add_GenomicAnnotation_userbed(feature_bed, 'result/', 'CRC0106')#
#
### PanCancer features#
filename='PanC0106.overlap.bed'#
feature_bed = PanC0106; feature_bed = feature_bed[, c('chr','start','end')]#
colnames(feature_bed) <- c('chromosome', 'start', 'end')#
add_GenomicAnnotation_userbed(feature_bed, 'result/', 'PanC0106')


### for mouse:
genome <- "mm10"#
session <- browserSession(url="https://genome.ucsc.edu/cgi-bin/")#
mySession  <- session#
genome(mySession) <- genome#
### add annotation: # genomeSession = mySession#
arms <- deal_mouseChr_arm(mySession)#
cyto.mm10 <- deal_mouseChr_band(mySession)#
knownGene.mm10 <- deal_mouseGeneAnno(mySession)#
### apc significantly different regions#
filename = 'ApcDataMatrix_SignedFindMarkers.bed'#
feature_bed = ApcDataMatrix_SignedFindMarkers#
feature_bed_anno = add_mm10GenomicAnnotation_userbed(feature_bed, 'result/', 'ApcDataMatrix_SignedFindMarkers')
### 40#
filename = 'ApcDataMatrix_SignedFindMarkers_40regions.bed'#
feature_bed = ApcDataMatrix_SignedFindMarkers_40regions#
feature_bed_anno = add_mm10GenomicAnnotation_userbed(feature_bed, 'result/', 'ApcDataMatrix_SignedFindMarkers_40')
