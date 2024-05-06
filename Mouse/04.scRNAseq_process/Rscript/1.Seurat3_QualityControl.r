#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(Seurat)
library(stringr)

option_list = list(
    make_option(c("-i", "--input_path"), type = "character", default = NULL, help = "File path (path of input dataset)"),
    make_option(c("-w", "--work_path"), type = "character", default = NULL, help = "File path (path of output dataset)"),
    make_option(c("-s", "--species"), type = "character", default = NULL, help = "File path (path of output dataset)"),
    make_option(c("-o", "--output_prefix"), type = "character", default = NULL, help = "Prefix of output filename")
    )

parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)


inputpath <- as.character(opt$input_path)
workpath <- as.character(opt$work_path)
species <- as.character(opt$species)
output_prename <- as.character(opt$output_prefix)
outputpath <- paste0(workpath, "/", output_prename, collapse = "")

# 反馈读参结果
print(str_c("-i: ", inputpath))
print(str_c("-w: ", workpath))
print(str_c("-o: ", output_prename))
# data1="/storage/gaoxiaofeiLab/yaoxingyun/cpu/yaoxingyun/singlecell/202212_scrnaseq_HumanPatient/1.rawdata/LC-X20221201013/D010101_504h_multi_result/outs/per_sample_outs/D010101_504h_multi_result/count/sample_filtered_feature_bc_matrix/"
# sample1="D010101_504h"
# mkdir $sample1
# workpath="/storage/gaoxiaofeiLab/yaoxingyun/cpu/yaoxingyun/singlecell/202212_scrnaseq_HumanPatient/Results/D010101_504h/"


################################################################
############################ Begin #############################
################################################################
dir.create(workpath)
setwd(workpath)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = str_c(inputpath, "/"))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = output_prename, min.cells = 3, min.features = 200)
pbmc$sample_label <- output_prename
print(str_c(output_prename, "_info : "))
print(pbmc)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
if(species=='Human'){
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
}else if(species=='Mus'){
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
}

# Visualize QC metrics as a violin plot
pdf(str_c(outputpath, "_1.Visualize_QC_metrics_as_a_violin_plot.pdf"), width = 10, height = 4)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(str_c(outputpath, "_2.FeatureScatter_QC.pdf"), width = 10, height = 5)
# CombinePlots(plots = list(plot1, plot2))
plot1 + plot2
dev.off()

# save the rds file
saveRDS(pbmc, file = str_c(outputpath, "_raw.rds"))
print(str_c(output_prename, " raw rds has saved!"))

