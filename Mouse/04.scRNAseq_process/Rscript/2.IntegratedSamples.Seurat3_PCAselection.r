#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(Seurat)
library(stringr)

option_list = list(
    make_option(c("-w", "--work_path"), type = "character", default = "/data2/xingyun/singlecell/202111_scrnaseq_Musspleen/4.CombinedWithOthers", help = "File path (path of output dataset)"),
    make_option(c("-l", "--sample_label_list"), type = "character", default = "rbcPD1,PD1,rbc,stressrbc_c,stressrbc_g",  help = "The list of sample_label using for every sample"),
    make_option(c("-f", "--input_file_list"),  type = "character", default = "spleen_rbcPD1_raw.rds,spleen_PD1_raw.rds,spleen_rbc_raw.rds,stressrbc_c_spleen_raw.rds,stressrbc_g_spleen_raw.rds",  help = "Path and full name of inputfile for disease"),
    make_option(c("-u", "--nFeature_RNA_upper_limit_list"), type = "character", default = "4000,4000,4500,4500,4500",  help = "Upper limit of nFeature_RNA"),
    make_option(c("-m", "--mt_upper_limit_list"), type = "character", default = "15,15,15,15,15",  help = "Upper limit of percent.mt"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "Combined_5samples_spleen", help = "Prefix of output filename")
    )
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)

sample_label_list <- as.character(opt$sample_label_list)
inputfile_list <- as.character(opt$input_file_list)
workpath <- as.character(opt$work_path)
output_prename <- as.character(opt$output_prefix)
nFeature_RNA_upper_limit_list <- as.character(opt$nFeature_RNA_upper_limit_list)
percent.mt_upper_limit_list <- as.character(opt$mt_upper_limit_list)

# 反馈读参结果
print(str_c("-w: ", workpath))
print(str_c("-o: ", output_prename))
print(str_c("-f: ", inputfile_list))
print(str_c("-u: ", nFeature_RNA_upper_limit_list))
print(str_c("-m: ", percent.mt_upper_limit_list))


############################################################
# 函数： filter_object 
single_sample <- function(inputfile, sample_label, nFeature_RNA_upper_limit, percent.mt_upper_limit){
	################# object ##################
	# Load the raw rds dataset
	pbmc <- readRDS(str_c(inputfile))
	pbmc$sample_label <- sample_label
	Idents(pbmc) <- sample_label
	pbmc@project.name <- sample_label
	pbmc@ meta.data$ orig.ident <- sample_label
	# Screening cells 
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_RNA_upper_limit & percent.mt < percent.mt_upper_limit)
	# print(str_c("### The original barcodes of pbmc object:"))
	# head(x = colnames(x = pbmc))
	pbmc <- RenameCells(object = pbmc, add.cell.id = sample_label)
	# print(str_c("### The new barcodes of pbmc object:"))
	# head(x = colnames(x = pbmc))
	pbmc <- NormalizeData(pbmc, verbose = FALSE)
	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
	# print("The pbmc object info: ")
	return(pbmc)
}


outputpath <- paste0(workpath, "/", output_prename, collapse = "")


################################################################
############################ Begin #############################
################################################################
setwd(workpath)
################# object input ##################
pbmc_list <- c()
for(i in 1:(length(strsplit(inputfile_list, split=",")[[1]]))){
	pbmc <- single_sample(strsplit(inputfile_list, split=",")[[1]][i], as.character(strsplit(sample_label_list, split=",")[[1]][i]), as.numeric(strsplit(nFeature_RNA_upper_limit_list, split=",")[[1]][i]), as.numeric(strsplit(percent.mt_upper_limit_list, split=",")[[1]][i]))
	print(str_c("===== pbmc ", i, " object is input! ====="))
	print(as.character(strsplit(sample_label_list, split=",")[[1]][i]))
	pbmc_list <- c(pbmc_list, pbmc)
}


################## Standard Workflow ###################
# Run the Standard workflow 
pbmc_remain <- pbmc_list
pbmc_remain[[1]] <- NULL
cca_merge <- merge(x = pbmc_list[[1]], y = pbmc_remain)
merge.list <- SplitObject(cca_merge, split.by = "sample_label")
merge.list <- lapply(X = merge.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merge.list)
# perform intergration
immune.anchors <- FindIntegrationAnchors(object.list = merge.list, anchor.features = features, dims = 1:10)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:10)
DefaultAssay(immune.combined) <- "integrated"

print("The immune.combined object info:")
print(immune.combined)
saveRDS(immune.combined, file = str_c(outputpath,"_Origin_Integrated.rds"))
print(str_c(output_prename, " origin rds has saved!"))

# visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

pdf(str_c(outputpath,"_3.1.DimPlot.pdf")) 
DimPlot(immune.combined, reduction = "pca")
dev.off()
pdf(str_c(outputpath,"_3.2.DimHeatmap.pdf"))
DimHeatmap(immune.combined, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()
pdf(str_c(outputpath,"_3.3.ElbowPlot.pdf"))
ElbowPlot(immune.combined)
dev.off()
