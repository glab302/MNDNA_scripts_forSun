#!/usr/bin/env Rscript
.libPaths('/storage/gaoxiaofeiLab/yaoxingyun/cpu/Software/Rlibrary_4.2.1/')
library(optparse)
library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(metap)
options(future.globals.maxSize = 8000 * 1024^2)

option_list = list(
    make_option(c("-w", "--work_path"),  type = "character", default = NULL,  help = "File path (path of output dataset)"),
    make_option(c("-o", "--output_prefix"), type = "character", default = NULL,  help = "Prefix of output filename"),
    make_option(c("-p", "--pca_selection"), type = "integer", default = NULL,  help = "PCA dim for CCA"),
    make_option(c("-f", "--input_file"),  type = "character", default = NULL,  help = "Path and full name of inputfile")
    )

parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)



inputfile <- as.character(opt$input_file)
workpath <- as.character(opt$work_path)
output_prename <- as.character(opt$output_prefix)
outputpath <- paste0(workpath, "/", output_prename, collapse = "")
pca <- as.numeric(opt$pca_selection)

# 反馈读参结果
print(str_c("-w: ", workpath))
print(str_c("-o: ", output_prename))
print(str_c("-f: ", inputfile))
print(str_c("-p: ", pca))



################################################################
############################ Begin #############################
################################################################
dir.create(workpath)
setwd(workpath)


immune.combined <- readRDS(str_c(inputfile))


####################################################################
###################### t-SNE and Clustering ########################
####################################################################
resolution_test <- seq(0.4, 1.2, by = 0.2)
for(resolution in resolution_test){
	# create the subdir for t-SNE
	subworkpath <- str_c(workpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution")
    dir.create(subworkpath)
  	# changes the assay object
	DefaultAssay(immune.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	print('--------------ScaleData--------------')
	immune.combined <- ScaleData(immune.combined, verbose = FALSE)
	print('--------------RunPCA--------------')
	immune.combined <- RunPCA(immune.combined, npcs = pca, verbose = FALSE)
	print('--------------RunTSNE--------------')
	immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:pca)
	print('--------------FindNeighbors--------------')
	immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:pca)
	print('--------------FindClusters--------------')
	immune.combined <- FindClusters(immune.combined, resolution = resolution)
	# QC plot
	p1 <- DimPlot(object = immune.combined, pt.size = .1, group.by = "sample_label") + ggplot2::theme(legend.position = "top")
	p2 <- VlnPlot(object = immune.combined, features = "PC_1", group.by = "sample_label", pt.size = .1) + ggplot2::theme(legend.position = "top")
	p3 <- VlnPlot(object = immune.combined, features = "PC_2", group.by = "sample_label", pt.size = .1) + ggplot2::theme(legend.position = "top")
	plot_grid(p1, p2, p3, nrow = 1)
	ggsave(str_c(subworkpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution_QCplot.pdf"), width=15, height=5)
	# Calculating Cell numbers and ratio
	original_idents <- Idents(immune.combined)
	cell_num_inclusters <- c()
	cellnumber <- c()
	cellratio_in_sample <- c()
	sample_label_list <- names(table(immune.combined$sample_label))
	for(clusters_i in 0:(length(levels(immune.combined))-1)){
		cellnum_in_cluster <- WhichCells(immune.combined, idents = clusters_i)
		sample_label_cell_num <- c()
		cellratio_in_sample_label <- c()
		for(sample_label_i in sample_label_list){
	    	sample_label_cell_num <- c(sample_label_cell_num, length(grep(str_c('^', sample_label_i, '_'), cellnum_in_cluster)))
			cellratio_in_sample_label <- c(cellratio_in_sample_label, length(grep(str_c('^', sample_label_i, '_'), cellnum_in_cluster))/table(immune.combined$sample_label)[sample_label_i])
	    }
		cellnumber <- rbind(cellnumber, sample_label_cell_num)
		cellratio_in_sample <- rbind(cellratio_in_sample, cellratio_in_sample_label)
    	cell_num_inclusters <- c(cell_num_inclusters, length(cellnum_in_cluster))
	}
	cells_statistics <- as.data.frame(cbind(cellnumber, cellnumber/cell_num_inclusters, cellratio_in_sample))
	colnames(cells_statistics) <- c(str_c('CellNum',sample_label_list), str_c('CellRatio_in_Clusters',sample_label_list), str_c('CellRatio_in_Samples',sample_label_list))
	rownames(cells_statistics) <- levels(immune.combined)
	write.table(cells_statistics, str_c(subworkpath, "/", output_prename, "_cellsratio.csv"), quote=FALSE, sep=",")

	new.cluster.ids <- str_c(levels(immune.combined), " :", cell_num_inclusters)
	# new.cluster.ids <- str_c(sample_label_list, ': ', table(immune.combined$sample_label))
	names(new.cluster.ids) <- levels(immune.combined)
	# Visualization
	immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
	p3 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE,label.size =3)+ NoLegend()
	Idents(immune.combined)=original_idents
	p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "sample_label",label.size =10)+ ggplot2::theme(legend.position = "top")
	p2 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size =10)+ ggplot2::theme(legend.position = "top")
	p4 <- DimPlot(immune.combined, reduction = "tsne", split.by = "sample_label",label.size =10)+ ggplot2::theme(legend.position = "right")
	plot_grid(plot_grid(p1, p2, p3, nrow=1, align="h", rel_widths = c(5,5,5)), p4, ncol=1, align="v", rel_heights = c(4, 4))
	# ggarrange(ggarrange(p1,p2,p3,ncol=3), p4, nrow=2)
	ggsave(str_c(subworkpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution.pdf"),width=13,height=10)
	# save the rds file
	saveRDS(immune.combined, file = str_c(subworkpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution.rds"))
	print("======================= Combined Object: ========================")
	print(immune.combined)

	# save the loupe result (only for tsne)
    tsne_result <- immune.combined@reductions$tsne@cell.embeddings
    # tsne_result <- immune.combined[['tsne']]@cell.embeddings
    tsne_result <- as.data.frame(cbind(rownames(tsne_result),tsne_result))
    colnames(tsne_result)[1] <- "Barcode"
    write.table(tsne_result, str_c(subworkpath, "/", output_prename,"_t-SNE_", pca, "PCA_", resolution, "Resolution_Loupe.csv"), sep = ",", quote = FALSE, row.names = FALSE)

  	# changes the assay object
    DefaultAssay(immune.combined) <- "RNA"
	# find diffgenes from the whole dataset and visualization
    immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE)#, min.pct = 0.25, logfc.threshold = 0.25)
    # immune.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10)#, wt = avg_logFC)
    # Find conserved cell type markers
    for(cluster in names(table(Idents(immune.combined)))){
      cluster.markers <- FindConservedMarkers(immune.combined, ident.1 = cluster, grouping.var = "sample_label", verbose = FALSE)
      write.table(cluster.markers, str_c(subworkpath, "/", output_prename,"_t-SNE_", pca, "PCA_", resolution, "Resolution_ConservedMarkers_in",as.character(cluster),".txt"), sep = "\t", quote = FALSE)
    }

	# changes the assay object
    DefaultAssay(immune.combined) <- "integrated"
    DoHeatmap(immune.combined, features = top10$gene) #+ NoLegend()
    ggsave(str_c(subworkpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution_DoHeatmap_top10.pdf"), width = 20, height = 20)
    write.table(top10, str_c(subworkpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution_DiffGenes_for_DoHeatmap_top10.csv"), sep = ",", quote = FALSE)
}


####################################################################
###################### UMAP and Clustering #########################
####################################################################
resolution_test <- seq(0.4, 1.2, by = 0.2)
for(resolution in resolution_test){
	# create the subdir for t-SNE
	subworkpath <- str_c(workpath, "/", output_prename, "_UMAP_", pca, "PCA_", resolution, "Resolution")
    dir.create(subworkpath)
	# get the two samples name
	control_name <- unique(immune.combined@meta.data$sample_label)[1]
	treat_name <- unique(immune.combined@meta.data$sample_label)[2]
  	# changes the assay object
	DefaultAssay(immune.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	immune.combined <- ScaleData(immune.combined, verbose = FALSE)
	immune.combined <- RunPCA(immune.combined, npcs = pca, verbose = FALSE)
	immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:pca)
	immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:pca)
	immune.combined <- FindClusters(immune.combined, resolution = resolution)
	# QC plot
	p1 <- DimPlot(object = immune.combined, pt.size = .1, group.by = "sample_label") + ggplot2::theme(legend.position = "top")
	p2 <- VlnPlot(object = immune.combined, features = "PC_1", group.by = "sample_label", pt.size = .1) + ggplot2::theme(legend.position = "top")
	p3 <- VlnPlot(object = immune.combined, features = "PC_2", group.by = "sample_label", pt.size = .1) + ggplot2::theme(legend.position = "top")
	plot_grid(p1, p2, p3, nrow = 1)
	ggsave(str_c(subworkpath, "/", output_prename, "_UMAP_", pca, "PCA_", resolution, "Resolution_QCplot.pdf"), width=15, height=5)
	# t-SNE and Clustering
	# Calculating Cell numbers and ratio
	original_idents <- Idents(immune.combined)
	cell_num_inclusters <- c()
	cellnumber <- c()
	cellratio_in_sample <- c()
	sample_label_list <- names(table(immune.combined$sample_label))
	for(clusters_i in 0:(length(levels(immune.combined))-1)){
		cellnum_in_cluster <- WhichCells(immune.combined, idents = clusters_i)
		sample_label_cell_num <- c()
		cellratio_in_sample_label <- c()
		for(sample_label_i in sample_label_list){
	    	sample_label_cell_num <- c(sample_label_cell_num, length(grep(sample_label_i, cellnum_in_cluster)))
			cellratio_in_sample_label <- c(cellratio_in_sample_label, length(grep(sample_label_i, cellnum_in_cluster))/table(immune.combined$sample_label)[sample_label_i])
	    }
		cellnumber <- rbind(cellnumber, sample_label_cell_num)
		cellratio_in_sample <- rbind(cellratio_in_sample, cellratio_in_sample_label)
    	cell_num_inclusters <- c(cell_num_inclusters, length(cellnum_in_cluster))
	}
	cells_statistics <- as.data.frame(cbind(cellnumber, cellnumber/cell_num_inclusters, cellratio_in_sample))
	colnames(cells_statistics) <- c(str_c('CellNum',sample_label_list), str_c('CellRatio_in_Clusters',sample_label_list), str_c('CellRatio_in_Samples',sample_label_list))
	rownames(cells_statistics) <- levels(immune.combined)
	write.table(cells_statistics, str_c(subworkpath, "/", output_prename, "_cellsratio.csv"), quote=FALSE, sep=",")

	new.cluster.ids <- str_c(levels(immune.combined), " :", cell_num_inclusters)
	# new.cluster.ids <- str_c(sample_label_list, ': ', table(immune.combined$sample_label))
	names(new.cluster.ids) <- levels(immune.combined)
	# Visualization
	immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
	p3 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,label.size =3)+ NoLegend()
	Idents(immune.combined)=original_idents
	p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample_label",label.size =10)+ ggplot2::theme(legend.position = "top")
	p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, label.size =10)+ ggplot2::theme(legend.position = "top")
	p4 <- DimPlot(immune.combined, reduction = "umap", split.by = "sample_label",label.size =10)+ ggplot2::theme(legend.position = "right")
	plot_grid(plot_grid(p1, p2, p3, nrow=1, align="h", rel_widths = c(5,5,5)), p4, ncol=1, align="v", rel_heights = c(4, 4))
	# ggarrange(ggarrange(p1,p2,p3,ncol=3), p4, nrow=2)
	ggsave(str_c(subworkpath, "/", output_prename, "_t-SNE_", pca, "PCA_", resolution, "Resolution.pdf"),width=13,height=10)
	# save the rds file
	saveRDS(immune.combined, file = str_c(subworkpath, "/", output_prename, "_UMAP_", pca, "PCA_", resolution, "Resolution.rds"))
	print("======================= Combined Object: ========================")
	print(immune.combined)

	# save the loupe result (only for umap)
    umap_result <- immune.combined@reductions$umap@cell.embeddings
    # umap_result <- immune.combined[['umap']]@cell.embeddings
    umap_result <- as.data.frame(cbind(rownames(umap_result),umap_result))
    colnames(umap_result)[1] <- "Barcode"
    write.table(umap_result, str_c(subworkpath,"/", output_prename,"_Loupe_UMAP_", pca, "PCA_", resolution, "Resolution_result.csv"), sep = ",", quote = FALSE, row.names = FALSE)

  	# changes the assay object
    DefaultAssay(immune.combined) <- "RNA"
	# find diffgenes from the whole dataset and visualization
    immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE)#, min.pct = 0.25, logfc.threshold = 0.25)
    # immune.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10)#, wt = avg_logFC)
    # Find conserved cell type markers
    for(cluster in names(table(Idents(immune.combined)))){
      cluster.markers <- FindConservedMarkers(immune.combined, ident.1 = cluster, grouping.var = "sample_label", verbose = FALSE)
      write.table(cluster.markers, str_c(subworkpath, "/", output_prename,"_UMAP_", pca, "PCA_", resolution, "Resolution_ConservedMarkers_in",as.character(cluster),".txt"), sep = "\t", quote = FALSE)
    }

  	# changes the assay object
    DefaultAssay(immune.combined) <- "integrated"
    DoHeatmap(immune.combined, features = top10$gene) #+ NoLegend()
    ggsave(str_c(subworkpath, "/", output_prename, "_UMAP_", pca, "PCA_", resolution, "Resolution_DoHeatmap_top10.pdf"), width = 20, height = 20)
    write.table(top10, str_c(subworkpath, "/", output_prename, "_UMAP_", pca, "PCA_", resolution, "Resolution_DiffGenes_for_DoHeatmap_top10.csv"), sep = ",", quote = FALSE)
}

