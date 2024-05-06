source('Rscript/stat.r')

color_defination <- list()
immune.combined <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.rds'))
immune.combined$sample_group <- gsub('_1$|_2$|_3$|_4$', '', immune.combined$sample_label)
print(table(immune.combined$sample_group))
color_defination['sample_group'] <- list(c('W_20w','A_20w'))
color_defination['sample_group_color'] <- list(c('#4c85bb','#FB6A4A'))
sample_labels <- c("W_20w_1","W_20w_2","W_20w_3","W_20w_4","A_20w_1","A_20w_2","A_20w_3","A_20w_4")
sample_groups <- color_defination$sample_group

#### 1. cellannotation by classic markers
# 1.1 original clusters:
p1 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, raster=FALSE, cols=c(paletteer_d("LaCroixColoR::PeachPear", n = length(table(Idents(immune.combined))), type = "continuous"))) #+ NoLegend()
ggsave('Anno1_SeuratClusters_withoutAnnotation.pdf', p1, width=7, height=5)

DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- factor(immune.combined$seurat_clusters, levels=c(12,9,24,8,21,23,10,3,2,0,15,6,25,18,16,14,26,4,11,1,7,29,5,22,19,17,13,28,20,27))
markers <- c('Kit','Cd34','Hlf','Mpo','Elane','Ms4a3','Fcgr3','Spi1','Camp','Ltf','Itgam','Ly6g','S100a8','S100a9','S100a6','Ly6c2','S100a4','Csf1r','Ly86','Adgre1','C1qb','Mrc1',#'Ms4a2','Lmo4',#,'Prtn3','Ctsg','Ms4a3''Fcer1a','Prss34',
  'Siglech','Ms4a2','Prss34','Pax5','Vpreb1', 'Ebf1','Cd79a','Il7r','Cd3g','Cd3e','Klrd1','Ptprc',
  'Apoe','Car2','Car1','Zfpm1','Gata1','Klf1','Tfrc','Hbb-bt','Hbb-bs','Fos','Fosb','Jun','Junb','Jund',
  'Cdh5','Col1a1','Kitl','Spp1')
p1 <- DotPlot(immune.combined, features = markers, dot.scale = 3, col.min=0, cols=c('white','red','red4')) + RotatedAxis() +theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave('Anno2_SeuratClusters_withoutAnnotation.dotplot.pdf', p1, width=11, height=6)

# 1.3 annotation clusters:
table(immune.combined$seurat_clusters)
immune.combined$celltype = dplyr::case_when(
  immune.combined$seurat_clusters %in% c(12) ~ "HSC/MPP",
  immune.combined$seurat_clusters %in% c(9) ~ "CMP", 
  immune.combined$seurat_clusters %in% c(17,19) ~ "MEP", 
  immune.combined$seurat_clusters %in% c(13,28) ~ "Erythroid",
  immune.combined$seurat_clusters %in% c(27,20) ~ "Unknown", 
  immune.combined$seurat_clusters %in% c(8,24,21) ~ 'GMP',
  immune.combined$seurat_clusters %in% c(23) ~ "Basophil", 
  immune.combined$seurat_clusters %in% c(0,2,3,10,21) ~ "Neutrophil",
  immune.combined$seurat_clusters %in% c(15,6,18,25) ~ "Monocyte",
  immune.combined$seurat_clusters %in% c(16) ~ "Macrophage",
  immune.combined$seurat_clusters %in% c(14) ~ "DC",
  immune.combined$seurat_clusters %in% c(22) ~ "NK",
  immune.combined$seurat_clusters %in% c(5) ~ "T",
  immune.combined$seurat_clusters %in% c(26) ~ "CLP", 
  immune.combined$seurat_clusters %in% c(4,11,1,7,29) ~ "B")
table(immune.combined$celltype)

### 注释list：
celltype_levels = c('HSC/MPP','CMP','GMP','Monocyte','Macrophage','DC','Neutrophil','Basophil','MEP','Erythroid','Unknown','CLP','B','NK','T')#,'Niche Cells')
celltype_colorlevels = c('#E9E4A6FF',paletteer_d("LaCroixColoR::PeachPear", n = 30, type = "continuous")[15:21],'#E9AE84FF','#F7582AFF','grey','#E9C593FF','#B09C85FF','#7E6148FF','#DB9550FF')#,'#E08D58FF')#,'grey')
color_defination['major_celltype'] <- list(celltype_levels)
color_defination['major_celltype_color'] <- list(celltype_colorlevels)

immune.combined$celltype = factor(immune.combined$celltype, levels=color_defination$major_celltype)
Idents(immune.combined) = immune.combined$celltype
p_cellanno <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel=TRUE, raster=FALSE, label.size = 3, cols=color_defination$major_celltype_color)+theme_bw()
ggsave('Anno1_SeuratClusters_withAnnotation.pdf', p_cellanno, width=6, height=4)

DefaultAssay(immune.combined) <- "RNA"
celltype_levels = c('HSC/MPP','CMP','GMP','Neutrophil','Monocyte','Macrophage','DC','Basophil','CLP','B','T','NK','MEP','Erythroid','Unknown')#,'Niche Cells')
Idents(immune.combined) <- factor(immune.combined$celltype, levels=celltype_levels)
p1 <- DotPlot(immune.combined, features = markers, col.min = -2.5, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis() +theme_bw()+
        theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
ggsave('Anno2_SeuratClusters_withAnnotation.dotplot.pdf', p1, width=10, height=4)

### 4. statistic summary:
sample_labels = c("W_20w_1","W_20w_2","W_20w_3","W_20w_4","A_20w_1","A_20w_2","A_20w_3","A_20w_4")
sample_groups = c('W_20w','A_20w')
celltypes = c('HSC/MPP','CMP','GMP','Neutrophil','Monocyte','Macrophage','DC','Basophil','CLP','B','T','NK','MEP','Erythroid','Unknown')
colors = c('#E9E4A6FF',paletteer_d("LaCroixColoR::PeachPear", n = 30, type = "continuous")[15:21],'#E9AE84FF','#F7582AFF','grey','#E9C593FF','#B09C85FF','#7E6148FF','#DB9550FF')#,'#E08D58FF')#,'grey')
cellnumber_summary(immune.combined, sample_labels, sample_groups, celltypes, colors, 'BM')


saveRDS(immune.combined, file = str_c(workpath,'/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))
list.save(color_defination, file = str_c(workpath,'/color_defination.json'))#"va1_result.txt")
write.table(cbind(immune.combined@meta.data, immune.combined@ reductions$ umap@ cell.embeddings), '/Users/xingyun/Documents/Data Report/DataReport/DataReport/【Mice BM scRNAseq】wt+apcmin/Final_20w_WTAPC/Supplementarytables_immune.combined.log', sep='\t')
print(str_c("\n===== CellAnnotation is finished! =====\n"))
