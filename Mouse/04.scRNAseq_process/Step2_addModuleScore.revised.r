source('Rscript/stat.r')

color_defination <- list.load(str_c(workpath,'/color_defination.json'))
immune.combined <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))

### 1. cellcycle in annotation clusters
color_defination['Phase'] <- list(c("G1","G2M","S"))
color_defination['Phase_color'] <- list(c("#0E5F98FF","#4CBBAAFF","#C6C7A0FF"))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- c(capitalize(tolower(cc.genes$s.genes)), 'Cenpu') # Cenpu: Mlf1ip(human)
g2m.genes <- c(capitalize(tolower(cc.genes$g2m.genes)), 'Pimreg','Jpt1') # Pimreg: Fam64a(Human); Jpt1: Hn1(Human)
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

### 2. il18相关(2019 paper)
library(readxl)
DefaultAssay(immune.combined) <- "RNA"
il18signal_gene <- readxl::read_xlsx('Downloads/12079_2019_544_MOESM6_ESM.xlsx')
il18signal_gene <- gsub('\r\n', '', setdiff(as.data.frame(il18signal_gene)[,1],'Gene name'))
il18signal_gene <- capitalize(tolower(il18signal_gene))
gene <- list()
gene$gene <- il18signal_gene#il18_mediated_signaling_pathway_mm10
immune.combined <- AddModuleScore(object=immune.combined, features=gene, ctrl=100, name='CD_features')
colnames(immune.combined@meta.data)[length(colnames(immune.combined@meta.data))] <- 'il18signal'
immune.combined[['il18signal_stage']] <- ifelse(immune.combined@meta.data[,'il18signal'] > mean(immune.combined@meta.data[,'il18signal']), 'High', 'Low')

saveRDS(immune.combined, file = str_c(workpath,'/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))
list.save(color_defination, file = str_c(workpath,'/color_defination.json'))#"va1_result.txt")
print(str_c("\n===== CellAnnotation is finished! =====\n"))