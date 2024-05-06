source('Rscript/stat.r')

color_defination <- list.load(str_c(workpath,'/color_defination.json'))
immune.combined <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))

########## GO/KEGG/REACTOME pathway #############
### 1. APC vs WT
dir.create('./APCvsWT_in_celltype/')
rds = immune.combined
Idents(rds) <- rds$celltype
rds$celltype.sample_group <- paste(Idents(rds), rds$sample_group, sep = "_")
rds$celltype <- Idents(rds)
Idents(rds) <- "celltype.sample_group"
GOpathway <- c()
KEGGpathway <- c()
REACTOMEpathway <- c()
GOpathway_plot <- c()
KEGGpathway_plot <- c()
REACTOMEpathway_plot <- c()
findmarkers_merge <- c()
for(clusteri in names(table(rds$celltype))){
    GO_df <- c()
    KEGG_df <- c()
    REACTOME_df <- c()    
    DefaultAssay(rds) <- "RNA"
    tmp <- FindMarkers(rds, ident.1 = str_c(clusteri, '_A_20w'), ident.2 = str_c(clusteri, '_W_20w'), verbose = FALSE, only.pos = TRUE)
    write.table(tmp, str_c('./APCvsWT_in_celltype/',gsub('/', '', clusteri), '_positiveDEGs_across_conditions.log'), sep='\t')
    tmp$celltype = clusteri
    findmarkers_merge <- rbind(findmarkers_merge, tmp)

    enrichment_list = enrichment(rownames(tmp), './APCvsWT_in_celltype/')
    if(length(enrichment_list[['GO']])!=0){
        GO_df = enrichment_list[['GO']]; GO_df$celltype = clusteri
        GOpathway <- rbind(GOpathway, GO_df)  
        GO_df_plot = enrichment_list[['GO_plot']]; GO_df_plot$celltype = clusteri
        GOpathway_plot <- rbind(GOpathway_plot, GO_df_plot)
    }
    print(head(GO_df))
    print(head(GOpathway))
    if(length(enrichment_list[['KEGG']])!=0){
        KEGG_df = enrichment_list[['KEGG']]; KEGG_df$celltype = clusteri
        KEGGpathway <- rbind(KEGGpathway, KEGG_df)
        KEGG_df_plot = enrichment_list[['KEGG_plot']]; KEGG_df_plot$celltype = clusteri
        KEGGpathway_plot <- rbind(KEGGpathway_plot, KEGG_df_plot)
    }
    if(length(enrichment_list[['REACTOME']])!=0){
        REACTOME_df = enrichment_list[['REACTOME']]; REACTOME_df$celltype = clusteri
        REACTOMEpathway <- rbind(REACTOMEpathway, REACTOME_df)
        REACTOME_df_plot = enrichment_list[['REACTOME_plot']]; REACTOME_df_plot$celltype = clusteri
        REACTOMEpathway_plot <- rbind(REACTOMEpathway_plot, REACTOME_df_plot)
    }
}

### 2. IL18+ APC vs WT
dir.create('./APCvsWT_in_IL-18pos_celltype.posDEGs/')
rds = subset(immune.combined, il18signal_stage=='High') #immune.combined # 
Idents(rds) <- rds$celltype
rds$celltype.sample_group <- paste(Idents(rds), rds$sample_group, sep = "_")
rds$celltype <- Idents(rds)
Idents(rds) <- "celltype.sample_group"
GOpathway <- c()
KEGGpathway <- c()
REACTOMEpathway <- c()
GOpathway_plot <- c()
KEGGpathway_plot <- c()
REACTOMEpathway_plot <- c()
findmarkers_merge <- c()
for(clusteri in names(table(rds$celltype))){
    GO_df <- c()
    KEGG_df <- c()
    REACTOME_df <- c()    
    DefaultAssay(rds) <- "RNA"
    tmp <- FindMarkers(rds, ident.1 = str_c(clusteri, '_A_20w'), ident.2 = str_c(clusteri, '_W_20w'), verbose = FALSE, only.pos = TRUE)
    write.table(tmp, str_c('./APCvsWT_in_IL-18pos_celltype.posDEGs/',gsub('/', '', clusteri), '_positiveDEGs_across_conditions.log'), sep='\t')
    tmp$celltype = clusteri
    findmarkers_merge <- rbind(findmarkers_merge, tmp)

    enrichment_list = enrichment(rownames(tmp), './APCvsWT_in_IL-18pos_celltype.posDEGs/')
    if(length(enrichment_list[['GO']])!=0){
        GO_df = enrichment_list[['GO']]; GO_df$celltype = clusteri
        GOpathway <- rbind(GOpathway, GO_df)  
        GO_df_plot = enrichment_list[['GO_plot']]; GO_df_plot$celltype = clusteri
        GOpathway_plot <- rbind(GOpathway_plot, GO_df_plot)
    }
    if(length(enrichment_list[['KEGG']])!=0){
        KEGG_df = enrichment_list[['KEGG']]; KEGG_df$celltype = clusteri
        KEGGpathway <- rbind(KEGGpathway, KEGG_df)
        KEGG_df_plot = enrichment_list[['KEGG_plot']]; KEGG_df_plot$celltype = clusteri
        KEGGpathway_plot <- rbind(KEGGpathway_plot, KEGG_df_plot)
    }
    if(length(enrichment_list[['REACTOME']])!=0){
        REACTOME_df = enrichment_list[['REACTOME']]; REACTOME_df$celltype = clusteri
        REACTOMEpathway <- rbind(REACTOMEpathway, REACTOME_df)
        REACTOME_df_plot = enrichment_list[['REACTOME_plot']]; REACTOME_df_plot$celltype = clusteri
        REACTOMEpathway_plot <- rbind(REACTOMEpathway_plot, REACTOME_df_plot)
    }
}
