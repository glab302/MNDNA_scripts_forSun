library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(scater)
library(metap)
library(paletteer)
library(ggsci)
library(Hmisc)
library(reshape2)
library(plyr)
library(org.Mm.eg.db)
library(rlist)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(GOplot)
library(ReactomeGSA.data)
library(ReactomeGSA)
library(ReactomePA)
mouse <- msigdbr(species = "Mus musculus")


data_summary <- function(data, varname, groupnames){
      require(plyr)
      summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum<-ddply(data, groupnames, .fun=summary_func,
                      varname)
      data_sum <- rename(data_sum, c("mean" = varname))
      return(data_sum)
    }



cellnumber_summary <- function(rds, sample_labels, sample_groups, celltypes, colors, name){
    library(reshape2)
    Idents(rds) <- factor(rds$celltype, levels=celltypes)
    cell_num_inclusters <- c()
    cellnumber <- c()
    cellratio_in_sample <- c()
    sample_label_list <- names(table(rds$sample_label))
    for(clusters_i in celltypes){
    	cellnum_in_cluster <- WhichCells(rds, idents = clusters_i)
    	sample_label_cell_num <- c()
    	cellratio_in_sample_label <- c()
    	for(sample_label_i in sample_label_list){
        # ## rename
        # cellnum_in_cluster <- gsub('BMT_20w_1', 'BMT_W_20w', cellnum_in_cluster)
        # cellnum_in_cluster <- gsub('BMT_20w_2', 'BMT_A_20w', cellnum_in_cluster)
            sample_label_cell_num <- c(sample_label_cell_num, length(grep(str_c('^', sample_label_i, '_'), cellnum_in_cluster)))
    		cellratio_in_sample_label <- c(cellratio_in_sample_label, length(grep(str_c('^', sample_label_i, '_'), cellnum_in_cluster))/table(rds$sample_label)[sample_label_i])
      }
    	cellnumber <- rbind(cellnumber, sample_label_cell_num)
    	cellratio_in_sample <- rbind(cellratio_in_sample, cellratio_in_sample_label)
    	cell_num_inclusters <- c(cell_num_inclusters, length(cellnum_in_cluster))
    }
    cells_statistics <- as.data.frame(cbind(cellnumber, cellnumber/cell_num_inclusters, cellratio_in_sample))
    colnames(cells_statistics) <- c(str_c('CellNum',sample_label_list), str_c('CellRatio_in_Clusters',sample_label_list), str_c('CellRatio_in_Samples',sample_label_list))
    rownames(cells_statistics) <- celltypes
    write.table(cells_statistics, str_c(name, 'cells_statistics.log'), sep='\t')
    # cells_statistics <- read.table('BoneMarrow_cellsratio.csv', sep=',')
    df3 <- melt(cbind(rownames(cells_statistics), cells_statistics[grep('CellRatio_in_Samples', colnames(cells_statistics))]))
    df3$label <- gsub("CellRatio_in_Samples", "", as.character(df3$variable))
    colnames(df3) <- c('celltype', 'variable', 'value', 'label')
    df3$label <- factor(df3$label, levels=sample_labels)
    df3$celltype <- factor(df3$celltype, levels=celltypes)
    df3$value100 <- df3$value*100
    df3$group <- gsub('_1$|_2$|_3$|_4$','', df3$label)
    df3$group <- factor(df3$group, sample_groups)
    p3 <- ggplot(data = df3, aes(x = group, y = value100, fill = group)) + ##
            geom_boxplot(outlier.shape = NA, lwd=0.1)+theme_bw()+
            geom_jitter(width = 0.3, size = 0.005)+#
            scale_fill_manual(values=c('#374E55FF','#E64B35FF'))+
            theme(legend.position='bottom',legend.title=element_blank(),axis.text.x = element_blank())+
            facet_grid(.~celltype)+
            stat_compare_means(label = "p.signif", method = "t.test",comparisons=list(c('W_20w','A_20w')))

    cells_statistics_insample <- cells_statistics[grep('CellRatio_in_Samples', colnames(cells_statistics))]
    colnames(cells_statistics_insample) <- gsub('CellRatio_in_Samples','',colnames(cells_statistics_insample))
    # cells_statistics_insample <- cells_statistics_insample[, c("W_5w_1","W_5w_2","W_20w_1","W_20w_2","W_20w_3","W_20w_4","A_5w_1","A_5w_2","A_20w_1","A_20w_2","A_20w_3","A_20w_4","BMT_W_20w","BMT_A_20w")]
    # cells_statistics_insample <- cells_statistics_insample[, c("W_20w_1","W_20w_2","A_20w_1","A_20w_2")]
    # pheatmap(cells_statistics_insample, scale='row')#, cluster_col=FALSE)
    # cells_statistics_insample <- cells_statistics_insample[, c("W_20w_3","W_20w_4","A_20w_3","A_20w_4")]
    cells_statistics_insample <- cells_statistics_insample[, sample_labels]
    df3_trans <- data_summary(df3, varname="value", groupnames=c("celltype", "group"))
    # Convert dose to a factor variable
    df3_trans$celltype=as.factor(df3_trans$celltype)
    colors=c('#374E55FF','#E64B35FF')
    p4 <- ggplot(df3_trans, aes(x=celltype, y=value, fill=group)) + 
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9))+
      scale_fill_manual(values=colors)+ theme(legend.position='bottom')+theme_classic()+theme(legend.position='right', axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
    g <- plot_grid(p3,p4,ncol=1)
    ggsave(str_c("./", name, "cellsratio_InAnnotationClusters.barplot2.pdf"), g, width=10, height=6)
}

enrichment <- function(genename, path2save){
    pathway_list <- list()
    genes <- genename
    trans.conserved <- bitr(genes, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Mm.eg.db'); rownames(trans.conserved) = trans.conserved$ENTREZID
    ego <- enrichGO(gene = trans.conserved$ENTREZID,OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 1,readable= TRUE)
    if(length(ego)!=0){
        if(any(ego@result$p.adjust<0.05)){
            p1_go <- dotplot(ego, x = "GeneRatio", color = "p.adjust", size = "Count", showCategory = 20, label_format = 70)
            ggsave(str_c(path2save, gsub('/','',clusteri), '._DEGs_il18_across_conditions.GOpathway.pdf'), p1_go, width=8, height=7)
            write.table(as.data.frame(ego@result), str_c(path2save, gsub('/','',clusteri), '._DEGs_il18_across_conditions.GOpathway.log'), sep='\t', row.names=FALSE)
            pathway_list[['GO']] = as.data.frame(ego@result)
            pathway_list[['GO_plot']] = p1_go$data[order(p1_go$data$GeneRatio, decreasing=TRUE),]
        }
    }
    ego <- enrichKEGG(gene = trans.conserved$ENTREZID, organism="mmu",pvalueCutoff = 0.05,qvalueCutoff = 1)
    if(length(ego)!=0){
        if(any(ego@result$p.adjust<0.05)){
            ego@result$Description <- gsub(' - Mus musculus \\(house mouse\\)','', ego@result$Description)
            p1_kegg <- dotplot(ego, x = "GeneRatio", color = "p.adjust", size = "Count", showCategory = 20, label_format = 70)
            ggsave(str_c(path2save, gsub('/','',clusteri), '._DEGs_il18_across_conditions.KEGG.pdf'), p1_kegg, width=8, height=6)
            SYMBOL_gene <- c()
            for(gene_id in ego@result$geneID){
                SYMBOL_gene <- c(SYMBOL_gene, gsub('/$','',paste0(trans.conserved[unlist(strsplit(gene_id,'/')), 'SYMBOL'],'/',collapse='')))
            }
            ego@result$geneID <- SYMBOL_gene
            write.table(as.data.frame(ego@result), str_c(path2save, gsub('/','',clusteri), '._DEGs_il18_across_conditions.KEGG.log'), sep='\t', row.names=FALSE)
            pathway_list[['KEGG']] = as.data.frame(ego@result)
            SYMBOL_gene <- c()
            for(gene_id in p1_kegg$data$geneID){
                SYMBOL_gene <- c(SYMBOL_gene, gsub('/$','',paste0(trans.conserved[unlist(strsplit(gene_id,'/')), 'SYMBOL'],'/',collapse='')))
            }
            p1_kegg$data$geneID <- SYMBOL_gene
            pathway_list[['KEGG_plot']] = p1_kegg$data[order(p1_kegg$data$GeneRatio, decreasing=TRUE),]
        }
    }
    ego <- enrichPathway(gene=trans.conserved$ENTREZID, organism="mouse",pvalueCutoff = 0.05, readable=TRUE)
    if(length(ego)!=0){
        if(any(ego@result$p.adjust<0.05)){
            p1_pathway <- dotplot(ego, x = "GeneRatio", color = "p.adjust", size = "Count", showCategory = 20, label_format = 70)
            ggsave(str_c(path2save, gsub('/','',clusteri), '._DEGs_il18_across_conditions.Reactome.pdf'), p1_pathway, width=7, height=12)
            write.table(as.data.frame(ego@result), str_c(path2save, gsub('/','',clusteri), '._DEGs_il18_across_conditions.Reactome.log'), sep='\t', row.names=FALSE)
            pathway_list[['REACTOME']] = as.data.frame(ego@result)
            pathway_list[['REACTOME_plot']] = p1_pathway$data[order(p1_pathway$data$GeneRatio, decreasing=TRUE),]
        }
    }
    return(pathway_list)
}

