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
library(corrplot)
load('data/ref.RData') # genome_size, feature_1000kb_anno, feature_100kb_anno
load('data/HD_deep_10sample.RData')
load('result/HD_deep_10sample_anno.RData')

### Plot: genomic distribution across the genome
    p0 <- ggplot() + #
        geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.1, linetype='dashed') + 
        geom_ribbon(data=HD12_60m_1000k_MN, aes(x=start, ymin = min, ymax = max), alpha=0.6, fill="#E64B35FF")+
        geom_line(data=HD12_60m_1000k_MN, aes(x=start, y=median, group=1), color="#E64B35FF", size=0.3)+
        theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0.5,2)
    p2 <- ggplot() + #
        geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.3, linetype='dashed') + 
        geom_ribbon(data=HD12_60m_1000k_gDNA, aes(x=start, ymin = min, ymax = max), alpha=0.6, fill="#3E60AA")+
        geom_line(data=HD12_60m_1000k_gDNA, aes(x=start, y=median, group=1), color="#3E60AA", size=0.1)+
        theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0.5,2)
    p1_anno <- ggplot(CN_smooth_r1_1000k)+
            # geom_tile(aes(x=start, y=Genome_scaled_size_y, fill = Genome_scaled_size) ,colour = NA)+
            # geom_tile(aes(x=start, y=genedensity_y, fill = genedensity_scale) ,colour = NA)+
            # geom_tile(aes(x=start, y=CFS_y, fill = CFS_region) ,colour = NA)+
            geom_tile(aes(x=start, y=broadPeak_y, fill = broadPeak) ,colour = NA)+
            scale_fill_gradient2(low = "white", high = "red4")+mytheme+
            facet_grid(.~arm,switch="x", space="free_x", scales="free")

    g <- plot_grid(p0, p2, p1_anno, ncol=1,rel_heights=c(1,1,0.6))
    ggsave('result/Fig1.60m_1000k.10samples.whole_genome.pdf', g, width=8, height=4)

    ### chr1:198000001-200100001
    p0 <- ggplot() + #
        geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.1, linetype='dashed') + 
        geom_ribbon(data=HD12_60m_100k_MN[(HD12_60m_100k_MN$chromosome==1)&(HD12_60m_100k_MN$start>198000001)&(HD12_60m_100k_MN$end<200100001), ], aes(x=start, ymin = min, ymax = max), alpha=0.5, fill="#E64B35FF")+
        geom_line(data=HD12_60m_100k_MN[(HD12_60m_100k_MN$chromosome==1)&(HD12_60m_100k_MN$start>198000001)&(HD12_60m_100k_MN$end<200100001), ], aes(x=start, y=median, group=1), color="#E64B35FF", size=0.3)+
        theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0,4)
    p2 <- ggplot() + #
        geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.3, linetype='dashed') + 
        geom_ribbon(data=HD12_60m_100k_gDNA[(HD12_60m_100k_gDNA$chromosome==1)&(HD12_60m_100k_gDNA$start>198000001)&(HD12_60m_100k_gDNA$end<200100001), ], aes(x=start, ymin = min, ymax = max), alpha=0.5, fill="#3E60AA")+
        geom_line(data=HD12_60m_100k_gDNA[(HD12_60m_100k_gDNA$chromosome==1)&(HD12_60m_100k_gDNA$start>198000001)&(HD12_60m_100k_gDNA$end<200100001), ], aes(x=start, y=median, group=1), color="#3E60AA", size=0.3)+
        theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0,4)
    p1_anno <- ggplot(CN_smooth_r1_100k[(CN_smooth_r1_100k$chr=='chr1')&(CN_smooth_r1_100k$start>198000001)&(CN_smooth_r1_100k$end<200100001), ])+
            geom_tile(aes(x=start, y=broadPeak_y, fill = broadPeak) ,colour = NA)+
            scale_fill_gradient2(low = "white", high = "red4")+mytheme+
            facet_grid(.~arm,switch="x", space="free_x", scales="free")

    g = plot_grid(p0, p2, p1_anno, ncol=1,rel_heights=c(1,1,0.2))
    ggsave('result/Fig1.60m_100k.chr1_198000001_200100001_zoomin.pdf', g, width=2.4, height=3)


### Plot: Figure 3b, pearson correlation
    a=merge(HD12_60m_1000k_MN[,c('feature','median')], CN_smooth_r1_1000k, by.x='feature', by.y='region')
    a=a[a$broadPeak==1, c('median','V2','overlappedGeneNumDensity')]; colnames(a) = c('MN-DNA\nread density','Chrom. size','Gene density')
    
    cor2genomic <- cor(a[,])
    cols <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(500))
    pdf('result/supFig3.60m_1000k.10samples.median_correlation2genomicregion.pdf', width=6, height=6)
    corrplot(cor2genomic, method = 'color', col=rev(COL2('RdBu', 50)), addCoef.col = 'black', cl.cex = 0.6, number.cex=0.7)
    dev.off()

