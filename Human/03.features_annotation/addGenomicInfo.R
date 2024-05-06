library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(tidyverse)
library(httr)

deal_Chr_arm <- function(genomeSession){
    ### hg38 gaps & blacklisted regions
    centromeres <- getTable(ucscTableQuery(genomeSession, table="centromeres"))
    centromeres.hg38 <- do.call(rbind,lapply(unique(centromeres$chrom),FUN = function(chr,centromeres){
                                data.frame(chrom=chr,
                                            chromStart=min(centromeres$chromStart[centromeres$chrom==chr]),
                                            chromEnd=max(centromeres$chromEnd[centromeres$chrom==chr]),
                                            type='centromere')},centromeres=centromeres))
    gaps <- getTable(ucscTableQuery(genomeSession, table="gap"))
    gap_centro <- rbind(gaps[,c('chrom','chromStart','chromEnd','type')], centromeres.hg38)
    gaps.hg38 <- GRanges(gap_centro$chrom, IRanges(gap_centro$chromStart, gap_centro$chromEnd), type=gap_centro$type)
    # 只保留chr1:22
    gaps.hg38 <- keepSeqlevels(gaps.hg38, paste0("chr", c(1:22)), pruning.mode="coarse")
    # seqinfo(gaps.hg38) <- seqinfo(hsapiens)[seqlevels(gaps.hg38),]
    tcmeres <- gaps.hg38[grepl("centromere|telomere", gaps.hg38$type)]
    # 在之前的染色体上去除这部分区间
    chromosomes <- GRanges(paste0("chr", 1:22), IRanges(0, seqlengths(Hsapiens)[1:22]))
    arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
    # arms <- arms[-c(25,27,29,41,43)]  ##删除"13p","14p","15p","21p",,"22p"
    # armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p", "19q","20p","20q","21q","22q")
    armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p","12q",
                "13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q",
                "19p", "19q","20p","20q","21p","21q","22p","22q")
    arms$arm <- armlevels
    
    tmp1 = as.data.frame(arms); colnames(tmp1)[6] = 'type'
    tmp2 = as.data.frame(tcmeres); colnames(tmp2)[6] = 'type'
    tmp = rbind(tmp1, tmp2)

    arms_tcmer = GRanges(tmp)
    return(arms_tcmer)
}
deal_mouseChr_arm <- function(genomeSession){
    gaps <- getTable(ucscTableQuery(genomeSession, table="gap"))
    # gap_centro <- rbind(gaps[,c('chrom','chromStart','chromEnd','type')], centromeres.mm10)
    gaps.mm10 <- GRanges(gaps$chrom, IRanges(gaps$chromStart, gaps$chromEnd), type=gaps$type)
    # 只保留chr1:22
    gaps.mm10 <- keepSeqlevels(gaps.mm10, paste0("chr", c(1:19)), pruning.mode="coarse")
    # seqinfo(gaps.mm10) <- seqinfo(hsapiens)[seqlevels(gaps.mm10),]
    tcmeres <- gaps.mm10[grepl("centromere|telomere", gaps.mm10$type)]
    # 在之前的染色体上去除这部分区间
    chromosomes <- GRanges(paste0("chr", 1:19), IRanges(0, seqlengths(Mmusculus)[1:19]))
    arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
    arms$arm <- 'chr'
    
    tmp1 = as.data.frame(arms); colnames(tmp1)[6] = 'type'
    tmp2 = as.data.frame(tcmeres); colnames(tmp2)[6] = 'type'
    tmp = rbind(tmp1, tmp2)

    arms_tcmer = GRanges(tmp)
    return(arms_tcmer)
}


deal_Chr_band <- function(genomeSession){
    ### cytoband # nolint # nolint # nolint
    cyto <- getTable(ucscTableQuery(genomeSession, table="cytoBand"))
    cyto.hg38 <- GRanges(cyto$chrom, IRanges(cyto$chromStart, cyto$chromEnd), type=cyto$name)
    # 只保留chr1:22
    cyto.hg38 <- keepSeqlevels(cyto.hg38, paste0("chr", c(1:22)), pruning.mode="coarse")
    return(cyto.hg38)
}
deal_mouseChr_band <- function(genomeSession){
    ### cytoband # nolint # nolint # nolint
    cyto <- getTable(ucscTableQuery(genomeSession, table="cytoBand"))
    cyto.mm10 <- GRanges(cyto$chrom, IRanges(cyto$chromStart, cyto$chromEnd), type=cyto$name)
    # 只保留chr1:22
    cyto.mm10 <- keepSeqlevels(cyto.mm10, paste0("chr", c(1:19)), pruning.mode="coarse")
    return(cyto.mm10)
}


deal_GeneAnno <- function(genomeSession){
    # `gene-poor region`; label gene with high length > 500kb
    # knownGene: hgnc: HUGO Gene Nomenclature
    knownGene <- getTable(ucscTableQuery(genomeSession, table="knownGene", genome='hg38'))
    # hgnc <- getTable(ucscTableQuery(genomeSession, table="hgnc"))
    ## region: $transcriptClass $geneType 
    knownGene <- knownGene[knownGene$transcriptClass!='pseudo', ]
    knownGene <- knownGene[-grep('_', knownGene$chrom), ]
    knownGene.hg38 <- do.call(rbind,lapply(unique(knownGene$geneName),FUN = function(genename,knownGene){
                                data.frame(genename,
                                            chrom=paste0(unique(knownGene$chrom[knownGene$geneName==genename]),collapse=','),
                                            chromStart=min(knownGene$chromStart[knownGene$geneName==genename]),
                                            chromEnd=max(knownGene$chromEnd[knownGene$geneName==genename]),
                                            transcriptClass=paste0(unique(knownGene$transcriptClass[knownGene$geneName==genename]),collapse=','),
                                            geneType=paste0(unique(knownGene$geneType[knownGene$geneName==genename]),collapse=','))},knownGene=as.data.frame(knownGene)))
    knownGene.hg38$width <- knownGene.hg38$chromEnd - knownGene.hg38$chromStart
    knownGene.hg38$LargeGene <- ifelse(knownGene.hg38$width>500000,'Y','N')
    knownGene.hg38gr <- GRanges(knownGene.hg38$chrom, IRanges(knownGene.hg38$chromStart, knownGene.hg38$chromEnd), 
                                genename=knownGene.hg38$genename, transcriptClass=knownGene.hg38$transcriptClass, 
                                geneType=knownGene.hg38$geneType, largeGene=knownGene.hg38$LargeGene)
    # chr1:22
    knownGene.hg38gr <- keepSeqlevels(knownGene.hg38gr, paste0("chr", c(1:22)), pruning.mode="coarse")
    knownGene.hg38gr$width <- width(ranges(knownGene.hg38gr))
    return(knownGene.hg38gr)
}

deal_mouseGeneAnno <- function(genomeSession){
    # `gene-poor region`; label gene with high length > 500kb
    # knownGene: hgnc: HUGO Gene Nomenclature
    knownGene <- getTable(ucscTableQuery(genomeSession, table="knownGene", genome='mm10'))
    knownGene$name2 <- gsub('\\.\\d+', '', knownGene$name)

    library(biomaRt)
    mart <- useMart("ensembl","mmusculus_gene_ensembl")##人类选择hsapiens_gene_ensembl
    gene_name <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = knownGene$name2, mart = mart)    
    
    ### merge
    knownGene_merge <- merge(knownGene, gene_name, by.x='name2', by.y='ensembl_transcript_id', all = TRUE)
    knownGene.mm10gr <- GRanges(knownGene_merge$chrom, IRanges(knownGene_merge$txStart, knownGene_merge$txEnd), 
                                genename=knownGene_merge$external_gene_name, exonCount=knownGene_merge$exonCount)
    # 只保留chr1:19
    knownGene.mm10gr <- keepSeqlevels(knownGene.mm10gr, paste0("chr", c(1:19)), pruning.mode="coarse")
    knownGene.mm10gr$width <- width(ranges(knownGene.mm10gr))
    return(knownGene.mm10gr)
}

deal_fragile <- function(){
    # `CFS / ERFS`
    cfs_original = read.table('refData/cfs_genomic_coordinates.txt', sep='\t', head=T, fill=T)
    cfs_original.grch38 = GRanges(cfs_original)
    cfs_merge = read.table('refData/cfs_genomic_corrdinates.unique.bed', sep='\t', head=F); colnames(cfs_merge)=c('chr','start','end'); cfs_merge$chr=str_c('chr',cfs_merge$chr)
    cfs_merge.grch38 = GRanges(cfs_merge)
    findoverlap = findOverlaps(cfs_merge.grch38, cfs_original.grch38)
    findoverlapd_df = as.data.frame(findoverlap)
    findoverlapd_df$symbol = cfs_original$Symbol[subjectHits(findoverlap)]
    if(length(table(findoverlapd_df[,1]))==nrow(cfs_merge)){
        label_add = c()
        for(row in 1:nrow(cfs_merge)){
            tmp = findoverlapd_df[findoverlapd_df$queryHits==row, ]
            if(nrow(tmp)>1){
                label_add <- c(label_add, paste0(tmp$symbol, collapse=', '))
            }else{
                label_add <- c(label_add, tmp$symbol)
            }
        }
        cfs_merge$symbol = label_add
    }

    cfs.hg38 <- GRanges(cfs_merge$chr, IRanges(cfs_merge$start, cfs_merge$end), type=cfs_merge$symbol)
    # 只保留chr1:22 (invalid seqlevels: chr21, chrY)
    cfs.hg38 <- keepSeqlevels(cfs.hg38, setdiff(names(table(cfs_merge$chr)),c('chrX','chrY')), pruning.mode="coarse")
    return(cfs.hg38)
}

deal_ERFS <- function(){
    erfs_original = read.table('refData/2013_cell_ERFS_tableS1.mm8_GeneRegion.map2GRCh38location.bed', sep='\t', head=F)
    erfs_original$V1 = str_c('chr',erfs_original$V1)
    colnames(erfs_original) = c('chr', 'start', 'end', 'Symbol','no','strand')
    erfs_original.grch38 = GRanges(erfs_original)

    erfs_merge = read.table('refData/erfs2grch38.sorted.unique.bed', sep='\t', head=F); colnames(erfs_merge)=c('chr','start','end'); erfs_merge$chr=str_c('chr',erfs_merge$chr)
    erfs_merge.grch38 = GRanges(erfs_merge)
    findoverlap = findOverlaps(erfs_merge.grch38, erfs_original.grch38)
    findoverlapd_df = as.data.frame(findoverlap)
    findoverlapd_df$symbol = erfs_original$Symbol[subjectHits(findoverlap)]
    if(length(table(findoverlapd_df[,1]))==nrow(erfs_merge)){
        label_add = c()
        for(row in 1:nrow(erfs_merge)){
            tmp = findoverlapd_df[findoverlapd_df$queryHits==row, ]
            if(nrow(tmp)>1){
                label_add <- c(label_add, paste0(tmp$symbol, collapse=', '))
            }else{
                label_add <- c(label_add, tmp$symbol)
            }
        }
        erfs_merge$symbol = label_add
    }

    erfs.hg38 <- GRanges(erfs_merge$chr, IRanges(erfs_merge$start, erfs_merge$end), type=erfs_merge$symbol)
    # 只保留chr1:22 (invalid seqlevels: chr21, chrY)
    erfs.hg38 <- keepSeqlevels(erfs.hg38, paste0("chr", c(1:22)), pruning.mode="coarse")
    return(erfs.hg38)
}

add_GenomicAnnotation_userbed <- function(userbed, path, name){
    ### add to userbed:
    if(length(grep('chr',userbed$chromosome))==0){
        userbed$chromosome = str_c('chr',userbed$chromosome)
    }
    userbed = GRanges(userbed[, setdiff(colnames(userbed),'feature')])

    ### find overlaps to the cytoband
    userbed$arm <- arms$type[(findOverlaps(userbed, arms, select='first'))] 
    userbed$cyto <- cyto.hg38$type[(findOverlaps(userbed, cyto.hg38, select='first'))]
    userbed$CFS <- cfs.hg38$type[(findOverlaps(userbed, cfs.hg38, select='first'))]
    userbed$ERFS <- erfs.hg38$type[(findOverlaps(userbed, erfs.hg38, select='first'))]

    ### find overlaps to the known genes
    overlap_with_genes <- findOverlaps(userbed, knownGene.hg38)
    userbed <- as.data.frame(userbed); userbed$id=1:dim(userbed)[1]
    userbed$region <- str_c(userbed$seqnames,':',userbed$start,'-',userbed$end)

    # gene region
    overlap_with_genes_df <- as.data.frame(overlap_with_genes)
    knownGene.hg38_df <- as.data.frame(knownGene.hg38); knownGene.hg38_df$id=1:dim(knownGene.hg38_df)[1]
    knownGene.hg38_df <- knownGene.hg38_df[-grep('ENSG00', knownGene.hg38_df$genename), ]
    overlap_with_genes_df_merge1 <- merge(overlap_with_genes_df, knownGene.hg38_df, by.x='subjectHits', by.y='id')
    userbed_merge_genes <- merge(userbed,overlap_with_genes_df_merge1, by.x='id', by.y='queryHits', all.x=T)
    userbed_merge <- do.call(rbind,lapply(unique(userbed_merge_genes$region),FUN = function(region,df){
        data.frame(region=region,
        chr=unique(df[df$region==region,'seqnames.x']),start=unique(df[df$region==region,'start.x']),end=unique(df[df$region==region,'end.x']),
        arm=unique(df[df$region==region,'arm']), cyto=unique(df[df$region==region,'cyto']), CFS=unique(df[df$region==region,'CFS']), ERFS=unique(df[df$region==region,'ERFS']),
        overlappedGeneNumDensity=length(unique(df[df$region==region,'genename'])),
        genename=paste0(unique(df[df$region==region,'genename']),',',collapse=''), 
        transcriptClass=paste0(unique(unlist(strsplit(df[df$region==region,'transcriptClass'],','))),',',collapse=''), 
        geneType=paste0(unique(unlist(strsplit(df[df$region==region,'geneType'],','))),',',collapse=''), 
        largeGene=paste0(unique(df[df$region==region,'largeGene']),',',collapse='')
    )},df=userbed_merge_genes))

    userbed_merge$genename = gsub('^NA,|^NA,$|\\,$','',userbed_merge$genename)
    userbed_merge$transcriptClass = gsub('^NA,|^NA,$|\\,$','',userbed_merge$transcriptClass)
    userbed_merge$geneType = gsub('^NA,|^NA,$|\\,$','',userbed_merge$geneType)
    userbed_merge$largeGene = gsub('^NA,|^NA,$|\\,$','',userbed_merge$largeGene)

    userbed_merge[is.na(userbed_merge$CFS), 'CFS'] = ''
    userbed_merge[is.na(userbed_merge$ERFS), 'ERFS'] = ''

    write.table(userbed_merge, str_c(path,'/GenomicAnnotationIn_',name,'_addtoRegion.log'), sep='\t', row.names=F, quote=F)
    return(userbed_merge)
}

add_mm10GenomicAnnotation_userbed <- function(userbed, path, name){
    ### add to userbed:
    if(length(grep('chr',userbed$chromosome))==0){
        userbed$chromosome = str_c('chr',userbed$chromosome)
    }
    userbed = GRanges(userbed[, setdiff(colnames(userbed),'feature')])
    ### find overlaps to the cytoband
    userbed$arm <- arms$type[(findOverlaps(userbed, arms, select='first'))] #subjectHits
    userbed$cyto <- cyto.mm10$type[(findOverlaps(userbed, cyto.mm10, select='first'))]

    ### find overlaps to the known genes
    overlap_with_genes <- findOverlaps(userbed, knownGene.mm10)
    userbed <- as.data.frame(userbed); userbed$id=1:dim(userbed)[1]
    userbed$region <- str_c(userbed$seqnames,':',userbed$start,'-',userbed$end)
    # gene region
    overlap_with_genes_df <- as.data.frame(overlap_with_genes)
    knownGene.mm10_df <- as.data.frame(knownGene.mm10); knownGene.mm10_df$id=1:dim(knownGene.mm10_df)[1]
    # knownGene.mm10_df <- knownGene.mm10_df[-grep('ENSG00', knownGene.mm10_df$genename), ]
    overlap_with_genes_df_merge1 <- merge(overlap_with_genes_df, knownGene.mm10_df, by.x='subjectHits', by.y='id')
    userbed_merge_genes <- merge(userbed,overlap_with_genes_df_merge1, by.x='id', by.y='queryHits', all.x=T)
    userbed_merge_genes$start.y[is.na(userbed_merge_genes$start.y)] = 0
    userbed_merge_genes$end.y[is.na(userbed_merge_genes$end.y)] = 0
    userbed_merge <- do.call(rbind,lapply(unique(userbed_merge_genes$region),FUN = function(region,df){
        data.frame(region=region,
        chr=unique(df[df$region==region,'seqnames.x']),start=unique(df[df$region==region,'start.x']),end=unique(df[df$region==region,'end.x']),
        arm=unique(df[df$region==region,'arm']), cyto=unique(df[df$region==region,'cyto']), # CFS=unique(df[df$region==region,'CFS']), ERFS=unique(df[df$region==region,'ERFS']),
        overlappedGeneNumDensity=length(unique(df[df$region==region,'genename'])),
        genename=paste0(unique(df[df$region==region,'genename']),',',collapse='')
    )},df=userbed_merge_genes))

    userbed_merge$genename = gsub('^NA,|^NA,$|\\,$','',userbed_merge$genename)

    write.table(userbed_merge, str_c(path,'/GenomicAnnotationIn_',name,'_addtoRegion.mouse.log'), sep='\t', row.names=F, quote=F)
    # print(head(userbed_merge[order(userbed_merge$MEP, decreasing=TRUE),]))
    return(userbed_merge)
}

