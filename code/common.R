library(ggthemes, viridis)

theme_Publication <- function(base_size=10, base_family="Arial") {

      (theme_foundation(base_size=base_size, base_family=base_family) + 
        theme(
              #plot.title = element_text(face = "bold",
              #                           size = rel(1.5), 
              #                           hjust = 0.5),
               text = element_text(family = base_family, size = base_size, color = "black"),

               panel.background = element_rect(fill = "transparent", colour = NA),
               #panel.border = element_rect(colour = "black"),
               panel.border = element_rect(colour = "black", fill=NA, size=1),
               #panel.border = element_blank(),

               plot.background = element_rect(fill = "transparent", colour = NA),
               
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               
               #axis.title.y = element_text(angle=90, vjust = 2),
               #axis.title.x = element_text(vjust = -0.2),
               axis.line = element_line(colour="black"),
               axis.text = element_text(size=base_size, color="black"),
               #axis.text.x =element_text(size = rel(1), angle = 90, vjust = 0.5, hjust=1),
               #axis.text.y =element_text(size = rel(1)),
               axis.ticks = element_line(size=0.3, color="black"),
               axis.title = element_text(size = rel(1.2), color="black"),
               
               strip.background=element_blank(),
               #strip.text = element_text(face="bold"),
               strip.text.x = element_text(size=base_size, color="black"),
               strip.text.y = element_text(size=base_size, color="black", angle = 0),
               legend.position="bottom",  
               legend.background = element_rect(
                    fill = alpha("white", 0),
                    color = alpha("white", 0)
                  ),
              legend.key = element_rect(color = NA, fill = NA),
              plot.margin=unit(c(2, 2, 2, 2),"mm"),
              legend.key.size = unit(1, "line"),
              legend.text = element_text(size = base_size),

              plot.tag=element_text(size=14, face="bold")

#               legend.key = element_rect(colour = NA),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.key.size= unit(0.2, "cm"),
#               legend.margin = margin(t = 0, unit = "cm"),
#               legend.title = element_text(face="italic"),
#               plot.margin=unit(c(5, 5, 5, 5),"mm")
               
          ))      
}


get_cor <- function(x, sampleTable) {
  #x <- countai
  #sampleTable <- sampleInfo
  ha = HeatmapAnnotation(
    type = sampleTable$cells, 
    treatment = sampleTable$treatment, 
    rep = sampleTable$replicate,
    col = list(#type = c("d" = brewer.pal(3, "Set1")[1] , "n" = brewer.pal(3, "Set1")[2], "wt" = brewer.pal(3, "Set1")[3]),
               rep = c("R1"= brewer.pal(6, "Set1")[4], "R2"=brewer.pal(6, "Set1")[5])
              )#,
    #annotation_legend_param = list(
    #    rep = list(nrow = 1),
    #    treatment=list(nrow = 1),
    #    type = list(nrow = 1))
                        )
  dt <- cor(x) 

  ht <- Heatmap(dt, name="correlation", 
    heatmap_legend_param = list(direction = "horizontal"),
              column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", dt[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              top_annotation = ha)   
  return(ht)
}

get_PCA <- function(x, sampleTable, title="PCA") {
  #x <- countai 
  #sampleTable <- sampleInfo      
  pca <-  x %>% t %>% prcomp()
  summ <- summary(pca)
  p1 <- pca$x %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    merge(., sampleTable, by.x="rn", by.y="sample") %>% 
    ggplot(aes(PC1, PC2, colour=Group, shape=Replicate, label=Barcode)) +
      geom_point(size=3) +
      theme_bw() +
      ggtitle(title) +
      geom_label_repel() +
      theme(legend.position="bottom") +
      xlab(paste0("PC 1: ", round(summ$importance[2,1]*100), "%")) +
      ylab(paste0("PC 2: ", round(summ$importance[2,2]*100), "%")) 
  return(p1)
}

get_DDS <- function(txi, sampleTable, fonas, pvalue_th=0.05, LFC_th=log(2), th_counts=-1, th_groups=-1) {
  #x <- countai 
  #sampleTable <- sampleInfo   
  #background <- "control"
  dds_o <-  DESeqDataSetFromTximport(
      txi,
      sampleTable,      
      design = ~Group
      )
  dds_o$Group <- relevel(dds_o$Group, ref = fonas)  
  #optional low count filtering 
  keep <- rowSums(counts(dds_o) >= th_counts) >= th_groups
  dds_i <- dds_o[keep,]
  dds <- DESeq(dds_o, quiet=TRUE)
  coef_name <- resultsNames(dds)[grepl("Group_", resultsNames(dds))]
  res <- results(dds, 
                 alpha=pvalue_th, 
                 lfcThreshold=LFC_th)
  vsd <- vst(dds, blind=FALSE)
  vsd_blind <- vst(dds, blind=TRUE)
  dds_shrunken <- lfcShrink(dds, 
                            coef = coef_name, 
                            type = "apeglm")
  return(list(dds=dds, 
              dds_o=dds_o, 
              res=res, 
              vsd=vsd, 
              vsd_blind=vsd_blind, 
              dds_shrunken=dds_shrunken))
}


get_heatmap <- function(x, countai, titlas) {
  #x <- dge_id
  #countai <- vst_counts
  tmp <- countai %>% as.data.table(., keep.rownames=TRUE) %>%  .[rn %in% x, ] %>% .[, rn := NULL]
  ht <- tmp %>% 
    t() %>% 
    scale() %>% 
    t() %>% 
    Heatmap(, 
      name="z-score", 
      column_title=titlas)
  return(ht)      
}

get_volcano <- function(x, pthres, lfc){
  #x <- dds_list[["res"]]
  tmp <- x %>% 
    as.data.table() %>% 
    .[, type := "NotSig"] %>% 
    .[padj <= pthres, type := "p. adjusted"] %>% 
    .[abs(log2FoldChange) >= lfc, type := "LFC"] %>% 
    .[abs(log2FoldChange) >= lfc & padj <= pthres, type := "LFC & p. adjusted"] 

  p <- ggplot(tmp, aes(log2FoldChange, y=(-1)*log10(padj), colour=type)) +  
    geom_point(size=0.5) +
    theme_bw()
  return(p)
}

replace_slashes <- function(s) {
  parts <- unlist(strsplit(s, "/"))
  new_string <- parts[1]
  for (i in 2:length(parts)) {
    # Use "/" every 10th separator, otherwise ";"
    sep <- ifelse((i - 1) %% 10 == 0, "/", ";")
    new_string <- paste0(new_string, sep, parts[i])
  }
  return(new_string)
}



make_table_genes <- function(egoX, captionas) {
  tmp <- sapply(egoX$geneID, replace_slashes)
  egoX <- as.data.table(egoX)
  egoX$geneID_modified <- tmp
  egoX$geneID <- NULL

  DT::datatable(
    as.data.table(egoX),
    options = list(
      pageLength = 10,
      autoWidth = TRUE
    ),
    filter = 'top',  # Add filter boxes on top of each column
    caption = captionas,
    rownames = FALSE
)
}

make_dot_plot <- function(x, captionas) {
  dotplot(x,
              font.size = 13,
              label_format = 30,
              showCategory  = 8) + 
  ggtitle(captionas) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")
}

get_DDS_featureCounts <- function(txi, sampleTable, fonas, ilgiai, pvalue_th=0.05, LFC_th=log(2)) {
  #x <- countai 
  #sampleTable <- sampleInfo   
  #background <- "control"
  m <- txi[, -(1:2)] %>% 
    as.matrix()
  rownames(m) <- txi$Geneid  
  dds_o <-  DESeqDataSetFromMatrix(
      m,
      sampleTable,      
      design = ~Group
      )
  mcols(dds_o)$basepairs <- ilgiai
  dds_o$Group <- relevel(dds_o$Group, ref = fonas)  
  coef_name <- resultsNames(dds_o)[grepl("Group_", resultsNames(dds_o))]
  dds <- DESeq(dds_o, quiet=TRUE)
  res <- results(dds_o, alpha=pvalue_th, lfcThreshold=LFC_th)
  vsd <- vst(dds, blind=FALSE)
  vsd_blind <- vst(dds, blind=TRUE)
  dds_shrunken <- lfcShrink(dds, 
                          coef = coef_name, 
                          type = "apeglm")
  return(list(dds=dds, 
              dds_o=dds_o, 
              res=res, 
              vsd=vsd, 
              vsd_blind=vsd_blind, 
              dds_shrunken=dds_shrunken))
}

do_gsea_ridge_LFC <- function(dds_list) {

  res <- dds_list[["dds_shrunken"]][order(dds_list[["dds_shrunken"]]$log2FoldChange, decreasing = TRUE), ] 

  # extract ranking vector
  geneList <- res$log2FoldChange
  names(geneList) <- rownames(res)   # these must be gene IDs (ENSEMBL or SYMBOL)
  geneList <- geneList[!is.na(geneList)]            # remove NAs
  names(geneList) <- gsub("\\.[0-9]+", "", names(geneList))

  gsea_go <- gseGO(geneList,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENSEMBL", 
                  pAdjustMethod = "BH")
 total_ridge <- ridgeplot(gsea_go, showCategory = 20) + theme(axis.text.y = element_text(size = 8))
  g <- gsea_go@result
  g_pos <- g[g$NES > 0, ]
  g_neg <- g[g$NES < 0, ]

  gsea_pos <- gsea_go
  gsea_pos@result <- g_pos

  gsea_neg <- gsea_go
  gsea_neg@result <- g_neg

  ridge_pos <- ridgeplot(gsea_pos, showCategory = 15) +
    ggtitle("Positively Enriched Pathways (NES > 0)") +
    theme(legend.position="bottom")

  ridge_neg <- ridgeplot(gsea_neg, showCategory = 15) +
    ggtitle("Positively Enriched Pathways (NES > 0)") +
    theme(legend.position="bottom")
 
  network <- enrichmentNetwork(gsea_go@result, drawEllipses = TRUE, fontSize = 2.5)
  
  return(list(total_ridge=total_ridge, 
         ridge_pos=ridge_pos,
         ridge_neg=ridge_neg,
         tinklas=network,
         results=gsea_go))
}

do_gsea_ridge_stat <- function(dds_list) {

  res <- dds_list[["res"]][order(dds_list[["res"]]$stat, decreasing = TRUE), ] 

  # extract ranking vector
  geneList <- res$stat
  names(geneList) <- rownames(res)   # these must be gene IDs (ENSEMBL or SYMBOL)
  geneList <- geneList[!is.na(geneList)]            # remove NAs
  names(geneList) <- gsub("\\.[0-9]+", "", names(geneList))

  gsea_go <- gseGO(geneList,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENSEMBL", 
                  pAdjustMethod = "BH")
 total_ridge <- ridgeplot(gsea_go, showCategory = 20) + theme(axis.text.y = element_text(size = 8))
  g <- gsea_go@result
  g_pos <- g[g$NES > 0, ]
  g_neg <- g[g$NES < 0, ]

  gsea_pos <- gsea_go
  gsea_pos@result <- g_pos

  gsea_neg <- gsea_go
  gsea_neg@result <- g_neg

  ridge_pos <- ridgeplot(gsea_pos, showCategory = 15) +
    ggtitle("Positively Enriched Pathways (NES > 0)") +
    theme(legend.position="bottom")

  ridge_neg <- ridgeplot(gsea_neg, showCategory = 15) +
    ggtitle("Positively Enriched Pathways (NES > 0)") +
    theme(legend.position="bottom")
 
  network <- enrichmentNetwork(gsea_go@result, drawEllipses = TRUE, fontSize = 2.5)
  
  return(list(total_ridge=total_ridge, 
         ridge_pos=ridge_pos,
         ridge_neg=ridge_neg,
         tinklas=network,
         results=gsea_go))
}