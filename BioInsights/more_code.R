library(ggthemes, viridis)
library(igraph)


get_pca <- function(matrica, info){
  # matrica <- filtruoti genai (yra pokytis)
  # info <- SampleInfo, meginiu informacija su tiriama grupe ir meginio pavadinimu
pca_result <- prcomp(t(matrica), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Barcode <- rownames(pca_df)
pca_df <- merge(pca_df, info, by = "Barcode")

var_exp <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
xlab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
ylab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Barcode)) +
         geom_point(size = 4) +
         geom_text_repel(size = 4,
                         box.padding = 1.0, 
                         point.padding = 0.8,
                         max.overlaps = Inf, 
                         show.legend = FALSE) +
         theme_minimal() +
         labs(x = xlab, y = ylab) +
         scale_color_brewer(palette = "Set1")+
         theme(legend.position = "bottom",
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 14, face = "bold"),
         plot.title = element_text(size = 18, face = "bold"),
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         axis.text.x = element_text(size = 16),
         axis.text.y = element_text(size = 16))
}

getpca <- function(matrica, info){
  # matrica <- filtruoti genai (yra pokytis)
  # info <- SampleInfo, meginiu informacija su tiriama grupe ir meginio pavadinimu
pca_result <- prcomp(t(matrica))
pca_df <- as.data.frame(pca_result$x)
pca_df$Barcode <- rownames(pca_df)
pca_df <- merge(pca_df, info, by = "Barcode")

var_exp <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
xlab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
ylab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Barcode)) +
         geom_point(size = 4) +
         geom_text_repel(size = 4,
                         box.padding = 1.0, 
                         point.padding = 0.8,
                         max.overlaps = Inf, 
                         show.legend = FALSE) +
         theme_minimal() +
         labs(x = xlab, y = ylab) +
         scale_color_brewer(palette = "Set1")+
         theme(legend.position = "bottom",
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 14, face = "bold"),
         plot.title = element_text(size = 18, face = "bold"),
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         axis.text.x = element_text(size = 16),
         axis.text.y = element_text(size = 16))
}


get_cor <- function(x, title= "PCA"){
# x <- counts
cor_mat <- cor(x)
cor_df <- reshape2::melt(cor_mat)
colnames(cor_df) <- c("Var1", "Var2", "value")

ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white") +
        labs(title="B", fill= " ")+
        scale_fill_gradientn(
        colors = c("#2626b6", "white", "#c52828"),
        limits = c(min(cor_df$value), max(cor_df$value))) +
        ggtitle(title)+
        theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16, face = "bold")) +
        coord_fixed()
}


get_near_table <- function(lncGR, geneGR){
  # GRanges objektai
  ## lncGR - atfiltruoti lncRNA
  ## geneGR - visi baltyma koduojantys genai (suzinoti visus kaimynus)
nearest_idx <- nearest(lncGR, geneGR)
nearest_gene <- geneGR$gene_name[nearest_idx]
distance_to_pc <- distance(lncGR, geneGR[nearest_idx])
genomic_context <- rep("intergenic", length(lncGR))

gene_upstream <- promoters(geneGR, upstream=1500, downstream=0)
hits <- findOverlaps(lncGR, gene_upstream, ignore.strand=FALSE)
qh <- queryHits(hits)
sh <- subjectHits(hits)

lnc_strand <- as.character(strand(lncGR[qh]))
gene_strand <- as.character(strand(geneGR[sh]))
lnc_tss <- start(lncGR[qh])
gene_tss <- ifelse(gene_strand == "+", start(geneGR[sh]), end(geneGR[sh]))
bidir_hit <- (gene_strand == "+" & lnc_tss < gene_tss & (gene_tss - lnc_tss) <= 2000) |
                        (gene_strand == "-" & lnc_tss > gene_tss & (lnc_tss - gene_tss) <= 2000)
bidir <- unique(qh[bidir_hit])
genomic_context[bidir] <- "bidirectional"


hits <- findOverlaps(lncGR, geneGR, ignore.strand=TRUE)
qh <- queryHits(hits)
sh <- subjectHits(hits)
lnc_strand <- as.character(strand(lncGR[qh]))
gene_strand <- as.character(strand(geneGR[sh]))

close_to_gene <- distance_to_pc[qh] <= 10000
is_antisense_hit <- (lnc_strand != gene_strand) & close_to_gene
antisense_idx <- setdiff(unique(qh[is_antisense_hit]), bidir)
genomic_context[antisense_idx] <- "antisense"

is_sense_hit <- lnc_strand == gene_strand
sense_idx <- setdiff(unique(qh[is_sense_hit]), c(bidir, antisense_idx))
genomic_context[sense_idx] <- "sense_overlapping"

final_table <- data.frame(
  lncRNA = lncGR$gene_name,
  genomic_context = genomic_context,
  nearest_protein_coding_gene = nearest_gene,
  distance = distance_to_pc)
final_table <- final_table[!duplicated(final_table[c("lncRNA", "nearest_protein_coding_gene", "distance")]), ]
final_table <- final_table %>%
  group_by(lncRNA, nearest_protein_coding_gene) %>%
  slice_min(distance) %>%
  ungroup()
}

plot_heatmap <- function(mat, title_name, annotation_col, ann_colors) {
  pheatmap(mat,
           scale = "row",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = title_name,
           silent = TRUE)
}

get_plot <- function(x) {
 ggplot(x, aes(x = genomic_context, y = count, fill = genomic_context)) +
  geom_col() +
  geom_label(aes(label = count),
             fill = "white",
             color = "black",
             size = 4,
             label.size = 0.7,
             show.legend = FALSE) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Genomic context", y = "Count", title= "Genominis kontekstas DE lncRNA") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.position = "none",
        plot.title = element_text(size=15))
}

cor_heatmap <- function(x, title_name) {
  pheatmap(x,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = FALSE,
         show_colnames = FALSE,
        main = title_name)
}

tinklas <- function(x) {
  # x <- tg
  ggraph(x, layout = "fr") +
  geom_edge_link(aes(color = edge_type), alpha = 0.3) +
  geom_node_point(aes(color = type, size = degree, shape = type)) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
  scale_color_manual(values = c("lncRNA" = "orange", "protein" = "skyblue")) +
  scale_edge_color_manual(values = c("lncRNA_corr" = "red", "PPI" = "#261ace")) +
  scale_shape_manual(values = c("lncRNA" = 17, "protein" = 16)) +
  theme_void() +
  ggtitle("lncRNA-PPI Network with Modules and Hubs")
}


TOP_tinklas <- function(x, top_names) {
  x <- x %>%
  activate(nodes) %>%
  mutate(type2 = case_when(
    name %in% top_names ~ "top_lncRNA",
    type == "lncRNA" ~ "lncRNA",
    TRUE ~ "protein"
  ))

  ggraph(x, layout = "fr") +
    geom_edge_link(aes(color = edge_type), alpha = 0.3) +
    geom_node_point(aes(color = type2, shape = type)) +
    geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
    scale_color_manual(values = c(
  "top_lncRNA" = "#ff2579",
  "lncRNA" = "orange",
  "protein" = "skyblue"
)) +
    scale_edge_color_manual(values = c("lncRNA_corr" = "#ff7f7f", "PPI" = "blue")) +
    scale_shape_manual(values = c("lncRNA" = 17, "protein" = 16)) +
    theme_void()
}

topVarGenes <- function(vsd, n=500) {
  mat <- assay(vsd)
  rv <- apply(mat, 1, var)
  mat[order(rv, decreasing=TRUE)[1:n], ]    }

get_dist <- function(counts, clone){
  dist_mat <- as.matrix(dist(t(counts)))

  idx_T <- which(clone=="TC1")
  idx_K <- which(clone=="K7")

  within_T <- dist_mat[idx_T, idx_T][upper.tri(dist_mat[idx_T, idx_T])]
  within_K <- dist_mat[idx_K, idx_K][upper.tri(dist_mat[idx_K, idx_K])]

  between <- dist_mat[idx_T, idx_K]
  
  mean_within <- mean(c(within_T, within_K))
  mean_between <- mean(between)
  
  ratio <- mean_between / mean_within
  return(list(ratio=ratio, dist_mat=dist_mat))  }


heatmap_ident <- function(mat, genes, gene_set, ann_colors, annotation_col){
  # mat -> mat_pc
  # genes -> all_genes
expr_subset <- mat[rownames(mat) %in% genes, ]
expr_subset <- expr_subset[apply(expr_subset, 1, function(x) {all(is.finite(x)) && sd(x) != 0}),]

gene_annotation <- stack(gene_set)
colnames(gene_annotation) <- c("Gene", "Program")

gene_annotation <- gene_annotation[gene_annotation$Gene %in% rownames(expr_subset), ]

gene_annotation <- gene_annotation[!duplicated(gene_annotation$Gene), ]
rownames(gene_annotation) <- gene_annotation$Gene
gene_annotation <- gene_annotation["Program"]

expr_subset_ordered <- expr_subset[rownames(gene_annotation), ] 

pheatmap(expr_subset_ordered,
         scale = "row",
         show_rownames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_row = gene_annotation,
         clustering_method = "complete",
         cluster_cols = TRUE,
         cluster_rows = FALSE)

}
