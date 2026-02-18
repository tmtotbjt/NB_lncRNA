library(ggthemes, viridis)
library(igraph)


get_long <- function(z, x) {z %>%
  filter(!is.na(.data[[x]])) %>%
  separate_rows(.data[[x]], sep = ",") %>%
  mutate(across(all_of(x), trimws))
}

get_fc <- function(a, b) {
  a %>%
  mutate(direction = case_when(
    log2FoldChange > 1  ~ "up",
    log2FoldChange < -1 ~ "down",
    TRUE ~ "ns")) %>%
  filter(direction != "ns") %>%
  group_by(.data[[b]]) %>%
  summarise(n = n_distinct(gene_id), 
            up = n_distinct(gene_id[direction == "up"]),
            down = n_distinct(gene_id[direction == "down"]),
            mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n))
}
get_dir <- function(a, b) {
  get_fc(a, b) %>%
    dplyr::select(.data[[b]], up, down) %>%
    pivot_longer(
      cols = c(up, down),
      names_to = "direction",
      values_to = "count"
    )
}

get_plot <- function(z, x) {ggplot(z,
       aes(x = reorder(.data[[x]], count), y = count, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size=15, hjust = 0))
}


get_prop <- function(sig, all, col){
    sig_prop <- table(sig[[col]]) / nrow(sig)
    bg_prop  <- table(all[[col]]) / nrow(all)

    prop_df <- merge(
    data.frame(temp = names(sig_prop), Sig = as.numeric(sig_prop)),
    data.frame(temp = names(bg_prop),  Bg  = as.numeric(bg_prop)),
    by = "temp",
    all = TRUE)

  names(prop_df)[names(prop_df) == "temp"] <- col

    prop_df$ratio <- prop_df$Sig / prop_df$Bg

  prop_df %>%
    filter(!is.na(Sig), Sig > 0) %>%
    arrange(desc(ratio))

}

get_pca <- function(matrica, info){
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
         axis.text.y = element_text(size = 16))+
        coord_fixed()
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