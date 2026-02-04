library(ggthemes, viridis, igraph)

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
         geom_text_repel(size = 6,
                         box.padding = 1.0, 
                         point.padding = 0.8,
                         max.overlaps = Inf, 
                         show.legend = FALSE) +
         theme_minimal() +
         labs(x = xlab, y = ylab) +
         scale_color_brewer(palette = "Set1")+
         theme(legend.position = "bottom",
         legend.text = element_text(size = 14),
         legend.title = element_text(size = 16, face = "bold"),
         plot.title = element_text(size = 20, face = "bold"),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18))+
        coord_fixed()
}




plot_net <- function(x){
modules <- split(x$Gene, x$Module)
modules <- modules[sapply(modules, length) >= 5]

module_expr <- lapply(modules, function(genes) {
  genes <- intersect(genes, rownames(mat_filt_pc_sym))
  
  if (length(genes) < 2) return(NULL)
  
  colMeans(mat_filt_pc_sym[genes, , drop = FALSE])
})

module_expr <- module_expr[!sapply(module_expr, is.null)]
module_expr <- do.call(rbind, module_expr)
rownames(module_expr) <- paste0("Module_", rownames(module_expr))
final_hits$Module <- paste0("Module_", final_hits$Module)

edges <- final_hits[abs(final_hits$rho) >= 0.9 & final_hits$perm_pval <= 0.01,]
lncGR$gene_id <- sub("\\..*$", "", lncGR$gene_id)
lnc_map <- setNames(
  lncGR$gene_name,
  lncGR$gene_id)
edges$lncRNA_symbol <- lnc_map[edges$lncRNA]


edge_df <- data.frame(
  from   = edges$lncRNA_symbol,
  to     = edges$Module,
  weight = abs(edges$rho),
  sign   = ifelse(edges$rho > 0, "positive", "negative"))

g_lnc_mod <- graph_from_data_frame(
  edge_df,
  directed = FALSE)

V(g_lnc_mod)$type <- ifelse(
  grepl("^Module_", V(g_lnc_mod)$name),
  "Module",
  "lncRNA")
V(g_lnc_mod)$degree <- degree(g_lnc_mod)

V(g_lnc_mod)$size <- ifelse(
  V(g_lnc_mod)$type == "lncRNA",
  6 + V(g_lnc_mod)$degree * 2,
  12)

V(g_lnc_mod)$color <- ifelse(
  V(g_lnc_mod)$type == "lncRNA",
  "pink",
  "darkcyan")

E(g_lnc_mod)$color <- ifelse(
  edge_df$sign == "positive",
  "lightgreen",
  "midnightblue"
)

E(g_lnc_mod)$width <- edge_df$weight * 3

plot(
  g_lnc_mod,
  layout = layout_with_fr,
  vertex.label.cex = 0.7,
  vertex.frame.color = "white",
  main = "lncRNA–PPI Module Network"
)
}