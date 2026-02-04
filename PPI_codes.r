

# PPI Network and identification of hub lncRNAs
```{r}
library(igraph)
library(ggraph)

cor_df_filtered <- cor_df_sig %>%
  filter(value > 0.90)

# Create network
network <- graph_from_data_frame(d = cor_df_sig, directed = FALSE)

# Add attributes
V(network)$type <- ifelse(V(network)$name %in% cor_df_sig$lncRNA, "lncRNA", "mRNA")
E(network)$color <- ifelse(E(network)$value > 0, "red", "blue")
E(network)$width <- abs(E(network)$value) * 2

# Plot
ggraph(network, layout = "fr") +
  geom_edge_link(aes(color = factor(E(network)$color), width = E(network)$width)) +
  scale_edge_color_manual(values = c("red", "blue")) +
  geom_node_point(aes(color = type), size = 3) +
  scale_color_manual(values = c("lncRNA" = "orange", "mRNA" = "green")) +
  theme_graph() +
  labs(title = "lncRNA-mRNA Correlation Network")
```




```{r}
geneGR_df <- as.data.frame(geneGR)
geneGR_df <- geneGR_df %>%
  distinct(gene_id, .keep_all = TRUE)
sig_with_symbol <- sig %>%
  left_join(geneGR_df %>% 
  dplyr::select(gene_id, gene_name), by = c("rn" = "gene_id"))
sig_with_symbol$rn <- sub("\\..*$","", sig_with_symbol$rn)
#sig_with_symbol <- sig_with_symbol[sig_with_symbol$rn %in% genes,]
```
```{r, eval = F, echo = F}
library(STRINGdb)
string_db <- STRINGdb$new(version = "12.0",
                          species = 9606,
                          score_threshold = 400,
                          network_type = "full")
mapped <- string_db$map(data.frame(gene_name = sig_with_symbol$gene_name), "gene_name", removeUnmappedRows = TRUE)
#mapped <- mapped[order(abs(mapped$log2FoldChange), decreasing = TRUE), ]
mapped_valid <- mapped[!is.na(mapped$STRING_id), ]

ppi_edges <- string_db$get_interactions(mapped_valid$STRING_id)
ppi_edges <- ppi_edges[ppi_edges$combined_score >= 700, ]
id2gene <- mapped_valid$gene_name
names(id2gene) <- mapped_valid$STRING_id

ppi_edges$gene1 <- id2gene[ppi_edges$from]
ppi_edges$gene2 <- id2gene[ppi_edges$to]

ppi_edges <- ppi_edges[!is.na(ppi_edges$gene1) & !is.na(ppi_edges$gene2), ]
library(igraph)

g <- graph_from_data_frame(
  ppi_edges[, c("gene1", "gene2")],
  directed = FALSE
)

g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
cl <- cluster_louvain(g)
modules <- split(V(g)$name, membership(cl))
modules <- modules[sapply(modules, length) >= 3]
ppi_modules <- data.frame(
  Module = rep(names(modules), sapply(modules, length)),
  Gene   = unlist(modules)
)

write.csv(ppi_modules, "~/NB_lncRNA/output/code/BioInsights/PPI_modules.csv", row.names = FALSE)
```

```{r, eval = F, echo = F}
ppi <- read.csv("~/NB_lncRNA/output/code/BioInsights/PPI_modules.csv")


modules <- split(ppi$Gene, ppi$Module)

ens2sym <- sig_with_symbol$gene_name
names(ens2sym) <- sig_with_symbol$rn
new_names <- ens2sym[rownames(mat_filt_pc)]
keep <- !is.na(new_names)

mat_filt_pc_sym <- mat_filt_pc[keep, ]
rownames(mat_filt_pc_sym) <- new_names[keep]


module_expr <- lapply(modules, function(genes) {
  genes <- intersect(genes, rownames(mat_filt_pc_sym))
  colMeans(mat_filt_pc_sym[genes, , drop = FALSE])
})

module_expr <- do.call(rbind, module_expr)

results <- data.frame()

for (lnc in rownames(mat_filt_lnc)) {
  for (mod in rownames(module_expr)) {

    cor_test <- cor.test(
      as.numeric(mat_filt_lnc[lnc, ]),
      as.numeric(module_expr[mod, ]),
      method = "spearman"
    )

    results <- rbind(results, data.frame(
      lncRNA = lnc,
      Module = mod,
      rho = cor_test$estimate,
      pval = cor_test$p.value
    ))
  }
}

perm_test <- function(x, y, nperm = 500) {

  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 5) return(NA_real_)

  # FORCE numeric scalar
  obs <- as.numeric(
    suppressWarnings(cor(x, y, method = "spearman"))
  )

  perm_rho <- numeric(nperm)

  for (i in seq_len(nperm)) {
    perm_rho[i] <- as.numeric(
      suppressWarnings(cor(sample(x), y, method = "spearman"))
    )
  }

  mean(abs(perm_rho) >= abs(obs))
}

results_prefilt <- subset(results, abs(rho) >= 0.6)
nrow(results_prefilt)

results_prefilt$perm_pval <- mapply(
  function(l, m) {
    perm_test(
      as.numeric(mat_filt_lnc[l, ]),
      as.numeric(module_expr[m, ])
    )
  },
  results_prefilt$lncRNA,
  results_prefilt$Module
)

final_hits <- subset(
  results_prefilt,
  abs(rho) >= 0.7 & perm_pval < 0.05
)

final_hits <- final_hits[order(-abs(final_hits$rho)), ]
saveRDS(final_hits, "~/NB_lncRNA/output/code/BioInsights/final_hits.RDS")
```