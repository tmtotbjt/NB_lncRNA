
tmp <-  readRDS("~/NB_lncRNA/output/code/lncRNA/DataQC/data2DEG.RDS")
sampleInfo <- tmp[[2]]
samples2do <- c("T1cH1", "T1cH2", "T1cN1", "T1cN2", "K7cH1", "K7cH2", "K7cN1", "K7cN2")
sampleInfo <- sampleInfo[Barcode %in% samples2do, ] 
sampleInfo$Group <- sampleInfo$Oxygen

Status <- factor(sampleInfo$Group, levels = c("Normox", "Hypox"))
design <- model.matrix(~Status)

#lncRNAs
list_of_4 <- tmp[[1]]
counts_lnc <- list_of_4[[2]]
counts_lnc <- counts_lnc[, samples2do]

dds <- readRDS("~/NB_lncRNA/output/code/BioInsights/Oxygen/dds_lnc.RDS")
vsd <- assay(vst(dds)) %>% as.data.table(., keep.rownames=TRUE) 

mat <- vsd[,-1]
mat <- as.matrix(mat)
rownames(mat) <- vsd$rn

res <- results(dds) %>% as.data.table(., keep.rownames=TRUE) 

sig_lncRNAs <- subset(res, padj < 0.05)
sig_lncRNAs <- sig_lncRNAs[order(sig_lncRNAs$padj), ]
sig_lncRNAs <- as.data.frame(sig_lncRNAs)
mat_filt_lnc <- mat[rownames(mat) %in% sig_lncRNAs$rn, ]
rownames(mat_filt_lnc) <- sub("\\..*$", "", rownames(mat_filt_lnc))
fit <- lmFit(mat_filt_lnc, design)
mat_lnc <- residuals.MArrayLM(fit, mat_filt_lnc)

#mRNAs
dds <- readRDS("~/NB_lncRNA/output/code/BioInsights/Oxygen/dds_pc.RDS")
vsd_pc <- vst(dds)
mat_pc <- assay(vsd_pc)
res_pc <- results(dds) %>% as.data.table(., keep.rownames=TRUE)

sig <- subset(res_pc, padj < 0.05)
sig <- sig[order(sig$padj), ]
sig <- as.data.frame(sig)
mat_filt_pc <- mat_pc[rownames(mat_pc) %in% sig$rn, ]
rownames(mat_filt_pc) <- sub("\\..*$", "", rownames(mat_filt_pc))
mat_lnc <- residuals.MArrayLM(fit, mat_filt_lnc)
fit <- lmFit(mat_filt_pc, design)
mat_pc <- residuals.MArrayLM(fit, mat_filt_pc)


# Cor
cor_matrix <- cor(t(mat_lnc), t(mat_pc), method = "spearman")
cor_df <- reshape2::melt(cor_matrix)
colnames(cor_df) <- c("lncRNA", "mRNA", "value")
cor_df_sig <- subset(cor_df, abs(value) > 0.75)

# GenomicRanges
lncGR <- readRDS("~/NB_lncRNA/output/code/BioInsights/Oxygen/lncGR_padj.RDS")
geneGR <- readRDS("~/NB_lncRNA/output/code/BioInsights/Oxygen/geneGR_padj.RDS")


# Traukiam
library(GenomicRanges)

lnc_ext <-  promoters(lncGR, upstream = 50000, downstream = 50000)
hits <- findOverlaps(lnc_ext, geneGR)

pairs <- data.frame(
  lncRNA = lncGR$gene_id[queryHits(hits)],
  mRNA = geneGR$gene_id[subjectHits(hits)])
pairs_unique <- unique(pairs)

pairs$lncRNA <- sub("\\..*$", "", pairs$lncRNA)
pairs$mRNA <- sub("\\..*$", "", pairs$mRNA)

cis_cor_df <- merge(cor_df_sig, pairs, by = c("lncRNA", "mRNA"))
cis_cor_df <- unique(cis_cor_df)
write.csv(cis_cor_df, "~/NB_lncRNA/output/code/cor_gc/HN_all.csv")



# ENS -> Symbol
geneGR_df <- as.data.frame(geneGR)
geneGR_df <- geneGR_df %>%
  distinct(gene_id, .keep_all = TRUE)
geneGR_df$gene_id <- sub("\\..*$", "", geneGR_df$gene_id)


sig_links <- merge(cis_cor_df,
                   geneGR_df[, c("gene_id", "gene_name")],
                   by.x = "mRNA",
                   by.y = "gene_id",
                   all.x = TRUE)
colnames(sig_links)[colnames(sig_links) == "gene_name"] <- "mRNA"

lncGR_df <- as.data.frame(lncGR)
lncGR_df <- lncGR_df %>%
  distinct(gene_id, .keep_all = TRUE)
lncGR_df$gene_id <- sub("\\..*$", "", lncGR_df$gene_id)

sig_links <- merge(sig_links,
                   lncGR_df[, c("gene_id", "gene_name")],
                   by.x = "lncRNA",
                   by.y = "gene_id",
                   all.x = TRUE)
colnames(sig_links)[colnames(sig_links) == "gene_name"] <- "lncRNA"

sig_links <- sig_links[,c(5, 4, 3)]
write.csv(sig_links, "~/NB_lncRNA/output/code/cor_gc/sym_HN.csv")
