

library(pacman)
p_load(readxl, dplyr, tidygraph, rtracklayer)


markers_from_articles <- read_excel("~/NB_lncRNA/input/markers_from_articles.xlsx")

mes <-  markers_from_articles[,1]
adrn <- markers_from_articles[,2]
NCC <-  markers_from_articles[,3]
NonA <- markers_from_articles[,4]
crmff <-  markers_from_articles[,5]
symblast <- markers_from_articles[,6]
SCPs <-  markers_from_articles[,7]
OSL <- markers_from_articles[,8]

adrn <- adrn %>% drop_na()
NCC <- NCC %>% drop_na()
NonA <- NonA %>% drop_na()
crmff <- crmff %>% drop_na()
symblast <- symblast %>% drop_na()
SCPs <- SCPs %>% drop_na()
OSL <- OSL %>% drop_na()


markers_from_articles <- read_xlsx("~/NB_lncRNA/input/cell_ident/Neuronal markers.xlsx", sheet = 2)
SCPs_k <-  markers_from_articles[,1]
Bridge <- markers_from_articles[,2]
SymChrmff <-  markers_from_articles[,3]
chrmff <- markers_from_articles[,4]
mes_k <-  markers_from_articles[,5]
endothel <- markers_from_articles[,6]
erytoblasts <-  markers_from_articles[,7]
cortex <- markers_from_articles[,8]
leuko <- markers_from_articles[,9]
some_cells <-  markers_from_articles[,10]
NB <- markers_from_articles[,11]




gtf <- import("~/NB_lncRNA/input/gencode.v49.annotation.gtf")
#gtf_lnc <- gtf[gtf$gene_type == "lncRNA"]
#gtf_pc <- gtf[gtf$gene_type == "protein_coding"]

map <- data.frame(
  gene_id = gtf$gene_id,
  gene_name = gtf$gene_name)
  
saveRDS(map, "~/NB_lncRNA/input/map_gene.RDS")

geneGR <- GRanges(
    seqnames = seqnames(gtf),
    ranges   = ranges(gtf),
    strand   = strand(gtf),
    gene_name = mcols(gtf)$gene_name,
    gene_id   = mcols(gtf)$gene_id)


SCPs_k <- geneGR[geneGR$gene_id %in% SCPs_k$SCPs , ]
Bridge <- geneGR[geneGR$gene_name %in%  Bridge$Bridge, ]
SymChrmff <- geneGR[geneGR$gene_name %in% SymChrmff$'Sympathoblastic/Chromaffin', ]
chrmff <- geneGR[geneGR$gene_name %in% chrmff$Chromaffin , ]
mes_k <- geneGR[geneGR$gene_name %in% mes_k$Mesenchyme , ]
endothel <- geneGR[geneGR$gene_name %in% endothel$Endothelium , ]
erytoblasts <- geneGR[geneGR$gene_name %in% erytoblasts$Erythoblasts , ]
cortex <-  geneGR[geneGR$gene_name %in% cortex$Cortex , ]
leuko <-  geneGR[geneGR$gene_name %in% leuko$Leukocytes , ]
some_cells <-   geneGR[geneGR$gene_name %in% some_cells$Other , ]
NB <- geneGR[geneGR$gene_name %in% NB$Neuroblastoma , ]

SCPs_k_ens <- unique(SCPs_k$gene_id)
Bridge_ens <- unique(Bridge$gene_id)
SymChrmff_ens <- unique(SymChrmff$gene_id)
chrmff_ens <- unique(chrmff$gene_id)
mes_k_ens <- unique(mes_k$gene_id)
endothel_ens <- unique(endothel$gene_id)
erytoblasts_ens <- unique(erytoblasts$gene_id)
cortex_ens <- unique(cortex$gene_id)
leuko_ens <- unique(leuko$gene_id)
some_cells_ens <- unique(some_cells$gene_id)
NB_ens <- unique(NB$gene_id)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SCPs")
writeData(wb, "SCPs", data.frame(gene_id = SCPs_k_ens))

addWorksheet(wb, "Bridge")
writeData(wb, "Bridge", data.frame(gene_id = Bridge_ens))

addWorksheet(wb, "Sympathoblastic/Chromaffin")
writeData(wb, "Sympathoblastic/Chromaffin", data.frame(gene_id = SymChrmff_ens))

addWorksheet(wb, "Chromaffin")
writeData(wb, "Chromaffin", data.frame(gene_id = chrmff_ens))

addWorksheet(wb, "Mesenchyme")
writeData(wb, "Mesenchyme", data.frame(gene_id = mes_k_ens))

addWorksheet(wb, "Endothelium")
writeData(wb, "Endothelium", data.frame(gene_id = endothel_ens))

addWorksheet(wb, "Erythoblasts")
writeData(wb, "Erythoblasts", data.frame(gene_id = erytoblasts_ens))

addWorksheet(wb, "Cortex")
writeData(wb, "Cortex", data.frame(gene_id = cortex_ens))

addWorksheet(wb, "Leukocytes")
writeData(wb, "Leukocytes", data.frame(gene_id = leuko_ens))

addWorksheet(wb, "Other")
writeData(wb, "Other", data.frame(gene_id = some_cells_ens))

addWorksheet(wb, "NB")
writeData(wb, "NB", data.frame(gene_id = NB_ens))

saveWorkbook(wb, "~/NB_lncRNA/input/cell_ident/gene_sets_lith.xlsx", overwrite = TRUE)


#geneGR$gene_id <- sub("\\..*$", "", geneGR$gene_id)
geneGR[geneGR$gene_id %in% "ENSG00000241860",]

