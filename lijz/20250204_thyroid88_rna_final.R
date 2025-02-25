packages <- c("tidyr", "dplyr", "plotly", "readr", "writexl", "readxl",
              "here", "pvclust", "DESeq2", "stringr", "ggplot2", "sesame",
              "Rtsne", "impute", "pheatmap", "limma", "AnnotationHub",
              "tidyverse", "ggrepel", "clusterProfiler", "msigdbr",
              "fgsea", "enrichplot", "CytoMethIC", "RColorBrewer")
invisible(lapply(packages, require, character.only = TRUE))

## PREPROCESSING__________________________________________________________
# txt <- as.data.frame(read_tsv(here("data", "20241203_rna_log2cpmfiltered.txt")))
# rownames(df) <- df$geneID
# df <- df[,-1] #22545   187
# write.csv(df, here("data", "20250113_rna_log2cpmfiltered.csv"))

## ________
mtx <- read.csv(here("data", "20250113_rna_log2cpmfiltered.csv"), row.names=1)
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
sum(ss$Sample_ID %in% colnames(mtx)) #78
# ss$Sample_ID[!ss$Sample_ID %in% colnames(mtx)]
# [1] "THY0150T" "THY0158T" "THY0164T" "THY0330T"

mtx = mtx[, colnames(mtx) %in% ss$Sample_ID] #22545    84
mtx <- mtx %>% unique() #22516    76

# code to replace names of columns in mtx
# colnames(mat)[colnames(mat) == "old_name"] <- "new_name"

## DIFFERENTIAL EXPRESSION __________________________________________________________
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx"))
mtx <- read.csv(here("data", "20250113_rna_log2cpmfiltered.csv"), row.names = 1) %>%
    .[, colnames(.) %in% ss$Sample_ID] %>%
    unique()  # 22516 genes × 84 samples
ss <- ss %>%
    filter(Include_In_Analysis == "1") %>%
    filter(Sample_ID %in% colnames(mtx)) #84 x 34
mtx <- mtx[,match(ss$Sample_ID, colnames(mtx))]

invasiveness <- factor(ss[match(colnames(mtx), ss$Sample_ID),]$Invasiveness,
                       levels = c("Low", "High"))
design <- model.matrix(~ invasiveness)
colnames(design) <- c("Intercept", "InvasivenessHigh")

fit <- mtx %>%
    lmFit(design) %>%
    eBayes()

deg_all <- topTable(fit, coef="InvasivenessHigh", n=Inf) #22516 x 6

# saveRDS(deg_all, here("data", "20250204_thyroid_deg_all.rds"))

##  _____________
deg_all <- readRDS(here("data", "20250204_thyroid_deg_all.rds"))

deg_stringent <- deg_all %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1)  # 1727 x 6
# saveRDS(deg_stringent, here("data", "20250204_thyroid_deg_stringent.rds"))

## DIFFERENTIAL METHYLATION ANALYSIS _______________________________________
res <- readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds")) #855375      7
res <- res %>% filter(
    p.adjust(Pval_InvasivenessHigh, method = "BH") < 0.05,  # Benjamini-Hochberg adjustment
    Eff_Invasiveness > 0.1
)
dim(res) # 45341 x 7

hyper_enrichment <- testEnrichment(
    res$Probe_ID[res$Est_InvasivenessHigh > 0.2],
    "TFBSconsensus",
    platform = "EPIC",
    universe = res$Probe_ID
) #783  15

# plotDot(hyper_enrichment, n_min = 20) # check
# saveRDS(hyper_enrichment, here("diff_meth", "20250118_hyper_tfbs.rds"))

hyper_enrichment <- readRDS(here("diff_meth", "20250118_hyper_tfbs.rds"))
# 783  15
tf_hyper <- prepareTF(hyper_enrichment, 783, 100, 0.05) %>%
    filter(FDR < 0.05) %>%
    pull(dbname) %>%
    unique()
# 2 lol

# hypomethylated TFs
hypo_enrichment <- testEnrichment(
    res$Probe_ID[res$Est_InvasivenessHigh < -0.2],
    "TFBSconsensus",
    platform = "EPIC",
    universe = res$Probe_ID
)
# plotDot(hypo_enrichment, n_min = 20)
# saveRDS(hypo_enrichment, here("diff_meth", "20250118_hypo_tfbs.rds"))
#
hypo_enrichment <- readRDS(here("diff_meth", "20250118_hypo_tfbs.rds"))
tf_hypo <- prepareTF(hypo_enrichment, 783, 100, 0.05) %>%
    filter(FDR < 0.05) %>%
    pull(dbname) %>%
    unique()
# 495 hypo TF

## TFBS HEATMAP LOL _____________________________________________________________
deg_hyper_expressed <- deg_stringent %>%
    filter(rownames(.) %in% tf_hyper)  #0
deg_hypo_expressed <- deg_stringent %>%
    filter(rownames(.) %in% tf_hypo) #25

deg_tf_expr <- mtx[unique(c(rownames(deg_hypo_expressed), rownames(deg_hyper_expressed))), ]
# 25 84
# saveRDS(deg_tf_expr, here("data", "20250204_deg_tf_expr.rds"))

deg_tf_expr <- readRDS(here("data", "20250204_deg_tf_expr.rds")) #25 84
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1") %>%
    filter(Sample_ID %in% colnames(deg_tf_expr)) #84 34
Invasiveness <- as.factor(ss$Invasiveness)
Cluster_Group <- as.factor(ss$Cluster_Group)

annotation_col <- data.frame(
    Cluster_Group = Cluster_Group,
    Invasiveness = Invasiveness,
    row.names = colnames(deg_tf_expr)
)

invasiveness_colors <- c(
    "High" = "#f8766dff",
    "Low" = "#00bfc4ff"
)
cluster_colors <- c(
    "HI" = "#f8766dff",
    "HIL" = "#f7b456ff",
    "LI" = "#00bfc4ff"
)

annotation_colors <- list(
    Cluster_Group = cluster_colors,
    Invasiveness = invasiveness_colors
)

pdf(here("figures", "20250202_thyroid_heatmap_rnaseq.pdf"), width=8, height=6, onefile=FALSE)
pheatmap(deg_tf_expr,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         color = colorRampPalette(rev(brewer.pal(11, "PiYG")))(50)
)
dev.off()
# color = colorRampPalette(brewer.pal(11, "PuOr"))(50)
# color = colorRampPalette(brewer.pal(11, "PiYG"))(50)

## TFBS VOLCANO PLOT_____________________________________________________________
deg_all <- readRDS(here("data", "20250204_thyroid_deg_all.rds")) # 22516 x 6
deg_all <- deg_all %>%
    mutate(
        significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                             "Significant", "Not Significant"),
        gene = rownames(.)
    )

hyper_enrichment <- readRDS(here("diff_meth", "20250118_hyper_tfbs.rds"))
hypo_enrichment <- readRDS(here("diff_meth", "20250118_hypo_tfbs.rds"))

tf_hyper <- prepareTF(hyper_enrichment, 783, 100, 0.05) %>%
    filter(FDR < 0.05) %>%
    distinct(dbname) %>%  # More efficient than unique()
    pull()

tf_hypo <- prepareTF(hypo_enrichment, 783, 100, 0.05) %>%
    filter(FDR < 0.05) %>%
    distinct(dbname) %>%
    pull()

deg_all <- deg_all %>%
    mutate(
        tf_status = case_when(
            gene %in% tf_hyper ~ "Hyper TF",
            gene %in% tf_hypo ~ "Hypo TF",
            TRUE ~ "Other"
        )
    )

pdf(here("figures", "20250204_volcano_tfbs_rnaseq.pdf"), width=6.5, height=6, onefile=FALSE)
volcano_plot <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(data = subset(deg_all, tf_status == "Other"),
               aes(color = tf_status), alpha = 0.25, size = 0.3) +
    geom_point(data = subset(deg_all, tf_status != "Other"),
               aes(color = tf_status), alpha = 1, size = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#757575ff") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#757575ff") +
    scale_color_manual(values = c("Hyper TF" = "#d61525ff", "Hypo TF" = "#0c66bcff", "Other" = "grey")) +
    theme_minimal() +
    labs(x = "log2(Fold Change)",
         y = "-log10(Adj. P-Val)",
         color = "TF Status")
significant_tfs <- deg_all[deg_all$tf_status != "Other" &
                               deg_all$adj.P.Val < 0.05 &
                               abs(deg_all$logFC) > 1, ]
volcano_plot = volcano_plot +
    geom_text_repel(data = significant_tfs,
                    aes(label = gene),
                    max.overlaps = 6,
                    color = "#4d4d4dff")
plot(volcano_plot)
dev.off()


## GSEA PLOTTING _____________________________________________________________
deg_all <- readRDS(here("data", "20250204_thyroid_deg_all.rds")) #22516     6
#  -log10(P) * sign(FC)
ranked_genes <- deg_all %>%
    mutate(
        ranking_metric = -log10(P.Value) * sign(logFC),
        gene_id = rownames(.)
    ) %>%
    arrange(desc(ranking_metric))

gene_ranks <- setNames(ranked_genes$ranking_metric, ranked_genes$gene_id)
h_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    mutate(collection = "H")

msig_list <- split(all_sets$gene_symbol, all_sets$gs_name)

gsea_results <- fgsea(
    pathways = msig_list,
    stats = gene_ranks,
    scoreType = 'std',
    minSize = 15,
    maxSize = 500,
    nPermSimple = 10000,
    eps = 0
) %>%
    as_tibble() %>%  #
    arrange(padj)
dim(gsea_results) #50  8

gsea_results_filtered <- gsea_results %>%
    filter(padj < 0.1) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = 20) %>%
    mutate(
        pathway = pathway %>%
            str_remove("^HALLMARK_") %>%
            str_replace_all("_", " ")
    )
# saveRDS(gsea_results_filtered, here("data", "20250205_gsea_results_filtered.rds"))
##
gsea_results_filtered <- readRDS(here("data", "20250205_gsea_results_filtered.rds"))

pdf(here("figures", "20250206_gsea_enrichment_dot_hallmark.pdf"), width=6.5, height=4, onefile=FALSE)
p <- ggplot(gsea_results_filtered,
            aes(x = NES,
                y = reorder(pathway, NES),
                size = -log10(padj),
                color = NES)) +
    geom_point() +
    scale_color_gradientn(colors = viridis(50)) +
    scale_size_continuous(name = "-log10(padj)") +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "Pathway")
plot(p)
dev.off()






