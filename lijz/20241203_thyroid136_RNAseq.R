packages <- c("tidyr", "dplyr", "plotly", "readr", "writexl", "readxl",
              "here", "pvclust", "DESeq2", "stringr", "ggplot2", "sesame",
              "Rtsne", "impute", "pheatmap", "limma", "AnnotationHub",
              "tidyverse", "ggrepel", "clusterProfiler", "msigdbr",
              "fgsea", "enrichplot")
invisible(lapply(packages, require, character.only = TRUE))

here()

##  LOG 2 NORMALIZED PREPROCESSING___________________________________________
txt <- as.data.frame(read_tsv(here("data", "20241203_rna_log2cpmfiltered.txt")))
rownames(df) <- df$geneID
df <- df[,-1] #22545   187
write.csv(df, here("data", "20250113_rna_log2cpmfiltered.csv"))

## ____________
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

## TSNE_____________________________________________________________
# set.seed(12345)
# pca_result <- prcomp(t(mtx), center = TRUE, scale. = TRUE)
#
# # Choose number of components to keep
# var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
# n_components <- which(var_explained >= 0.8)[1] #25
#
# # Extract PCA-reduced data
# pca_data <- pca_result$x[, 1:n_components]
#
# set.seed(12346)
# tsne <- Rtsne(
#     pca_data,
#     dims = 2,
#     perplexity = 8,
#     max_iter = 2000,
#     pca = FALSE  # No need for further PCA
# )
#
# tsne = Rtsne(t(mtx), dims=2, perplexity=9)
# df = as.data.frame(tsne$Y)
# colnames(df) = c("tSNE1", "tSNE2")
# rownames(df) <- colnames(mtx)
#
# df$Invasiveness = ss$Invasiveness[match(rownames(df), ss$Sample_ID)]
# df$Driver = ss$Driver_Group[match(rownames(df), ss$Sample_ID)]
#
# pdf('~/Documents/HPC_share/thyroid/figures/20250113_tsne_rna_invasiveness.pdf', width=6, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = df, aes(x = tSNE1, y = tSNE2, color = Invasiveness)) +
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()
#
# pdf('~/Documents/HPC_share/thyroid/figures/20250113_tsne_rna_driver.pdf', width=6, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = df, aes(x = tSNE1, y = tSNE2, color = Driver)) +
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()

## DIFFERENTIALLY EXPRESSED GENES _____________________________________________________________
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
mtx <- read.csv(here("data", "20250113_rna_log2cpmfiltered.csv"), row.names = 1) %>%
    .[, colnames(.) %in% ss$Sample_ID] %>%
    unique()  # 22516 genes × 84 samples

invasiveness <- factor(ss[match(colnames(mtx), ss$Sample_ID),]$Invasiveness,
                       levels = c("Low", "High"))
design <- model.matrix(~ invasiveness)
colnames(design) <- c("Intercept", "InvasivenessHigh")

# Fit linear model and calculate differential expression
fit <- mtx %>%
    lmFit(design) %>%
    eBayes()

# Get DEGs with different thresholds
deg_all <- topTable(fit, coef="InvasivenessHigh", n=Inf) #22516 genes

# write.csv(deg_all, here("data", "20250119_thyroid_deg_all.csv"))

##  _______________________________________
deg_all <- read.csv(here("data", "20250119_thyroid_deg_all.csv"))

deg_stringent <- deg_all %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1)  # 1727 genes

# # summary stats
# summary_stats <- list(
#     total_genes = nrow(deg_all),
#     total_sig = sum(deg_all$adj.P.Val < 0.05),
#     up_sig = sum(deg_all$adj.P.Val < 0.05 & deg_all$logFC > 0),
#     down_sig = sum(deg_all$adj.P.Val < 0.05 & deg_all$logFC < 0),
#     stringent_total = nrow(deg_stringent),
#     stringent_up = sum(deg_stringent$logFC > 0),
#     stringent_down = sum(deg_stringent$logFC < 0)
# )
# #
# cat("Differential Expression Analysis Summary:\n")
# cat("Total genes analyzed:", summary_stats$total_genes, "\n")
# cat("Significantly differential (adj.P.Val < 0.05):", summary_stats$total_sig, "\n")
# cat("- Upregulated:", summary_stats$up_sig, "\n")
# cat("- Downregulated:", summary_stats$down_sig, "\n")
# cat("\nStringent criteria (adj.P.Val < 0.05 & |logFC| > 1):", summary_stats$stringent_total, "\n")
# cat("- Upregulated:", summary_stats$stringent_up, "\n")
# cat("- Downregulated:", summary_stats$stringent_down, "\n")


top_up <- deg_stringent %>%
    filter(logFC > 0) %>%
    arrange(desc(logFC)) %>%
    head(20) %>%
    mutate(Gene = rownames(.),
           Direction = "Upregulated")
top_down <- deg_stringent %>%
    filter(logFC < 0) %>%
    arrange(logFC) %>%
    head(20) %>%
    mutate(Gene = rownames(.),
           Direction = "Downregulated")
# top_genes_plot <- bind_rows(top_up, top_down) %>%
#     mutate(Gene = factor(Gene, levels = Gene[order(logFC)]))  # Order genes by logFC

# TOP UP/DOWNREGULATED GENES PLOT ____________
pdf(here("figures", "20250119_thyroid_rnaseq_topgenes.pdf"), width=5, height=7, onefile=FALSE)
p <- ggplot(top_genes_plot,
            aes(x = Gene, y = logFC, fill = Direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    coord_flip() +  # Flip coordinates for horizontal bars
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top"
    ) +
    labs(
        title = "Top 20 Up and Down Regulated Genes",
        y = "log2 Fold Change",
        fill = "Direction"
    )
plot(p)
dev.off()

## GSEA _____________________________________________________________
#  -log10(P) * sign(FC)
ranked_genes <- deg_all %>%
    mutate(
        ranking_metric = -log10(P.Value) * sign(logFC)
    ) %>%
    arrange(desc(ranking_metric))

gene_ranks <- setNames(ranked_genes$ranking_metric, rownames(ranked_genes))

# save ranking for GSEA app verison !

# gene_ranks <- as.data.frame(gene_ranks)
# gene_ranks$genes <- rownames(gene_ranks)
# write.table(gene_ranks,
#             file = here("data", "20250118_thyroid_ranked_genes.rnk"),
#             sep = "\t",
#             row.names = TRUE,
#             col.names = TRUE,
#             quote = FALSE)

# 2. Get MSigDB collections
h_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%  # Hallmark pathways
    dplyr::select(gs_name, gene_symbol)
c2_cp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%  # Canonical pathways (KEGG, Reactome)
    dplyr::select(gs_name, gene_symbol)
c6_sets <- msigdbr(species = "Homo sapiens", category = "C6") %>%  # Oncogenic signatures
    dplyr::select(gs_name, gene_symbol)
c7_sets <- msigdbr(species = "Homo sapiens", category = "C7") %>%  # Immunologic signatures
    dplyr::select(gs_name, gene_symbol)
c8_sets <- msigdbr(species = "Homo sapiens", category = "C8") %>%  # Cell type signatures
    dplyr::select(gs_name, gene_symbol)

all_sets <- bind_rows(
    # h_sets %>% mutate(collection = "H")
    # c2_cp %>% mutate(collection = "C2_CP")
    # c6_sets %>% mutate(collection = "C6")
    # c7_sets %>% mutate(collection = "C7")
    c8_sets %>% mutate(collection = "C8")
)
# Convert to format needed for GSEA
msig_list <- split(all_sets$gene_symbol, all_sets$gs_name)

gsea <- fgsea(
    pathways = msig_list,
    stats = gene_ranks,
    scoreType = 'std',
    minSize = 15,
    maxSize = 500,
    nPermSimple = 1000
)
# 5953    8
head(gsea[order(pval), ])
# #
# collapsedPathways <- collapsePathways(gsea[order(padj)][padj < 0.05], h_sets, gene_ranks)
# mainPathways <- gsea[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

# topPathwaysUp <- gsea[ES > 0][head(order(padj), n = 20), pathway]
# topPathwaysDown <- gsea[ES < 0][head(order(padj), n = 20), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

## GSEA PLOTTING_____________________________________________________________
top_gsea <- gsea %>%
    filter(padj < 0.05) %>%
    # Add absolute NES for ranking
    mutate(absNES = abs(NES)) %>%
    # Sort by absolute NES
    arrange(desc(absNES)) %>%
    # Take top 20 pathways
    slice_head(n = 20)

pdf(here("figures", "20250118_gsea_enrichment_dot_c8.pdf"), width=9, height=7, onefile=FALSE)
p <- ggplot(top_gsea,
                           aes(x = NES,
                               y = reorder(pathway, NES),
                               size = -log10(padj),
                               color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0,
                          name = "NES") +
    scale_size_continuous(name = "-log10(padj)") +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "Pathway",
         title = "Top 20 Enriched Pathways")
plot(p)
dev.off()

## GENE ONTOLOGY_____________________________________________________________
library(clusterProfiler)
library(org.Hs.eg.db)

deg_genes <- deg_stringent %>%
    rownames()

# Convert gene names to ENTREZ IDs for GO analysis
gene_ids <- bitr(deg_genes,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
# Warning message:
#     In bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) :
#     10.19% of input gene IDs are fail to map...

go_results <- enrichGO(gene          = gene_ids$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",  # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

# dot plot of top enriched
# pdf(here("figures", "20250118_GO_dotplot.pdf"), width=10, height=8)
# dotplot(go_results, showCategory=20)
# dev.off()
#
# go_df <- as.data.frame(go_results)
# pdf(here("figures", "20250118_thyroid_GO_dotplot.pdf"), width=8, height=10)
# ggplot(go_df[1:20,],
#        aes(x = Count, y = reorder(Description, Count))) +
#     geom_point(aes(size = Count, color = p.adjust)) +
#     scale_color_gradient(low = "red", high = "blue") +
#     theme_minimal() +
#     labs(x = "Gene Count", y = "GO Term",
#          title = "Top 20 Enriched GO Terms")
# dev.off()

## TFBS CONFIRMATION_____________________________________________________________
res <- readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds")) #855375      7
# res <- res %>% dplyr::filter(Pval_InvasivenessHigh < 0.05, Eff_Invasiveness > 0.1) #55026     7
res <- res %>% filter(
    p.adjust(Pval_InvasivenessHigh, method = "BH") < 0.05,  # Benjamini-Hochberg adjustment
    Eff_Invasiveness > 0.1
)
# 45341 cpgs


# hypermethylated TFs
hyper_enrichment <- testEnrichment(
    res$Probe_ID[res$Est_InvasivenessHigh > 0.2],
    "TFBSconsensus",
    platform = "EPIC",
    universe = res$Probe_ID
) #783  15

plotDot(hyper_enrichment, n_min = 20) # check
saveRDS(hyper_enrichment, here("diff_meth", "20250118_hyper_tfbs.rds"))
#
hyper_enrichment <- readRDS(here("diff_meth", "20250118_hyper_tfbs.rds"))
# 9 hyper TF
tf_hyper <- preparePlotDF(hyper_enrichment, 783, 100, 0.05) %>%
    filter(FDR < 0.05) %>%
    pull(dbname) %>%
    unique()
# 2 hyper lol

# hypomethylated TFs
hypo_enrichment <- testEnrichment(
    res$Probe_ID[res$Est_InvasivenessHigh < -0.2],
    "TFBSconsensus",
    platform = "EPIC",
    universe = res$Probe_ID
)
plotDot(hypo_enrichment, n_min = 20)
saveRDS(hypo_enrichment, here("diff_meth", "20250118_hypo_tfbs.rds"))
#
hypo_enrichment <- readRDS(here("diff_meth", "20250118_hypo_tfbs.rds"))
tf_hypo <- preparePlotDF(hypo_enrichment, 783, 100, 0.05) %>%
    filter(FDR < 0.05) %>%
    pull(dbname) %>%
    unique()
# 495 hypo TF

# 5. Expression Analysis of Methylated TFs ----------------------------------
# check_expression <- function(genes, expression_matrix, min_prop = 0.2, min_exp = 1) {  # Add min_exp parameter
#     genes <- genes[!is.na(genes)]
#     matches <- match(genes, rownames(expression_matrix))
#     valid_matches <- !is.na(matches)
#     expr_mat <- expression_matrix[matches[valid_matches], ]
#
#     # Add minimum expression threshold along with proportion filter
#     expressed <- rownames(expr_mat)[
#         apply(expr_mat, 1, median, na.rm=TRUE) > min_exp &  # Check median expression
#             rowSums(expr_mat > 0, na.rm=TRUE) >= (min_prop * ncol(expr_mat))
#     ]
#     return(expressed)
# }
#
# expressed_hypo_tfs <- check_expression(tf_hypo, mtx, min_prop = 0.2, min_exp = 1)  #413
# expressed_hyper_tfs <- check_expression(tf_hyper, mtx, min_prop = 0.2, min_exp = 1) #1

# 6. Identify Differentially Expressed TFs ------------------------------------------------------
# x <- deg_all %>%
#     filter(rownames(.) %in% tf_hyper) #2/2 TF in the RNAseq
# x <- deg_all %>%
#     filter(rownames(.) %in% tf_hypo) #469/495 TF in the RNAseq

deg_hyper_expressed <- deg_stringent %>%
    filter(rownames(.) %in% tf_hyper)  #0
deg_hypo_expressed <- deg_stringent %>%
    filter(rownames(.) %in% tf_hypo) #25

## TFBS VOLCANO PLOT_____________________________________________________________
deg_all$significant <- ifelse(deg_all$adj.P.Val < 0.05 & abs(deg_all$logFC) > 1,
                              "Significant", "Not Significant")
deg_all$gene <- rownames(deg_all)

# Add color for TFs of interest
deg_all$tf_status <- "Other"
deg_all$tf_status[deg_all$gene %in% tf_hypo] <- "Hypo TF"
deg_all$tf_status[deg_all$gene %in% tf_hyper] <- "Hyper TF"

pdf(here("figures", "20250118_volcano_tfbs_rnaseq.pdf"), width=6, height=7, onefile=FALSE)
volcano_plot <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(data = subset(deg_all, tf_status == "Other"),
               aes(color = tf_status), alpha = 0.2, size = 0.5) +  # Lower alpha for grey points
    geom_point(data = subset(deg_all, tf_status != "Other"),
               aes(color = tf_status), alpha = 1, size = 0.5) +    # Full opacity for TFs
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    scale_color_manual(values = c("Hyper TF" = "red", "Hypo TF" = "blue", "Other" = "grey")) +
    theme_minimal() +
    labs(x = "log2 Fold Change",
         y = "-log10 Adjusted P-value",
         color = "TF Status")

# Add labels for significant TFs
significant_tfs <- deg_all[deg_all$tf_status != "Other" &
                               deg_all$adj.P.Val < 0.05 &
                               abs(deg_all$logFC) > 1, ]
volcano_plot +
    geom_text_repel(data = significant_tfs,
                    aes(label = gene),
                    max.overlaps = 20)
plot(volcano_plot)
dev.off()

## TFBS HEATMAP LOL _____________________________________________________________
deg_tf_expr <- mtx[unique(c(rownames(deg_hypo_expressed), rownames(deg_hyper_expressed))), ]
# saveRDS(deg_tf_expr, here("data", "20250204_deg_tf_expr.rds"))
deg_tf_expr <- readRDS(here("data", "20250204_deg_tf_expr.rds"))

annotation_col <- data.frame(
    Invasiveness = invasiveness,
    row.names = colnames(deg_tf_expr)
)

invasiveness_colors <- list(
    Invasiveness = c("Low" = "#00BFC4", "High" = "#F8766D")
)

# Create the heatmap
pdf(here("figures", "20250202_thyroid_heatmap_rnaseq.pdf"), width=8, height=6, onefile=FALSE)
pheatmap(deg_tf_expr,
         scale = "row",  # Scale by row (gene)
         annotation_col = annotation_col,
         annotation_colors = invasiveness_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         color = colorRampPalette(brewer.pal(11, "PuOr"))(20)
)

dev.off()

## RANDOM CODE IDEK BRUH _____________________________________________________________
# meth_data <- data.frame(
#     methylation_diff = res$Est_InvasivenessHigh,
#     significant = res$Pval_InvasivenessHigh < 0.05
# )
#
# # Plot density
# pdf(here("figures", "20250118_distribution_meth_rnaseq.pdf"), width=6, height=6, onefile=FALSE)
# p <- ggplot(meth_data, aes(x = methylation_diff, fill = significant)) +
#     geom_density(alpha = 0.5) +
#     scale_fill_manual(values = c("grey", "darkred")) +
#     theme_minimal() +
#     labs(x = "Methylation Difference (High - Low Invasiveness)",
#          y = "Density",
#          title = "Distribution of Methylation Changes",
#          fill = "Significant")
# plot(p)
# dev.off()


# tf_validation <- data.frame(
#     TF = tf_list,
#     Mean_Expression = NA,
#     Is_Expressed = NA,
#     Is_Differential = NA,
#     LogFC = NA,
#     Padj = NA
# )
#
# for(i in seq_along(tf_list)) {
#     tf <- tf_list[i]
#
#     if(tf %in% rownames(mtx)) {
#         tf_validation$Mean_Expression[i] <- rowMeans(mtx[tf,]) # not sure if this is helpful btwn invaisveness?
#         tf_validation$Is_Expressed[i] <- tf_validation$Mean_Expression[i] > 1  # Set appropriate threshold
#
#         # Check if differentially expressed but i should do all of them
#         if(tf %in% rownames(deg)) {
#             tf_row <- deg[tf,]
#             tf_validation$Is_Differential[i] <- tf_row$P.Value < 0.05
#             tf_validation$LogFC[i] <- tf_row$logFC
#             tf_validation$Padj[i] <- tf_row$P.Value
#         }
#     }
# }
#
# sum(!is.na(tf_validation$Is_Differential))
#
# ggplot(tf_validation, aes(x = reorder(TF, Mean_Expression), y = Mean_Expression)) +
#     geom_bar(stat = "identity") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = "Expression Levels of Enriched TFs",
#          x = "Transcription Factor",
#          y = "Mean Expression")
#
# # Color by differential expression
# ggplot(tf_validation %>% filter(!is.na(LogFC)),
#        aes(x = LogFC, y = -log10(Padj), label = TF)) +
#     geom_point(aes(color = abs(LogFC) > 1 & Padj < 0.05)) +
#     theme_minimal() +
#     labs(title = "Differential Expression of TFs",
#          x = "Log2 Fold Change",
#          y = "-log10(adjusted p-value)")
#
# # Summary statistics
# summary_stats <- tf_expression %>%
#     summarize(
#         Total_TFs = n(),
#         Expressed = sum(Is_Expressed, na.rm = TRUE),
#         Differential = sum(Is_Differential, na.rm = TRUE)
#     )

ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
mtx <- read.csv(here("data", "20250113_rna_log2cpmfiltered.csv"), row.names=1)
mtx = mtx[, colnames(mtx) %in% ss$Sample_ID]
mtx <- mtx %>% unique()

invasiveness <- factor(ss[match(colnames(mtx), ss$Sample_ID),]$Invasiveness,
                       levels = c("Low", "High"))  # Makes "Low" the reference level

design <- model.matrix(~ invasiveness)
colnames(design) <- c("Intercept", "InvasivenessHigh")

lm <- lmFit(mtx, design)
fit <- eBayes(lm)
deg <- topTable(fit, coef="InvasivenessHigh", n=Inf) %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1)

# diff_expr <- deg %>%
#     mutate(
#         expr_status = case_when(
#             logFC > 1 ~ "Upregulated",
#             logFC < -1 ~ "Downregulated"
#         )
#     )
# 1414    7

# 2. Methylation analysis
res <- readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness_uncondensed.rds"))
diff_meth <- res %>%
    filter(FPval_Invasiveness < 0.05, abs(Est_InvasivenessHigh) > 0.2) %>%
    mutate(
        meth_status = case_when(
            Est_InvasivenessHigh > 0.2 ~ "Hypermethylated",
            Est_InvasivenessHigh < -0.2 ~ "Hypomethylated",
            TRUE ~ NA_character_  # Handle any other cases
        )
    )

# 3. Expression analysis
diff_expr <- deg %>%
    mutate(
        expr_status = case_when(
            logFC > 1 ~ "Upregulated",
            logFC < -1 ~ "Downregulated",
            TRUE ~ NA_character_  # Handle any other cases
        )
    ) %>%
    rownames_to_column("Gene_Symbol")
# [1] 1414    8


# 4. Get and process manifest
ah <- AnnotationHub()
manifest <- ah[["AH116484"]]
probe_gene_map <- manifest %>%
    as.data.frame() %>%
    dplyr::select(
        ProbeID = IlmnID,
        Gene_Symbol = UCSC_RefGene_Name,
        Region = UCSC_RefGene_Group,
        Chromosome = CHR
    ) %>%
    filter(!is.na(Gene_Symbol)) %>%
    separate_rows(Gene_Symbol, sep = ";")

integrated_results <- diff_meth %>%
    inner_join(probe_gene_map, by = c("Probe_ID" = "ProbeID")) %>%
    inner_join(diff_expr, by = "Gene_Symbol") %>%
    filter(
        (meth_status == "Hypomethylated" & expr_status == "Upregulated") |
            (meth_status == "Hypermethylated" & expr_status == "Downregulated")
    )
# 1002   18
significant_genes <- unique(integrated_results$Gene_Symbol)

# 6. Analysis summaries
# Basic counts
print("Summary Statistics:")
print(paste("Number of differential methylation probes:", nrow(diff_meth)))
print(paste("Number of differential expression genes:", nrow(diff_expr)))
print(paste("Number of integrated results:", nrow(integrated_results)))
# [1] "Number of differential methylation probes: 13813"
# > print(paste("Number of differential expression genes:", nrow(diff_expr)))
# [1] "Number of differential expression genes: 1414"
# > print(paste("Number of integrated results:", nrow(integrated_results)))
# [1] "Number of integrated results: 1002"

# Methylation-Expression relationship table
summary_table <- table(integrated_results$meth_status, integrated_results$expr_status)
print("\nMethylation-Expression Relationships:")
print(summary_table)
# Downregulated Upregulated
# Hypermethylated            72           0
# Hypomethylated              0         930


# region_summary <- integrated_results %>%
#     group_by(Region) %>%
#     summarise(
#         n_probes = n(),
#         n_genes = n_distinct(Gene_Symbol),
#         avg_meth_change = mean(Est_InvasivenessHigh),  # Fixed from meth_status
#         avg_expr_change = mean(logFC)                  # Fixed from AveExpr
#     )
#
# integrated_results_promoter <- integrated_results %>%
#     filter(Region %like% "TSS" | Region %like% "Promoter")
# #
# # Add summary by genomic region
# region_summary <- integrated_results %>%
#     group_by(Region, meth_status, expr_status) %>%
#     summarize(count = n(), .groups = 'drop')

# 7. Optional: Save results
write.csv(integrated_results,
          here("data", "20250116_methylation_expression_integration.csv"),
          row.names = FALSE)

# Using clusterProfiler for enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Get genes showing inverse relationship

# Run GO enrichment
ego <- enrichGO(gene = significant_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP")  # Biological Process

# KEGG pathway analysis
# First convert gene symbols to ENTREZ IDs
genes_entrez <- bitr(significant_genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
kegg_results <- enrichKEGG(gene = genes_entreSz$ENTREZID)

###
# hyper_down_genes <- integrated_results %>%
#     filter(meth_status == "Hypermethylated",
#            expr_status == "Downregulated") %>%
#     pull(Gene_Symbol) %>%
#     unique()
#
# hypo_up_genes <- integrated_results %>%
#     filter(meth_status == "Hypomethylated",
#            expr_status == "Upregulated") %>%
#     pull(Gene_Symbol) %>%
#     unique()
#
# # Analyze each set
# ego_hyper_down <- enrichGO(gene = hyper_down_genes,
#                            OrgDb = org.Hs.eg.db,
#                            keyType = "SYMBOL",
#                            ont = "BP")
#
# ego_hypo_up <- enrichGO(gene = hypo_up_genes,
#                         OrgDb = org.Hs.eg.db,
#                         keyType = "SYMBOL",
#                         ont = "BP")
#
# hyper_down_df <- as.data.frame(ego_hyper_down)
# hypo_up_df <- as.data.frame(ego_hypo_up)
#
# top_pathways <- rbind(
#     hyper_down_df %>%
#         top_n(10, -p.adjust) %>%
#         mutate(group = "Hyper/Down"),
#     hypo_up_df %>%
#         top_n(10, -p.adjust) %>%
#         mutate(group = "Hypo/Up")
# )
#
#
# ggplot(top_pathways,
#        aes(x = -log10(p.adjust),
#            y = reorder(Description, -log10(p.adjust)),
#            fill = group)) +
#     geom_bar(stat = "identity") +
#     theme_minimal() +
#     labs(x = "-log10(adjusted p-value)",
#          y = "Pathway",
#          title = "Top Enriched Pathways by Methylation/Expression Pattern")
#
#
#
# invasion_terms <- c("adhesion", "migration", "invasion", "EMT",
#                     "matrix", "metastasis")
#
# invasion_pathways <- ego@result %>%
#     as.data.frame() %>%
#     filter(grepl(paste(invasion_terms, collapse = "|"), Description,
#                  ignore.case = TRUE))
#
# # Plot invasion-specific pathways
# ggplot(invasion_pathways,
#        aes(x = -log10(p.adjust),
#            y = reorder(Description, -log10(p.adjust)))) +
#     geom_bar(stat = "identity", fill = "steelblue") +
#     theme_minimal() +
#     labs(x = "-log10(adjusted p-value)",
#          y = "Invasion-related Pathway",
#          title = "Enriched Invasion-related Pathways")
#
#

ranked_list <- integrated_results %>%
    arrange(desc(logFC)) %>%
    pull(logFC, name = Gene_Symbol)

# Method 2: Combined metric (if you want to incorporate both methylation and expression)
ranked_list_combined <- integrated_results %>%
    mutate(
        combined_score = logFC * -sign(Est_InvasivenessHigh)
    ) %>%
    arrange(desc(combined_score)) %>%
    pull(combined_score, name = Gene_Symbol)

# 2. Get gene sets
# Hallmark gene sets
h_gene_sets <- msigdbr(species = "Homo sapiens",
                       category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(x = .$gene_symbol, f = .$gs_name)

# Cancer and EMT specific sets
c6_gene_sets <- msigdbr(species = "Homo sapiens",
                        category = "C6") %>%  # Oncogenic signatures
    dplyr::select(gs_name, gene_symbol) %>%
    split(x = .$gene_symbol, f = .$gs_name)

# 3. Run GSEA
gsea_hallmark <- GSEA(
    geneList = ranked_list_combined,
    TERM2GENE = msigdbr(species = "Homo sapiens",
                        category = "H") %>%
        dplyr::select(gs_name, gene_symbol),
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)

# 4. Visualizations
# A. Enrichment plot for specific pathways
gseaplot2(gsea_hallmark,
          geneSetID = 1:3,  # Plot top 3 pathways
          pvalue_table = TRUE,
          title = "GSEA Enrichment Plot")

# B. Ridgeplot for multiple pathways
ridgeplot(gsea_hallmark,
          showCategory = 10) +
    labs(x = "Enrichment Distribution")

# C. GSEA table visualization
dotplot(gsea_hallmark,
        showCategory = 20,
        split = ".sign") +  # Split by direction of enrichment
    ggtitle("Hallmark Pathways GSEA")

# D. Network of enriched terms
emapplot(gsea_hallmark,
         showCategory = 15)

# 5. Focus on invasion-related pathways
invasion_terms <- c("EMT", "INVASION", "METASTASIS", "MIGRATION",
                    "ADHESION", "MAPK", "RAS")

invasion_results <- gsea_hallmark@result %>%
    as.data.frame() %>%
    filter(grepl(paste(invasion_terms, collapse = "|"),
                 Description, ignore.case = TRUE))

# 6. Extract leading edge genes
# These are the genes driving the enrichment
leading_edge <- getLeadingEdge(gsea_hallmark)

# Create summary of leading edge genes
leading_edge_summary <- data.frame(
    pathway = names(leading_edge),
    genes = sapply(leading_edge, paste, collapse = ", ")
)



