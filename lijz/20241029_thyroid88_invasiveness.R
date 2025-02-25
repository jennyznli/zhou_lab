packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here", "pvclust",
              "stringr", "ggplot2", "sesame", "Rtsne", "impute", "pheatmap", "gprofiler")
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

# PREPROCESSING 88________________________________________________________________________

betas = readRDS(here("data", "20240320_thyroid136_betas_condensed.rds"))
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
betas = betas[, colnames(betas) %in% ss$IDAT] #930659   88

betas <- betas %>%
    cleanMatrixForClusterW() %>%
    impute.knn(k = 10, rng.seed = 1234)
# Filter rows with >0.50 missingness and columns with >0.50 missingness.
# Before:  930659 rows and  136 columns.
# After:  856603 rows and  136 columns.

betas_knn <- t(betas$data) #88 855437
rs_probes <- grep('rs', colnames(betas_knn))
betas_filtered = betas_knn[,-rs_probes] #88 855375

saveRDS(betas_filtered, here("data", "20241029_thyroid88_betas_processed.rds"))

# 88 TSNE ANALYSIS________________________________________________________________________
betas = t(readRDS(here("data", "20241029_thyroid88_betas_processed.rds")))
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")

mtx = t(bSubMostVariable(betas, 3000))

set.seed(12345678)
pca_result <- prcomp(mtx, center = TRUE, scale. = TRUE)
var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
n_components <- which(var_explained >= 0.8)[1] #25
pca_data <- pca_result$x[, 1:n_components]

set.seed(123467)
tsne <- Rtsne(
    pca_data,
    dims = 2,
    perplexity = 8,
    max_iter = 2000,
    pca = FALSE
)

df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2

# saveRDS(df, here("data", "20250204_thyroid88_tsne_coords.rds"))

# INVASIVENESS 88 PLOT _______________________________________________________________
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
df <- readRDS(here("data", "20250204_thyroid88_tsne_coords.rds"))
ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2

ss$Classifier_Accuracy <- as.factor(ss$Classifier_Accuracy)
ss$Confidence <- as.factor(ss$Confidence)

cluster_colors <- c(
    "HI" = "#f8766dff",
    "HIL" = "#f7b456ff",
    "LI" = "#00bfc4ff"
)

invasiveness_colors <- c(
    "High" = "#f8766dff",
    "Low" = "#00bfc4ff"
)

driver_colors <- c(
    "BRAF V600E" = "#ff6cc3ff",
    "Kinase Fusion" = "#20bb20ff",
    "Ras-like" = "#00b4f0ff",
    "DICER1" = "#b47cffff"
)

sex_colors <- c(
    "Female" = "#8c44f6ff",
    "Male" = "#00cd92ff"
)

custom_shapes <- c("0" = 16,
                   "1" = 17)

PLOT_PATH <- "~/Documents/HPC_share/thyroid/figures/"
PLOT_DATE <- "20250204"
PLOT_SIZE <- list(width = 6, height = 5)

common_theme <- theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "bottom"
)

create_tsne_plot <- function(data, color_var, color_values, shape_var = NULL, shape_values = NULL) {
    p <- ggplot(data = data, aes(x = tSNE1, y = tSNE2)) +
        coord_fixed() +
        common_theme

    if (!is.null(shape_var)) {
        p <- p + geom_point(aes(color = !!sym(color_var),
                                shape = !!sym(shape_var)),
                            size = 2) +
            scale_shape_manual(values = shape_values)
    } else {
        p <- p + geom_point(aes(color = !!sym(color_var)),
                            size = 2, shape = 16)
    }

    p + scale_color_manual(values = color_values)
}

plot_configs <- list(
    list(name = "invasiveness_accuracy",
         color_var = "Invasiveness",
         color_values = invasiveness_colors,
         shape_var = "Classifier_Accuracy",
         shape_values = c(16, 17)),
    list(name = "drivergroup",
         color_var = "Driver_Group",
         color_values = driver_colors),
    list(name = "sex",
         color_var = "Sex",
         color_values = sex_colors)
)

for (config in plot_configs) {
    plot <- create_tsne_plot(
        data = ss,
        color_var = config$color_var,
        color_values = config$color_values,
        shape_var = config$shape_var,
        shape_values = config$shape_values
    )

    filename <- file.path(PLOT_PATH,
                          sprintf("%s_tsne88_%s.pdf",
                                  PLOT_DATE,
                                  config$name))

    ggsave(filename = filename,
           plot = plot,
           width = PLOT_SIZE$width,
           height = PLOT_SIZE$height,
           units = "in",
           device = "pdf")
}


## HEATMAP UNSUPERVISED CLUSTERING DIFF METHYLATED___________________________________

res3 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_invasiveness.rds")

dim(res3 %>% dplyr::filter(FPval_Invasiveness < 0.05, Est_InvasivenessHigh > 0.2)) #1145
dim(res3 %>% dplyr::filter(FPval_Invasiveness < 0.05,  Est_InvasivenessHigh < -0.2)) #12500

pos = (res3 %>% dplyr::filter(FPval_Invasiveness < 0.05, Est_InvasivenessHigh > 0.2))$Probe_ID #1145
neg = (res3 %>% dplyr::filter(FPval_Invasiveness < 0.05,  Est_InvasivenessHigh < -0.2))$Probe_ID #12500

sel = c(pos, neg) #13645

annotation_col <- data.frame(
    Methylation = factor(c(rep("Hyper", 1145), rep("Hypo", 12500)))
)
rownames(annotation_col) = sel

sel_betas <- betas1[,colnames(betas1) %in% sel] #88 13645
rownames(sel_betas) = ss$Sample_ID

annotation_row <- data.frame(
    Invasiveness = factor(ss$Invasiveness)
)
rownames(annotation_row) <- ss$Sample_ID
annotation_row = annotation_row %>%
    mutate(Invasiveness) %>%
    arrange(Invasiveness)

sel_betas_reord = sel_betas[match(rownames(annotation_row), rownames(sel_betas)),]
sel_betas_reord = sel_betas_reord[,match(rownames(annotation_col), colnames(sel_betas))]

## HEATMAP UNSUPERVISED CLUSTERING DIFF METHYLATED___________________________________
#
# pdf('~/Documents/HPC_share/thyroid/figures/20240904_thyroid88_dm_heatmap.pdf', width=10, height=13)
# pheatmap(sel_betas_reord, cluster_cols = FALSE, cluster_rows = FALSE,
#          color = colorRampPalette(c("blue", "white", "red"))(50),
#          annotation_row = annotation_row,
#          annotation_col = annotation_col,
#          show_rownames = TRUE, show_colnames = FALSE # Show row names and hide column names
#          # main = "Differentially Methylated Probes Heatmap"
#          )
# dev.off()

# clustered version
pdf('~/Documents/HPC_share/thyroid/figures/20240904_thyroid88_dm_heatmap_clustered.pdf', width=10, height=13)
pheatmap(sel_betas_reord, cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         show_rownames = TRUE, show_colnames = FALSE # Show row names and hide column names
         # main = "Differentially Methylated Probes Heatmap"
)
dev.off()

# HEATMAP UNSUPERVISED CLUSTERING MOST VARIABLE (UNUSED)___________________________________
betas = readRDS(file = "~/Documents/HPC_share/thyroid/20240320_thyroid136_betas_condensed.rds") #930659
ss = read_excel("~/Documents/HPC_share/thyroid/20231102_thyroid_master.xlsx")

ss = ss %>% dplyr::filter(INCLUDE_IN_ANALYSIS == "1")
betas = betas[, colnames(betas) %in% ss$IDAT] #930659   88

betas <- cleanMatrixForClusterW(betas)
# Filter rows with >0.50 missingness and columns with >0.50 missingness.
# Before:  930659 rows and  136 columns.
# After:  856603 rows and  136 columns.

betas1 <- impute.knn(betas, k = 10, rng.seed=1234)
betas_knn <- betas1$data
betas_knn = t(betas_knn) #88 855437

rs <- grep('rs', colnames(betas_knn))
betas1 = betas_knn[,-rs] #88 855375
betas1 = t(betas1)

bSubMostVariable <- function(betas, n=3000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}
# take 1% most variable
mtx = bSubMostVariable(betas1, 8553)
mtx = t(mtx)

pv_model <- pvclust(mtx, method.hclust = "ward.D2", method.dist = "euclidean")

pdf('~/Documents/HPC_share/thyroid/figures/20240917_thyroid88_unsupervised_heatmap.pdf', width=10, height=13)
plot(pv_model)
pvrect(pv_model, alpha = 0.95)  # Add rectangles around significant clusters
dev.off()

# VOLCANO PLOT - DOUBLE CHECK THIS_______________________________________________________________

res3 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_invasiveness.rds")
# sel_vol = res3 %>% dplyr::filter(Eff_Invasiveness > 0.3)

pdf('~/Documents/HPC_share/thyroid/figures/20240703_volcano_global.pdf', width=6, height=5, onefile=FALSE)
sel_vol = res3[sample(nrow(res3), 80000), ]
g <- ggplot(sel_vol) + geom_point(aes(Est_InvasivenessHigh, -log10(Pval_InvasivenessHigh)), size=.2) +
    # ggtitle("Volcano Plot") +
    theme(panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(g)
dev.off()

## INVASIVENESS ENRICHMENT GRAPHS ________________________________________________
# x = testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh > 0.2], "TFBSconsensus", platform="EPIC", universe=res3$Probe_ID)


res3 = readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds"))
# Low v High TFBS HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_invasiveness_TFBS_enrichment_hyper.pdf', width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh > 0.2], "TFBSconsensus", platform="EPIC", universe=res3$Probe_ID), n_min = 40)
dev.off()


# preparePlotDF <- function(df, n_min, n_max, max_fdr) {
#     # filters out nD is 0 or less
#     df <- df[df$nD >0,]
#     # zero values in the FDR column are replaced with the smallest positive floating-point number representable in R
#     df$FDR[df$FDR==0] <- .Machine$double.xmin
#     # sort by FDR - show stronger stat significance
#     df <- df[order(df$FDR),]
#     # filtering which ones show up
#     if (sum(df$FDR < max_fdr) < n_min) {
#         df1 <- head(df, n=n_min)
#     } else {
#         df <- df[df$estimate > 0,] # enrichment only, exclude depletion
#         df1 <- df[df$FDR < max_fdr,]
#         df1 <- head(df1, n=n_max)
#     }
#     # extract group names
#     gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
#     if ("Target" %in% colnames(df1)) {
#         df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$dbname)
#     } else {
#         df1$db1 <- sprintf("%s: %s", gp, df1$dbname)
#     }
#     df1$db1 <- factor(df1$db1, levels=rev(df1$db1))
#     df1$dbname <- factor(df1$dbname, levels=rev(df1$dbname))
#
#     # Add pathway grouping
#     pathway_groups <- list(
#         "Hippo" = c("TEAD1", "YAP1", "WWTR1"),
#         "AP1" = c("JUN", "JUNB", "FOS", "FOSL1", "FOSL2"),
#         "Other" = NULL  # will catch all others
#     )
#
#     # Assign pathway groups
#     df1$pathway <- "Other"
#     for (path_name in names(pathway_groups)) {
#         if (!is.null(pathway_groups[[path_name]])) {
#             df1$pathway[df1$dbname %in% pathway_groups[[path_name]]] <- path_name
#         }
#     }
#
#     # Order factors to group related TFs together
#     df1$pathway <- factor(df1$pathway, levels=c("Hippo", "AP1", "Other"))
#     df1$dbname <- reorder(df1$dbname, as.numeric(df1$pathway))
#
#     return(df1)
# }
#
# plotDot <- function(df, n_min = 10, n_max = 30, max_fdr = 0.05) {
#     db1 <- FDR <- overlap <- estimate <- NULL
#     df1 <- preparePlotDF(df, n_min, n_max, max_fdr)
#
#     ggplot(df1) +
#         # Add subtle background for pathway groups
#         geom_rect(aes(fill = pathway),
#                   xmin = -Inf, xmax = Inf,
#                   ymin = -Inf, ymax = Inf,
#                   alpha = 0.1) +
#         # Main points
#         geom_point(aes(x = dbname, y = -log10(FDR),
#                        size = overlap,
#                        color = estimate)) +
#         # Facet by pathway
#         facet_grid(pathway ~ ., scales = "free_y", space = "free") +
#         # Styling
#         coord_flip() +
#         scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
#         scale_fill_manual(values = c(
#             "Hippo" = "#E64B35",
#             "AP1" = "#4DBBD5",
#             "Other" = "grey90"
#         )) +
#         theme_minimal() +
#         theme(
#             panel.grid.major.y = element_blank(),
#             strip.text = element_text(face = "bold"),
#             legend.position = "right"
#         ) +
#         labs(
#             title = "Enriched TFBS by Pathway",
#             y = "-log10(FDR)",
#             x = "",
#             caption = "Point size indicates overlap; color indicates enrichment estimate"
#         )
# }
## PLOT ESTIMATE
preparePlotDF1 <- function(df, n_min, n_max, max_fdr) {
    df <- df[df$nD > 0,]
    df$FDR[df$FDR == 0] <- .Machine$double.xmin

    # Sort by estimate instead of FDR
    df <- df[order(abs(df$FDR), decreasing = TRUE),]

    if (sum(df$FDR < max_fdr) < n_min) {
        df1 <- head(df, n = n_min)
    } else {
        df1 <- df[df$FDR < max_fdr,]
        df1 <- head(df1, n = n_max)
    }

    gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
    if ("Target" %in% colnames(df1)) {
        df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$dbname)
    } else {
        df1$db1 <- sprintf("%s: %s", gp, df1$dbname)
    }

    # Factor levels based on estimate order
    df1$db1 <- factor(df1$db1, levels = rev(df1$db1))
    df1$dbname <- factor(df1$dbname, levels = rev(df1$dbname))
    df1
}

plotDot1 <- function(df, n_min = 20, n_max = 20, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF1(df, n_min, n_max, max_fdr)

    ggplot(df1) +
        geom_point(aes(x = estimate, y = dbname, size = overlap, color = -log10(FDR))) +
        scale_color_gradient(low = "#FF0000", high = "#0000FF") +
        xlab("Estimate") +
        ylab("") +
        theme_minimal()
}

## PLOTTING LOG FDR
preparePlotDF2 <- function(df, n_min, n_max, max_fdr) {
    df <- df[df$nD >0,]
    df$FDR[df$FDR==0] <- .Machine$double.xmin
    df <- df[order(df$FDR),]
    if (sum(df$FDR < max_fdr) < n_min) {
        df1 <- head(df, n=n_min)
    } else {
        df <- df[df$estimate > 0,] # enrichment only, exclude depletion
        df1 <- df[df$FDR < max_fdr,]
        df1 <- head(df1, n=n_max)
    }

    gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
    if ("Target" %in% colnames(df1)) {
        df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$dbname)
    } else {
        df1$db1 <- sprintf("%s: %s", gp, df1$dbname)
    }
    df1$db1 <- factor(df1$db1, levels=rev(df1$db1))
    df1$dbname <- factor(df1$dbname, levels=rev(df1$dbname))

    df1
}

plotDot2 <- function(df, n_min = 20, n_max = 20, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(dbname, estimate, size=overlap, color=-log10(FDR))) +
        coord_flip() +
        scale_color_gradient(low="#FF0000",high="#0000FF") +
        ylab("Estimate (OR)") + xlab("") +
        theme_minimal()  # Adjust the base font size if needed
}

res3 = readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds"))

x <- testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh > 0.2],
                    "TFBSconsensus",
                    platform = "EPIC",
                    universe = res3$Probe_ID)

y <- testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh < -0.2],
                    "TFBSconsensus",
                    platform = "EPIC",
                    universe = res3$Probe_ID)
res3$FDR <- p.adjust(res3$Pval_InvasivenessHigh, method = "BH")



#ESTIMATE
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_invasiveness_TFBS_enrichment_hypo_est.pdf',
    width = 3.5, height = 3.75, onefile = FALSE)
plotDot(x, n_max = 20)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_invasiveness_TFBS_enrichment_hyper_est.pdf',
    width = 3.5, height = 3.75, onefile = FALSE)
plotDot(y,n_min = 20)
dev.off()


#FDR
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_invasiveness_TFBS_enrichment_hypo.pdf',
    width = 3.25, height = 3.75, onefile = FALSE)
plotDot(x, n_max = 20)
dev.off()


pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_invasiveness_TFBS_enrichment_hyper.pdf',
    width = 3.25, height = 3.75, onefile = FALSE)
plotDot(y,n_min = 20)
dev.off()







# Low v High TFBS HYPER OLD
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_invasiveness_TFBS_enrichment_hypo.pdf', width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh < -0.2], "TFBSconsensus", platform="EPIC", universe=res3$Probe_ID), n_min = 40)
dev.off()

# Low v High HM HYPO OLD
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_invasiveness_HM_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh > 0.2], "HMconsensus", platform="EPIC", universe=res3$Probe_ID) ) + ggtitle("Low v High HM Hypo")
dev.off()

# Low v High HM HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_invasiveness_HM_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh < -0.2], "HMconsensus", platform="EPIC", universe=res3$Probe_ID) ) + ggtitle("Low v High HM Hyper")
dev.off()

# Low v High TISSUE HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_invasiveness_tissue_enrichment_hypo.pdf', family="ArialMT", width=6, height=7, onefile=FALSE)
plotDot(testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh > 0.2], "tissueSignature", platform="EPIC", universe=res3$Probe_ID) ) + ggtitle("Low v High Tissue Signature Hypo")
dev.off()

# Low v High TISSUE HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_invasiveness_tissue_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res3$Probe_ID[res3$Est_InvasivenessHigh < -0.2], "tissueSignature", platform="EPIC", universe=res3$Probe_ID) ) + ggtitle("Low v High Tissue Signature Hyper")
dev.off()

## KYCG EXPERIMENTATION _________________________________________________________________________

res3 = readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds"))

query <- names(sesameData_getProbesByGene("TTF1", "EPICv2"))
results <- testEnrichment(res3,
                          KYCG_buildGeneDBs(query, max_distance=100000, platform="EPICv2"),
                          platform="EPICv2")
results[,c("dbname","estimate","gene_name","FDR", "nQ", "nD", "overlap")]

## gene ontology analysis
betas = readRDS(here("data", "20240320_thyroid136_betas.rds"))

## go/pathway enrichment
query <- res3$Probe_ID
regs <- sesameData_getTxnGRanges("hg38", merge2gene = TRUE)
genes <- sesameData_annoProbes(query, regs, platform="EPICv2", return_ov_features=TRUE)
genes

# select_probes <- function(probes1, probes2) {
#     # Check each rowname against all probes
#     matches <- sapply(probes1, function(x) {
#         any(startsWith(x, probes2))
#     })
#     return(matches)
# }
# x <- rownames(betas)
# selected_probes <- select_probes(query, x)

## use gene name
gostres <- gost(genes$gene_name, organism = "mmusculus")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)

## use Ensembl gene ID, note we need to remove the version suffix
gene_ids <- sapply(strsplit(names(genes),"\\."), function(x) x[1])
gostres <- gost(gene_ids, organism = "mmusculus")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)

## set enrichment analysis
res <- testEnrichmentSEA(query, "MM285.seqContextN")
KYCG_plotSetEnrichment(res[[1]])


## ENRICHMENT OF MOST IMPORTANT?
imp = read_tsv(col_names = FALSE, "~/Documents/HPC_share/thyroid/features/20240610_features.tsv") #855437 4
colnames(imp) = c("Model", "CpG", "Accuracy", "Gini")
imp = imp %>% arrange(desc(Accuracy))

pos = imp %>% filter(Accuracy > 0)
rs <- grep('rs', pos$CpG)
# print(length(rs))
mostImp = pos[1:10000,]$CpG

mtx = betas_knn
mtx = mtx[rownames(mtx) %in% ss$IDAT, colnames(mtx) %in% mostImp]


## KYCG EXPERIMENTATION _________________________________________________________________________
res3 = readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness_uncondensed.rds"))
query <- res3$Probe_ID[res3$Est_InvasivenessHigh < -0.2]
query <- res3$Probe_ID[res3$Est_InvasivenessHigh > 0.2]

regs <- sesameData_getTxnGRanges("hg38", merge2gene = TRUE)
genes <- sesameData_annoProbes(query, regs, platform="EPICv2", return_ov_features=TRUE)
genes

## GENE ONTOLOGY_______________________________________________________________
library(gprofiler2)

## use Ensembl gene ID, note we need to remove the version suffix
gene_ids <- sapply(strsplit(names(genes),"\\."), function(x) x[1])
gostres <- gost(gene_ids, organism = "hsapiens")
gostres$result[order(gostres$result$p_value),]
p <- gostplot(gostres, interactive = FALSE)
p <- gostplot(gostres)
ggsave(here("figures" ,"20241207_thyroid136_go.pdf"), p, width=15, height=15)

## SET ENRICHMENT_______________________________________________________________










