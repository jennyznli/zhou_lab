packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here", "pvclust",
              "stringr", "ggplot2", "sesame", "Rtsne", "impute", "pheatmap", "gprofiler",
              "RColorBrewer", "ggpubr", "CytoMethIC")
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

# PREPROCESSING________________________________________________________________________

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

# DIFF METH ________________________________________________________________________
betas <- t(readRDS(here("diff_meth", "20241029_thyroid88_betas_processed.rds")))

ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")

ss$Actual_Age <- as.numeric(ss$Actual_Age)
ss$Sex <- as.factor(ss$Sex)
ss$Invasiveness <- as.factor(ss$Invasiveness)

se = SummarizedExperiment(betas, colData = ss)

se_ok = (checkLevels(assay(se), colData(se)$Sex) &
             checkLevels(assay(se), colData(se)$Invasiveness))
sum(se_ok) #855375 same as betas
se = se[se_ok,]

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Invasiveness <- relevel(factor(colData(se)$Invasiveness), "Low")

# se = se[1:10,]

smry = DML(se, ~Invasiveness + Sex + Actual_Age)
res = summaryExtractTest(smry)
head(res)

saveRDS(smry, here("diff_meth", "20250217_thyroid88_smry_invasiveness_condensed.rds"))
saveRDS(res, here("diff_meth", "20250217_thyroid88_res_invasiveness_condensed.rds"))

## COLOR KEYS________________________________________________________________________
methylation_colors <- c(
    "Hypo" = "#0c66bcff",
    "Hyper" = "#d61525ff"
)
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

age_colors <- c(
    "Female" = "#8c44f6ff",
    "Male" = "#00cd92ff"
)

t_colors <- c(
    "T0" = "#98ff66ff",
    "T1" = "#82f128ff",
    "T1a" = "#01ed00ff",
    "T1b" = "#01cc00ff",
    "T2" = "#05e0c3ff",
    "T3" = "#00bbffff",
    "T3a" = "#147bffff",
    "T3b" = "#4337ffff",
    "T4" = "#343399ff",
    "T4a" = "#c30cecff"
)

n_colors <- c(
    "N0" = "#ffe93fff",
    "N1a" = "#ff7e04ff",
    "N1b" = "#fc4211ff"
)

m_colors <- c(
    "M0" = "#fc59bcff",
    "M1" = "#a708ddff"
)



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

# FIG 1A-C. TSNE PLOT _______________________________________________________________
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
df <- readRDS(here("data", "20250204_thyroid88_tsne_coords.rds"))
ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2
ss$Classifier_Accuracy <- as.factor(ss$Classifier_Accuracy)
ss$Actual_Age <- as.numeric(ss$Actual_Age)
ss$Confidence <- as.factor(ss$Confidence)
custom_shapes <- c("0" = 16,
                   "1" = 17)
PLOT_PATH <- "~/Documents/HPC_share/thyroid/figures/"
PLOT_DATE <- "20250224"
PLOT_SIZE <- list(width = 5.5, height = 4.5)

common_theme <- theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "lightgray"),
    aspect.ratio = 1,
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "none"  # Changed this line to remove legend
)

create_tsne_plot <- function(data, color_var, color_values = NULL, shape_var = NULL, shape_values = NULL, continuous = FALSE) {
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

    if (continuous) {
        p <- p + scale_color_gradientn(colors=gnuplot(20))
    } else {
        p <- p + scale_color_manual(values = color_values)
    }

    return(p)
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
         color_values = sex_colors),
    list(name = "leukocyte_fraction",
         color_var = "Leukocyte_Fraction",
         continuous = TRUE),
    list(name = "T",
         color_var = "T",
         color_values = t_colors),
    list(name = "N",
         color_var = "N",
         color_values = n_colors),
    list(name = "M",
         color_var = "M",
         color_values = m_colors),
    list(name = "age",
         color_var = "Actual_Age",
         continuous = TRUE)
)

for (config in plot_configs) {
    plot <- create_tsne_plot(
        data = ss,
        color_var = config$color_var,
        color_values = config$color_values,
        shape_var = config$shape_var,
        shape_values = config$shape_values,
        continuous = if (!is.null(config$continuous)) config$continuous else FALSE
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

## CALCULATE LEUKOCYTE ___________________
# betas = t(readRDS(here("data", "20241029_thyroid88_betas_processed.rds")))
# ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
#     filter(Include_In_Analysis == "1")
# cmi_model <- readRDS(here("prev_data", "LeukoFrac_EPIC_20240614.rds"))
# leuko = cmi_predict(betas, cmi_model)
#
# df <- as.data.frame(sapply(leuko, function(x) x))
# rownames(df) <- colnames(betas)
# colnames(df) <- c("Leukocyte_Fraction")
# write.csv(df, here("data", "20250203_thyroid88_leuko.csv"))

# FIG 2. DIFF METHYLATION ANALYSIS _______________________________________________________________


# FIG 2A. VOLCANO PLOT _______________________________________________________________
# 20250214_thyroid88_res_invasiveness_uncondensed.rds

res <- readRDS(here("data", "20250214_thyroid88_res_invasiveness_uncondensed.rds")) #855375      9
# res <- readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds"))
res$Pval_InvasivenessHigh <- p.adjust(res$Pval_InvasivenessHigh, method = "BH")
res$threshold <- "NS"
res$threshold[res$Pval_InvasivenessHigh < 0.05 & res$Est_InvasivenessHigh > 0.2] <- "Up"
res$threshold[res$Pval_InvasivenessHigh < 0.05 & res$Est_InvasivenessHigh < -0.2] <- "Down"

sig_points <- subset(res, threshold != "NS")
nonsig_points <- subset(res, threshold == "NS")

p <- ggplot() +
    geom_bin2d(data = nonsig_points,
               aes(x = Est_InvasivenessHigh, y = -log10(Pval_InvasivenessHigh)),
               bins = 250) +  # Reduced number of bins to make each bin more visible
    scale_fill_gradient(low = "lightgray",
                        high = "black",  # Darker grey for better contrast
                        trans = "log10") + # Add log transform to better show density differences
    geom_point(data = sig_points,
               aes(x = Est_InvasivenessHigh, y = -log10(Pval_InvasivenessHigh),
                   color = threshold),
               size = 0.1,
               alpha = 0.3) +
    scale_color_manual(values = c(
        "Down" = "#0c66bcff",
        "Up" = "#d61525ff"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#575757ff") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "#575757ff") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#575757ff") +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey90"),
        legend.position = "right",
        legend.title = element_blank()
    ) +
    labs(
        x = "Estimate",
        y = "-log10(FDR)"
    ) +
    guides(fill = "none")

pdf(here("figures", "20250214_thyroid88_dm_volcano_covariates.pdf"), width=5, height=5, onefile=FALSE)
plot(p)
dev.off()

## FIG 2B. HEATMAP UNSUPERVISED CLUSTERING DIFF METHYLATED___________________________________
betas = readRDS(here("data", "20241029_thyroid88_betas_processed.rds")) #88 855375
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
res <- readRDS(here("diff_meth", "20240704_thyroid88_res_invasiveness.rds")) #855375      7
res$Pval_InvasivenessHigh <- p.adjust(res$Pval_InvasivenessHigh, method = "BH")

sig_probes <- res %>%
    dplyr::filter(Pval_InvasivenessHigh < 0.05,
                  abs(Est_InvasivenessHigh) > 0.2) %>%
    mutate(Methylation = ifelse(Est_InvasivenessHigh > 0.2, "Hyper", "Hypo"))

annotation_row <- data.frame(
    Methylation = factor(sig_probes$Methylation), #13644     8
    row.names = sig_probes$Probe_ID
)
annotation_col <- data.frame(
    Invasiveness = factor(ss$Invasiveness),
    Cluster_Group = factor(ss$Cluster_Group),
    row.names = ss$Sample_ID
) %>% arrange(Invasiveness)

rownames(betas) <- ss$Sample_ID[match(rownames(betas), ss$IDAT)]
sel_betas <- t(betas[,sig_probes$Probe_ID ])  #  13644    88
sel_betas_reord <- sel_betas[,match(rownames(annotation_col), colnames(sel_betas))]

annotation_colors <- list(
    Methylation = methylation_colors,
    Invasiveness = invasiveness_colors,
    Cluster_Group = cluster_colors
)

pdf(here("figures", "20250208_thyroid88_dm_heatmap.pdf"), width = 13, height = 12, onefile = FALSE)
pheatmap(sel_betas_reord,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#0070ffff", "white", "#f70000ff"))(50),
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_legend = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         legend = FALSE)
dev.off()


## FIG S2A. GLOBAL MEANS BOXPLOT - CLUSTERS________________________________________________
betas = readRDS(here("data", "20241029_thyroid88_betas_processed.rds")) #88 855375
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")
ss$Global_Means = rowMeans(betas)

wcc <- pairwise.wilcox.test(ss$Global_Means, ss$Cluster_Group,
                            p.adjust.method = "BH")

pdf(here("figures", "20250202_thyroid88_globalmeans_group_wc.pdf"), width=4.5, height=4.5, onefile=FALSE)
p <- ggplot(ss, aes(x = Cluster_Group, y = Global_Means, fill = Cluster_Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("HI", "HIL"),
                                          c("HI", "LI"),
                                          c("HIL", "LI")),
                       label = "p.format") +
    scale_fill_manual(values = cluster_colors) +
    labs(x = "Cluster",
         y = "Global Mean Methylation",
         fill = "Cluster") +
    theme_minimal()
plot(p)
dev.off()

## FIG S2B. BARGRAPH DIFF METH TOTAL PROBES________________________________________________
# 20250210_thyroid88_smry_invasiveness_uncondensed.rds

# ADJUSTED P VALUES (BH)
res1 <- readRDS(here("diff_meth", "20240704_thyroid88_res_cluster1_group.rds")) #855375      9
res2 <- readRDS(here("diff_meth", "20240704_thyroid88_res_cluster2_group.rds")) #855375      9

combined_pvals_res1 <- c(res1$Pval_CLUSTER_GROUPGROUP2,
                         res1$Pval_CLUSTER_GROUPGROUP3)
adjusted_pvals_res1 <- p.adjust(combined_pvals_res1, method = "BH")
n <- nrow(res1)
res1$BH_GROUP2 <- adjusted_pvals_res1[1:n]
res1$BH_GROUP3 <- adjusted_pvals_res1[(n+1):(2*n)]

combined_pvals_res2 <- c(res2$Pval_CLUSTER_GROUPGROUP1,
                         res2$Pval_CLUSTER_GROUPGROUP3)
adjusted_pvals_res2 <- p.adjust(combined_pvals_res2, method = "BH")
res2$BH_GROUP1 <- adjusted_pvals_res2[1:n]
res2$BH_GROUP3 <- adjusted_pvals_res2[(n+1):(2*n)]

# saveRDS(res1, here("diff_meth", "20250209_thyroid88_res_cluster1_group_BH.rds"))
# saveRDS(res2, here("diff_meth", "20250209_thyroid88_res_cluster2_group_BH.rds"))

## START HERE
res1 <- readRDS(here("diff_meth", "20250209_thyroid88_res_cluster1_group_BH.rds"))
res2 <- readRDS(here("diff_meth", "20250209_thyroid88_res_cluster2_group_BH.rds"))

get_sig_probes <- function(res, bh_col, est_col) {
    res %>%
        dplyr::filter(!!sym(bh_col) < 0.05,
                      abs(!!sym(est_col)) > 0.2) %>%
        mutate(Methylation = ifelse(!!sym(est_col) > 0.2, "Hyper", "Hypo"))
}

sig_results <- list(
    "HI / HIL" = get_sig_probes(res1, "BH_GROUP2", "Est_CLUSTER_GROUPGROUP2"),
    "HI / LI" = get_sig_probes(res1, "BH_GROUP3", "Est_CLUSTER_GROUPGROUP3"),
    "HIL / LI" = get_sig_probes(res2, "BH_GROUP3", "Est_CLUSTER_GROUPGROUP3")
)

data <- bind_rows(
    lapply(names(sig_results), function(group) {
        data.frame(
            Group = group,
            table(sig_results[[group]]$Methylation)
        ) %>%
            rename(DMPs = Freq) %>%
            arrange(desc(Var1))
    })
)

total_dmps <- data.frame(
    Group = names(sig_results),
    Total = sapply(sig_results, nrow),
    stringsAsFactors = FALSE
) %>%
    mutate(Label = format(Total, big.mark = ","))

data$Methylation <- factor(data$Var1, levels = c("Hyper", "Hypo"))

pdf(here("figures", "20250202_thyroid88_dm_clusters_combined.pdf"), width=4, height=4, onefile=FALSE)
p <- ggplot(data, aes(x = Group, y = DMPs/10000)) +
    geom_bar(aes(fill = Var1), stat = "identity", position = "stack") +
    scale_fill_manual(values = methylation_colors) +
    geom_text(data = total_dmps,
              aes(y = Total/10000, label = Label),
              vjust = -0.3,
              size = 3.5) +
    scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"),
                       limits = c(0, max(total_dmps$Total/10000) * 1.1),
                       breaks = seq(0, ceiling(max(total_dmps$Total/10000)), by = 2)) +
    theme_minimal() +
    theme(legend.title = element_text(size = 10),
          panel.grid.major.x = element_blank()) +
    labs(fill = "Methylation",
         x = "Cluster Comparison",
         y = expression("DMPs [" * 10^4 * "]"))
plot(p)
dev.off()

## INSERT REST HERE________________________________________________



