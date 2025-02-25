x<-c("tidyr", "dplyr", "plotly", "readr", "readxl", "pvclust", "stringr", "ggplot2", "sesame", "Rtsne", "impute", "pheatmap")
lapply(x, require, character.only = TRUE)

## PREPROCESSING _______________________________________________________________
betas = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240320_thyroid136_betas_condensed.rds") #930659
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")

ss = ss %>% dplyr::filter(Include_In_Analysis == "1")

betas = betas[, colnames(betas) %in% ss$IDAT] #930659   88

cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
    cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
                f_row, f_col))
    cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
    namtx = is.na(mtx)
    good_row = rowSums(namtx) <= ncol(mtx) * f_row
    good_col = colSums(namtx) <= nrow(mtx) * f_col
    cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
    mtx[good_row, good_col]
}
betas <- cleanMatrixForClusterW(betas)
# Filter rows with >0.50 missingness and columns with >0.50 missingness.
# Before:  930659 rows and  136 columns.
# After:  856603 rows and  136 columns.

betas1 <- impute.knn(betas, k = 10, rng.seed=1234)
betas_knn <- betas1$data
betas_knn = t(betas_knn) #88 855437

rs <- grep('rs', colnames(betas_knn))
betas1 = betas_knn[,-rs] #88 855375


# NUMBER OF DIFF METH BETWEEN CLUSTERS graph_____________________________________________________________
res1 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster1_group.rds")
res2 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster2_group.rds")

# significant_cpgs <- results[results$FDR < 0.05 & abs(results$beta) > 0.2, ]

# GROUP 1 V 2
(length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2]) + length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2]))
# 117179

# GROUP 1 V 3
(length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 > 0.2]) + length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 < -0.2]))
# 21806

# GROUP 2 V 3
(length(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2]) + length(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2]))
# 122841
data <- data.frame(
    Group = c("HI / HIL", "HI / LI", "HIL / LI"),
    DMPs = c(117179, 21806, 122841),
    Label = c("117,179", "21,806", "122,841")
)

# Create the plot
p <- ggplot(data, aes(x = Group, y = DMPs / 10000)) +
    geom_bar(stat = "identity", fill = "#00BFC4") +
    geom_text(aes(label = Label), vjust = -0.3, size = 3.5) +
    scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"), limits = c(0, 13)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 12))

pdf('~/Documents/HPC_share/thyroid/figures/20240703_thyroid88_dm_clusters_overall.pdf', width=3, height=5, onefile=FALSE)
plot(p)
dev.off()

# DIFF METH BETWEEN CLUSTERS OVERALL enrichment graph_____________________________________________________________
res1 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster1_group.rds")
res2 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster1_group.rds")

# 1 v 2
res = testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v2_enrichment_hypo.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res, n_label=30) + ggtitle("Cluster 1 v. 2 Hypo")
dev.off()

res = testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2], platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v2_enrichment_hyper.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res, n_label=30) + ggtitle("Cluster 1 v. 2 Hyper")
dev.off()

# 1 v 3
res = testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 > 0.2], platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_enrichment_hypo.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res, n_label=30) + ggtitle("Cluster 1 v. 3 Hypo")
dev.off()

res = testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 < -0.2], platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_enrichment_hyper.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res, n_label=30) + ggtitle("Cluster 1 v. 3 Hyper")
dev.off()

# 2 v 3
res = testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_enrichment_hypo.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res, n_label=30) + ggtitle("Cluster 2 v. 3 Hypo")
dev.off()

res = testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_enrichment_hyper.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res, n_label=30) + ggtitle("Cluster 2 v. 3 Hyper")
dev.off()

# DIFF METH BETWEEN CLUSTERS SPECIFIC enrichment graph_____________________________________________________________

plotDot <- function(df, n_min = 10, n_max = 30, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(dbname, -log10(FDR), size=Overlap, color=Estimate)) +
        coord_flip() + ggtitle("Enriched Databases") +
        scale_color_gradient(low="blue",high="red") +
        ylab("-log10(FDR)") + xlab("")
}

preparePlotDF <- function(df, n_min, n_max, max_fdr) {
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
## _______________________________________

db = KYCG_listDBGroups("EPIC")

# 1 v 3 TFBS HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_TFBS_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 > 0.2], "TFBSconsensus", platform="EPIC", universe=res1$Probe_ID), n_min = 40) + ggtitle("Cluster 1 v. 3 TFBS Hypo")
dev.off()

# 1 v 3 TFBS HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_TFBS_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 < -0.2], "TFBSconsensus", platform="EPIC", universe=res1$Probe_ID), n_min = 40) + ggtitle("Cluster 1 v. 3 TFBS Hyper")
dev.off()

# 1 v 3 HM HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_HM_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 > 0.2], "HMconsensus", platform="EPIC", universe=res1$Probe_ID) ) + ggtitle("Cluster 1 v. 3 HM Hypo")
dev.off()

# 1 v 3 HM HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_HM_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 < -0.2], "HMconsensus", platform="EPIC", universe=res1$Probe_ID) ) + ggtitle("Cluster 1 v. 3 HM Hyper")
dev.off()

# 1 v 3 TISSUE HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_tissue_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 > 0.2], "tissueSignature", platform="EPIC", universe=res1$Probe_ID), n_min = 10) + ggtitle("Cluster 1 v. 3 Tissue Signature Hypo")
dev.off()

# 1 v 3 TISSUE HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v3_tissue_enrichment_hyper.pdf', family="ArialMT", width=6, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 < -0.2], "tissueSignature", platform="EPIC", universe=res1$Probe_ID), n_min = 10) + ggtitle("Cluster 1 v. 3 Tissue Signature Hyper")
dev.off()
## _______________________________________


# 1 v 2 TFBS HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v2_TFBS_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], "TFBSconsensus", platform="EPIC", universe=res1$Probe_ID), n_min = 40) + ggtitle("Cluster 1 v. 2 TFBS Hypo")
dev.off()


# 1 v 2 TFBS HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v2_TFBS_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2], "TFBSconsensus", platform="EPIC", universe=res1$Probe_ID), n_min = 40 ) + ggtitle("Cluster 1 v. 2 TFBS Hyper")
dev.off()

# 1 v 2 HM HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v2_HM_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], "HMconsensus", platform="EPIC", universe=res1$Probe_ID) ) + ggtitle("Cluster 1 v. 2 HM Hypo")
dev.off()

# 1 v 2 HM HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster1v2_HM_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2], "HMconsensus", platform="EPIC", universe=res1$Probe_ID) ) + ggtitle("Cluster 1 v. 2 HM Hyper")
dev.off()


## _______________________________________

# 2 v 3 TFBS HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_TFBS_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], "TFBSconsensus", platform="EPIC", universe=res2$Probe_ID), n_min = 40 ) + ggtitle("Cluster 2 v. 3 TFBS Hypo")
dev.off()

# 2 v 3 TFBS HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_TFBS_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], "TFBSconsensus", platform="EPIC", universe=res2$Probe_ID), n_min = 40 ) + ggtitle("Cluster 2 v. 3 TFBS Hyper")
dev.off()

# 2 v 3 HM HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_HM_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], "HMconsensus", platform="EPIC", universe=res2$Probe_ID) ) + ggtitle("Cluster 2 v. 3 HM Hypo")
dev.off()

# 2 v 3 HM HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_HM_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], "HMconsensus", platform="EPIC", universe=res2$Probe_ID) ) + ggtitle("Cluster 2 v. 3 HM Hyper")
dev.off()

# 2 v 3 TISSUE HYPO
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_tissue_enrichment_hypo.pdf', family="ArialMT", width=6, height=7, onefile=FALSE)
plotDot(testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], "tissueSignature", platform="EPIC", universe=res2$Probe_ID) ) + ggtitle("Cluster 2 v. 3 Tissue Signature Hypo")
dev.off()

# 2 v 3 TISSUE HYPER
pdf('~/Documents/HPC_share/thyroid/figures/20240902_thyroid88_cluster2v3_tissue_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], "tissueSignature", platform="EPIC", universe=res2$Probe_ID) ) + ggtitle("Cluster 2 v. 3 Tissue Signature Hyper")
dev.off()


# NUMBER OF DIFF METH INVASIVENESS graph_____________________________________________________________

res3 = readRDS(file = "~/Documents/HPC_share/thyroid/20240704_thyroid88_res_invasiveness.rds")


## GLOBAL MEANS FIGURE - INVASIVENESS AND BY GROUP _____________________________________________________________
library(ggplot2)
library(ggpubr) # for adding statistical test results

betas1 = t(betas1)
global_means = rowMeans(betas1)
ss$Global_Means <- global_means

# FOR INVASIVENESS HIGH AND LOW AND WILCOX TEST
wc <- wilcox.test(Global_Means ~ Invasiveness, data = ss)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_globalmeans_invasiveness.pdf', width=4, height=4, onefile=FALSE)
p <- ggplot(ss, aes(x = Invasiveness, y = GLOBAL_MEANS)) +
    geom_boxplot(aes(fill = Invasiveness), alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.3) +  # Add individual points
    scale_fill_manual(values = c("High" = "red", "Low" = "blue")) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",
                       label.y = max(ss$GLOBAL_MEANS) + 0.02) + # Adjust position of p-value
    theme_minimal()
plot(p)
dev.off()

# Test    P_value W_statistic
# W Wilcoxon 0.00311334         564

# FOR CLUSTERS
wcc <- pairwise.wilcox.test(ss$Global_Means, ss$Cluster_Group,
                     p.adjust.method = "BH")  # Benjamini-Hochberg adjustment
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_globalmeans_group_wc.pdf', width=4.5, height=4.5, onefile=FALSE)
p <- ggplot(ss, aes(x = Cluster_Group, y = Global_Means, fill = Cluster_Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
    stat_compare_means(method = "wilcox.test",
                       comparisons = list(c("HI", "HIL"),
                                          c("HI", "LI"),
                                          c("HIL", "LI")),
                       label = "p.signif") +  # Shows significance stars
    scale_fill_manual(values = c("#f8766dff", "#f7b456ff", "#00bfc4ff")) +
    theme_minimal()
plot(p)
dev.off()

#       HI      HIL
# HIL 0.00061 -
# LI  0.00022 0.49348


# TSNE ANALYSIS FOR 136 ________________________________________________________________________

betas1 = t(betas1)
bSubMostVariable <- function(betas, n=3000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

mtx = bSubMostVariable(betas1, 10000)
mtx = t(mtx)

set.seed(12345)
tsne = Rtsne(mtx, dims=2, perplexity=12)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")

ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2

## DIFF METH BETWEEN CLUSTERS______________________________________________________

res1 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster1_group.rds")
res2 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster2_group.rds")

# GROUP 1 V 2
a = c(
    sum(res1$Est_CLUSTER_GROUPGROUP2 > 0.2)/855375, #74159
    sum(res1$Est_CLUSTER_GROUPGROUP2 >= -0.2 & res1$Est_CLUSTER_GROUPGROUP2 <= 0.2)/855375, #738196
    sum(res1$Est_CLUSTER_GROUPGROUP2 < -0.2)/855375) #43020)

# GROUP 1 V 3
b = c(
    sum(res1$Est_CLUSTER_GROUPGROUP3 > 0.2)/855375,  #19380
    sum(res1$Est_CLUSTER_GROUPGROUP3 >= -0.2 & res1$Est_CLUSTER_GROUPGROUP3 <= 0.2)/855375,  #833569
    sum(res1$Est_CLUSTER_GROUPGROUP3 < -0.2)/855375) #2426

# GROUP 2 V 3
c = c(
    sum(res2$Est_CLUSTER_GROUPGROUP3 > 0.2)/855375,  #57556
    sum(res2$Est_CLUSTER_GROUPGROUP3 >= -0.2 & res2$Est_CLUSTER_GROUPGROUP3 <= 0.2)/855375,  #732534
    sum(res2$Est_CLUSTER_GROUPGROUP3 < -0.2)/855375) #65285

x <- c(a, b, c)

data <- data.frame(
    Group = rep(c("HI / HIL", "HI / LI", "HIL / LI"), each = 3),
    Methylation = rep(c("(1, 0.2)", "(0.2, -0.2)", "(-0.2, -1)"), 3),
    DMPs = c(74159, 738196, 43020, 19380, 833569, 2426, 57556, 732534, 65285),
    DMPs_percent = x
)

# Convert 'Methylation' to factor to maintain order
data$Methylation <- factor(data$Methylation, levels = c("(1, 0.2)", "(0.2, -0.2)", "(-0.2, -1)"))

# Create the plot
p <- ggplot(data, aes(x = Group, y = DMPs_percent, fill = Methylation)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("(1, 0.2)" = "#F8766D", "(0.2, -0.2)" = "#d3d3d3", "(-0.2, -1)" = "#00BFC4")) +
    labs(y = "DMP Ratio") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 12))

pdf('~/Documents/HPC_share/thyroid/figures/20240703_thyroid88_dm_clusters.pdf', width=4, height=5, onefile=FALSE)
plot(p)
dev.off()


res1 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster1_group.rds")
res2 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster2_group.rds")

# significant_cpgs <- results[results$FDR < 0.05 & abs(results$beta) > 0.2, ]

# GROUP 1 V 2
(length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2]) + length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2]))
# 117179

# GROUP 1 V 3
(length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 > 0.2]) + length(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP3 < -0.2]))
# 21806

# GROUP 2 V 3
(length(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2]) + length(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2]))
# 122841
data <- data.frame(
    Group = c("HI / HIL", "HI / LI", "HIL / LI"),
    DMPs = c(117179, 21806, 122841),
    Label = c("117,179", "21,806", "122,841")
)

# Create the plot
p <- ggplot(data, aes(x = Group, y = DMPs / 10000)) +
    geom_bar(stat = "identity", fill = "#00BFC4") +
    geom_text(aes(label = Label), vjust = -0.3, size = 3.5) +
    scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"), limits = c(0, 13)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 12))

pdf('~/Documents/HPC_share/thyroid/figures/20240703_thyroid88_dm_clusters_overall.pdf', width=3, height=5, onefile=FALSE)
plot(p)
dev.off()


## DIFF METH BETWEEN CLUSTERS (FIG S2C)______________________________________________________
res1 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster1_group.rds")
res2 = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240704_thyroid88_res_cluster2_group.rds")

data <- data.frame(
    Group = rep(c("HI / HIL", "HI / LI", "HIL / LI"), each = 2),  # Changed to 2 categories
    Methylation = rep(c("(0.2, 1)", "(-1, -0.2)"), 3),  # Removed middle category
    DMPs = c(
        # HI / HIL
        74159, 43020,  # Removed middle value
        # HI / LI
        19380, 2426,   # Removed middle value
        # HIL / LI
        57556, 65285   # Removed middle value
    )
)

total_dmps <- data.frame(
    Group = c("HI / HIL", "HI / LI", "HIL / LI"),
    Total = c(117179, 21806, 122841),
    Label = c("117,179", "21,806", "122,841")
)

data$Methylation <- factor(data$Methylation,
                           levels = c("(0.2, 1)", "(-1, -0.2)"))  # Removed middle level

pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_dm_clusters_combined.pdf',
    width = 4, height = 4, onefile = FALSE)

p <- ggplot(data, aes(x = Group, y = DMPs/10000)) +  # Changed back to DMPs/10000
    geom_bar(aes(fill = Methylation), stat = "identity", position = "stack") +
    scale_fill_manual(values = c("(0.2, 1)" = "#ea7df0ff",
                                 "(-1, -0.2)" = "#20bb20ff")) +  # Removed middle color
    geom_text(data = total_dmps,
              aes(y = Total/10000, label = Label),
              vjust = -0.3,
              size = 3.5) +
    scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"),
                       limits = c(0, 14),
                       breaks = seq(0, 14, by = 2)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 12))

plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_dm_clusters_combined2.pdf',
    width = 4, height = 4, onefile = FALSE)
p <- ggplot(data, aes(x = Group, y = DMPs/10000)) +  # Changed back to DMPs/10000
    geom_bar(aes(fill = Methylation), stat = "identity", position = "stack") +
    scale_fill_manual(values = c("(0.2, 1)" = "red",
                                 "(-1, -0.2)" = "blue")) +  # Removed middle color
    geom_text(data = total_dmps,
              aes(y = Total/10000, label = Label),
              vjust = -0.3,
              size = 3.5) +
    scale_y_continuous(name = expression("DMPs [" * 10^4 * "]"),
                       limits = c(0, 14),
                       breaks = seq(0, 14, by = 2)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 12))

plot(p)
dev.off()


## DIFF METH BTWN CLUSTERS TISSUE (FIG 2F, S2E)_____________________________________________________

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

plotDot2 <- function(df, n_min = 10, n_max = 10, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF2(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(dbname, estimate, size=overlap, color=-log10(FDR))) +
        coord_flip() +
        scale_color_gradient(low="#FF0000",high="#0000FF") +
        ylab("Estimate (OR)") + xlab("") +
        theme_minimal()  # Adjust the base font size if needed
}

# 1 v 2 / HI V HIL TISSUE HYPO
c <- testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], "tissueSignature", platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster1v2_tissue_enrichment_hyper_est.pdf',  width = 4, height = 2.3,  onefile = FALSE)
plotDot2(c, n_max = 10, n_min = 10)
dev.off()

# 1 v 2 / HI V HIL TISSUE HYPER
d <- testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2], "tissueSignature", platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster1v2_tissue_enrichment_hypo_est.pdf',  width = 5.2, height = 2.3,  onefile = FALSE)
plotDot2(d, n_max = 10, n_min = 10)
dev.off()

# 2 v 3 / HIL V LI TISSUE HYPO
a <- testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], "tissueSignature", platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster2v3_tissue_enrichment_hyper_est.pdf', width = 5.2, height = 2.3,  onefile=FALSE)
plotDot2(a, n_max = 10, n_min = 10)
dev.off()

# 2 v 3 / HIL V LI TISSUE HYPER
b <- testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], "tissueSignature", platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster2v3_tissue_enrichment_hypo_est.pdf', width = 4, height = 2.3, onefile=FALSE)
plotDot2(b, n_max = 10, n_min = 10)
dev.off()

## DIFF METH BTWN CLUSTERS TISSUE (FIG S2E-F)_____________________________________________________

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

plotDot2 <- function(df, n_min = 10, n_max = 10, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))

    df1 <- preparePlotDF2(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(dbname, estimate, size=overlap, color=-log10(FDR))) +
        coord_flip() +
        scale_color_gradient(low="#FF0000",high="#0000FF") +
        ylab("Estimate (OR)") + xlab("") +
        theme_minimal()  # Adjust the base font size if needed
}

# 1 v 2 / HI v HIL TFBS HYPER
c <- testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], "TFBS", platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster1v2_TFBS_enrichment_hyper.pdf', width = 3.5, height = 3.75,  onefile=FALSE)
plotDot(c, n_max = 20, n_min = 20)
dev.off()

# 1 v 2 / HI v HIL TFBS HYPO
d <- testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2], "TFBS", platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster1v2_TFBS_enrichment_hypo.pdf',  width = 3.5, height = 3.75, onefile=FALSE)
plotDot(d, n_max = 20, n_min = 20)
dev.off()

# 2 V 3 / HIL v LI TFBS HYPO
a <- testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], "TFBS", platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster2v3_TFBS_enrichment_hyper.pdf',width = 3.5, height = 3.75,  onefile=FALSE)
plotDot(a, n_max = 20, n_min = 20)
dev.off()

# 2 V 3 / HIL v LI TFBS HYPER
b <- testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], "TFBS", platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster2v3_TFBS_enrichment_hypo.pdf',  width = 3.5, height = 3.75,  onefile=FALSE)
plotDot(b, n_max = 20, n_min = 20)
dev.off()


## SAME EXCEPT TOP 40 INSTEAD

# 1 v 2 / HI v HIL TFBS HYPER
c <- testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], "TFBS", platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster1v2_TFBS_enrichment_hyper40.pdf', width = 3.5, height = 5.75,  onefile=FALSE)
plotDot(c, n_max = 40, n_min = 40)
dev.off()

# 1 v 2 / HI v HIL TFBS HYPO
d <- testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 < -0.2], "TFBS", platform="EPIC", universe=res1$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster1v2_TFBS_enrichment_hypo40.pdf',  width = 3.5, height = 5.75, onefile=FALSE)
plotDot(d, n_max = 40)
dev.off()

# 2 V 3 / HIL v LI TFBS HYPO
a <- testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 > 0.2], "TFBS", platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster2v3_TFBS_enrichment_hyper40.pdf',width = 3.5, height = 5.75,  onefile=FALSE)
plotDot(a, n_max = 40, n_min = 40)
dev.off()

# 2 V 3 / HIL v LI TFBS HYPER
b <- testEnrichment(res2$Probe_ID[res2$Est_CLUSTER_GROUPGROUP3 < -0.2], "TFBS", platform="EPIC", universe=res2$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20250202_thyroid88_cluster2v3_TFBS_enrichment_hypo40.pdf',  width = 3.5, height = 5.75,  onefile=FALSE)
plotDot(b, n_max = 40, n_min = 40)
dev.off()



