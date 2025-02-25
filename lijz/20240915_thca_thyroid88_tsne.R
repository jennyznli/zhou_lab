packages <- c("tidyr", "dplyr", "plotly", "readr", "writexl", "readxl",
              "here", "pvclust", "DESeq2", "stringr", "ggplot2", "sesame",
              "Rtsne", "impute", "pheatmap", "limma", "AnnotationHub",
              "tidyverse", "ggrepel", "clusterProfiler", "msigdbr",
              "fgsea", "enrichplot", "CytoMethIC")
invisible(lapply(packages, require, character.only = TRUE))

# PREPROCESSING THCA BETAS _________________________________________________________________________

ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")
thca = load(file = "/Users/jennyzli/Documents/HPC_share/thyroid/THCA.rda")
betas <- get("betas")

# rownames are weird, truncate them
colnames(betas) <- substr(colnames(betas), 1, 15)
betas = betas[, colnames(betas) %in% ss$Sample_ID] #485577    496
saveRDS(betas, "~/Documents/HPC_share/thyroid/20240915_thca_betas.rds")

# PREPROCESSING _________________________________________________________________________

betas <- readRDS("~/Documents/HPC_share/thyroid/20240915_thca_betas.rds")

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
# Before:  485577 rows and  496 columns.
# After:  395802 rows and  496 columns.

betas1 <- impute.knn(betas, k = 10, rng.seed=1234)
betas_knn <- betas1$data
betas_knn = t(betas_knn) #88 855437

rs <- grep('rs', colnames(betas_knn))
betas1 = betas_knn[,-rs] # 496 395737
betas1 = t(betas1)
saveRDS(betas1, "~/Documents/HPC_share/thyroid/20240915_thca_betas_processed.rds")

# TSNE THCA_________________________________________________________________________
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")
betas <- readRDS("~/Documents/HPC_share/thyroid/data/20240915_thca_betas_processed.rds")

bSubMostVariable <- function(betas, n=3000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

mtx = bSubMostVariable(betas, 10000) #395737  496
mtx = t(mtx)
# mtx88 = mtx[rownames(mtx) %in% ss$Sample_ID,]
mtx_reord <- mtx[match(ss$Sample_ID, rownames(mtx)), ]

##
set.seed(12345)
tsne = Rtsne(mtx_reord, dims=2, perplexity=23)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")

ss$tSNE1 <- df$tSNE1
ss$tSNE2 = df$tSNE2

pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_invasiveness.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Invasiveness, text = paste("Sample:", Patient_ID))) +
    # ggtitle("THCA tSNE") +
    # scale_shape_manual(values = c(20, 17)) +  # Use circle (16) and diamond (18) shapes
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray") +
    viridis::scale_color_viridis(discrete = TRUE)
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_tcga_subgroup.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = TCGA_subgroup, text = paste("Sample:", Patient_ID))) +
    # ggtitle("THCA tSNE") +
    # scale_shape_manual(values = c(20, 17)) +  # Use circle (16) and diamond (18) shapes
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )

plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_tcgacluster.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Methylation_Cluster_TCGA, text = paste("Sample:", Patient_ID))) +
    # ggtitle("THCA tSNE") +
    # scale_shape_manual(values = c(20, 17)) +  # Use circle (16) and diamond (18) shapes
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )

plot(p)
dev.off()

custom_colors <- c("BRAF V600E" = "#FF5733",
                   "Kinase Fusion" = "#3498DB",
                   "Ras-like" = "#3448DB",
                   "Indeterminate" = "lightgray")


pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_driver.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Genetic_Driver_simplified, text = paste("Sample:", Patient_ID))) +
    scale_color_manual(values = custom_colors) +  # Use custom colors
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )

plot(p)
dev.off()

## CONTINUOUS ONES

ss$ERK_Score <- as.numeric(ss$ERK_Score)
ss$Differentiation_Score <- as.numeric(ss$Differentiation_Score)
ss$BRS_Score <- as.numeric(ss$BRS_Score)

pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_ERK.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = ERK_Score, text = paste("Sample:", Patient_ID))) +
    scale_color_gradient(low = "blue", high = "red") +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )

plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_differentiation.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Differentiation_Score, text = paste("Sample:", Patient_ID))) +
    scale_color_gradient(low = "blue", high = "red") +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )

plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20240915_thca_tsne_BRS.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = BRS_Score, text = paste("Sample:", Patient_ID))) +
    scale_color_gradient(low = "blue", high = "red") +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )

plot(p)
dev.off()


# THCA AND TCGA CO EMBEDDING________________________________________________________________________
## GET OVERLAPPING PROBES
probes = list(
    HM450 = sesameAnno_buildAddressFile("~/Documents/HPC_share/thyroid/prev_data/HM450.hg38.manifest.tsv.gz")$Probe_ID,
    EPICv2 = sesameAnno_buildAddressFile("~/Documents/HPC_share/thyroid/prev_data/EPICv2.hg38.manifest.tsv.gz")$Probe_ID)

probes$EPICv2_prefix = sapply(strsplit(probes$EPICv2,"_"), function(x) x[1])

probesOv = with(probes, intersect(HM450, EPICv2_prefix))
table(substr(probesOv,1,2))
saveRDS(probesOv, file="~/Documents/HPC_share/thyroid/prev_data/20240917_Overlap_HM450_EPICv2.rds")

## OVERLAPPING CGS ________________________________
ovs = readRDS("~/Documents/HPC_share/thyroid/prev_data/20240917_Overlap_HM450_EPICv2.rds") #394442
ovs = sort(ovs[substr(ovs,1,2)=="cg"]) #391466
saveRDS(ovs, file="~/Documents/HPC_share/thyroid/prev_data/20240917_Overlap_HM450_EPICv2_CG.rds")
# write_tsv(tibble(Probe_ID=ovs), file="~/Documents/HPC_share/thyroid/EPICv2_HM450_CG_391466.tsv")

## OVERLAPPING BETAS ________________________________
probesOv = readRDS("~/Documents/HPC_share/thyroid/prev_data/20240917_Overlap_HM450_EPICv2_CG.rds")  #391466

betas88 = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240320_thyroid136_betas_condensed.rds")[probesOv,]
betas496 = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240915_thca_betas.rds")[probesOv,]
ss88 = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx") %>% dplyr::filter(Include_In_Analysis == "1")
ss496 = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx") %>% dplyr::filter(!(Invasiveness == "NA")) #449  19

sel_cols <- c("Source", "Sample_ID", "Cluster_Group", "Invasiveness", "Driver_Group")
# sel_cols <- c("Source", "Sample_ID", "Cluster_Group", "Invasiveness", "Driver_Group")
ss88 <- ss88[, sel_cols]
ss496 <- ss496[, sel_cols]
ss_comb <- rbind(ss88, ss496)
# write_xlsx(ss_comb, "~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_537.xlsx")
######
ss_comb <- read_excel("~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_537.xlsx")

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
betas88 <- cleanMatrixForClusterW(betas88)
betas496 <- cleanMatrixForClusterW(betas496)

betas881 <- impute.knn(betas88, k = 10, rng.seed=1234)$data
# betas881 = t(betas881)
betas4961 <- impute.knn(betas496, k = 10, rng.seed=1234)$data
# betas4961 = t(betas4961)

common_probes <- intersect(rownames(betas881), rownames(betas4961)) #334430

betas882 <- betas881[common_probes, colnames(betas881) %in% ss88$IDAT] #333047     88
betas4962 <- betas4961[common_probes, colnames(betas4961) %in% ss496$Sample_ID] #333047    449

betas_comb = cbind(betas882, betas4962) #333047    537
# saveRDS(betas_comb, here("data", "20250203_thyroid88_thca_combined_processed_betas.rds"))

## LEUKOCYTE FRACTION _______________________________________________________________
betas_comb <- readRDS(here("data", "20250203_thyroid88_thca_combined_processed_betas.rds"))
leuko_comb = cmi_predict(betas_comb, cmi_model)

df_comb <- as.data.frame(sapply(leuko_comb, function(x) x))
rownames(df_comb) <- ss_comb$Sample_ID
colnames(df_comb) <- c("Leukocyte_Fraction")
ss_comb$Leukocyte_Fraction = df_comb$Leukocyte_Fraction

saveRDS(df_comb, here("data", "20250203_thyroid88_thca_leuko.rds"))
df_comb = readRDS(here("data", "20250203_thyroid88_thca_leuko.rds"))

## df = readRDS("~/zhoulab_zhouw3/20241004_OS_with_leuko.rds")
## df_tsne = readRDS("~/zhoulab_zhouw3/20241004_OS_tSNE.rds") |> mutate(prognosis_necrosis = as.factor(as.integer(prognosis_necrosis)), perc_leukocyte2 = df$leuko[match(sample, df$name)])
## pdf('~/gallery/20250104_OS_leuko.pdf', family="ArialMT", width=5, height=3.5, onefile=FALSE)
## ggplot(df_tsne) + geom_point(aes(tSNE1, tSNE2, color=perc_leukocyte2)) + scale_color_gradientn(colors=coolwarm(20))
## dev.off()

## MAKE TSNE ________________________________
betas_comb <- readRDS("~/Documents/HPC_share/thyroid/data/20250203_thyroid88_thca_combined_processed_betas.rds")
mtx = bSubMostVariable(betas_comb, 3000)
mtx = t(mtx) #537 3000

set.seed(12345)
tsne = Rtsne(mtx, dims=2, perplexity=24)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
rownames(df) = rownames(mtx)
saveRDS(df, here("data", "2025_thyroid88_thca_tsne_coords.rds")) #537   2
ss_comb <- cbind(ss_comb, df)
# saveRDS(ss_comb, here("data", "2025_thyroid88_thca_tsne537_df.rds")) #537   8




# df$Invasiveness = ss88$Invasiveness[match(rownames(df), ss88$IDAT)]
# df$Driver_Group = ss88$Driver_Group[match(rownames(df), ss88$IDAT)]
# df$Cluster_Group = ss88$Cluster_Group[match(rownames(df), ss88$IDAT)]
#
# df$Study = ifelse(is.na(df$Invasiveness),"THCA","Pediatric")
# df$Invasiveness[is.na(df$Invasiveness)] = ss496$Invasiveness[match(rownames(df), ss496$Sample_ID)][is.na(df$Invasiveness)]
# df$Driver_Group[is.na(df$Driver_Group)] = ss496$Driver_Group[match(rownames(df), ss496$Sample_ID)][is.na(df$Driver_Group)]
# df$Cluster_Group[is.na(df$Cluster_Group)] = ifelse(
#     df$Invasiveness[is.na(df$Cluster_Group)] == "High",
#     "HI",
#     ifelse(df$Invasiveness[is.na(df$Cluster_Group)] == "Low",
#            "LI",
#            NA)
# )

## PLOT TSNE ________________________________
ss_comb <- readRDS(here("data", "2025_thyroid88_thca_tsne537_df.rds")) #537   8
# df = df %>% filter(!Driver_Group == "Indeterminate") %>% filter(!is.na(Cluster_Group)) #496   6
## sort through for braf only !!
# df_rev <- df %>%
#     filter(Driver %in% c("BRAF V600E", "Ras-like"))
# sum(ss88 %>% filter(Driver_Group == c("Ras-like", "BRAF V600E")))

custom_shapes <- c("THCA" = 16,  # Circle
                   "PED" = 1)  # Diamond
cluster_colors <- c(
    "HI" = "#f8766dff",
    "HIL" = "#f7b456ff",
    "LI" = "#00bfc4ff"
)

driver_colors <- c(
    "BRAF V600E" = "#ff6cc3ff",    # red
    "Kinase Fusion" = "#20bb20ff",   # blue
    "Ras-like" = "#00b4f0ff",   # blue
    "DICER1" = "#b47cffff",
)

pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_invasiveness.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss_comb, aes(x = tSNE1, y = tSNE2, color = Cluster_Group, shape = Source)) +
    scale_color_manual(values = cluster_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_drivers.pdf', width=6.3, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss_comb, aes(x = tSNE1, y = tSNE2, color = Driver_Group, shape = Source)) +
    scale_color_manual(values = driver_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_leuko.pdf', width=6.1, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss_comb, aes(x = tSNE1, y = tSNE2, color = Leukocyte_Fraction, shape = Source)) +
    scale_shape_manual(values = custom_shapes) +
    scale_color_gradientn(colors=gnuplot(20)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", color = NA),  # Added color = NA
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

# pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_leuko.pdf', width=6.3, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = ss_comb, aes(x = tSNE1, y = tSNE2, color = Leukocyte_Fraction, shape = Source)) +
#     scale_shape_manual(values = custom_shapes) +
#     scale_color_gradientn(colors=gnuplot(20)) +
#     theme_minimal() +  # This gives a clean white background
#     theme(
#         panel.background = element_rect(fill = "white", color = "black"),
#         plot.background = element_rect(fill = "white", color = NA),
#         legend.title = element_blank(),
#         panel.grid = element_blank()  # This removes the grid lines completely
#         # If you want subtle grid lines, use:
#         # panel.grid.major = element_line(color = "gray90", size = 0.2),
#         # panel.grid.minor = element_blank()
#     )
# plot(p)
# dev.off()



# pdf('~/Documents/HPC_share/thyroid/figures/20250112_merged_tsne_drivers_braf_ras.pdf', width=6, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = df_rev, aes(x = tSNE1, y = tSNE2, color = Driver, shape = Study)) +
#     scale_shape_manual(values = custom_shapes) +  # Set custom shapes for each study
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()
#
# pdf('~/Documents/HPC_share/thyroid/figures/20250112_merged_tsne_invasiveness.pdf', width=6, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = df_rev, aes(x = tSNE1, y = tSNE2, color = Invasiveness, shape = Study)) +
#     scale_shape_manual(values = custom_shapes) +  # Set custom shapes for each study
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()


## AIMES REVISED VERSION

# probesOv = readRDS("~/Documents/HPC_share/thyroid/prev_data/20240917_Overlap_HM450_EPICv2_CG.rds") #394442
#
# betas88 = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240320_thyroid136_betas_condensed.rds")[probesOv,]
# betas496 = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240915_thca_betas.rds")[probesOv,]
# ss88 = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")
# ss496 = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")
#
# ss88_rev <- ss88 %>%
#     filter(Driver_Group %in% c("BRAF V600E", "Ras-like"))
# ss496_rev <- ss496 %>%
#     filter(Driver_Group %in% c("BRAF V600E", "Ras-like"))
#
# betas88_rev = betas88[, colnames(betas88) %in% ss88_rev$IDAT] #379648     44
# betas496_rev = betas496[, colnames(betas496) %in% ss496_rev$Sample_ID] #391466    387
#
# cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
#     cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
#                 f_row, f_col))
#     cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
#     namtx = is.na(mtx)
#     good_row = rowSums(namtx) <= ncol(mtx) * f_row
#     good_col = colSums(namtx) <= nrow(mtx) * f_col
#     cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
#     mtx[good_row, good_col]
# }
# betas88_rev <- cleanMatrixForClusterW(betas88_rev)
# betas496_rev <- cleanMatrixForClusterW(betas496_rev)
#
# betas88_rev1 <- impute.knn(betas88_rev, k = 10, rng.seed=1234)$data
# betas496_rev1 <- impute.knn(betas496_rev, k = 10, rng.seed=1234)$data
#
# common_probes <- intersect(rownames(betas88_rev1), rownames(betas496_rev)) #334430
#
# betas88_rev2 <- betas88_rev1[common_probes,]
# betas496_rev2 <- betas496_rev1[common_probes,]
#
# betas_comb = cbind(betas88_rev2, betas496_rev2) #33047 now    632
# mtx = bSubMostVariable(betas_comb, 3000)
# mtx = t(mtx)
#
# ##
# set.seed(12345)
# tsne = Rtsne(mtx, dims=2, perplexity=21)
# df = as.data.frame(tsne$Y)
# colnames(df) = c("tSNE1", "tSNE2")
#
# # ss_comb = cbind(ss88, ss496)
# # ss_comb = cbind(ss88, ss496)
# rownames(df) = rownames(mtx)
#
# df$Invasiveness = ss88$Invasiveness[match(rownames(df), ss88$IDAT)]
# df$Driver = ss88$Driver_Group[match(rownames(df), ss88$IDAT)]
# df$Study = ifelse(is.na(df$Invasiveness),"THCA","Pediatric")
# df$Invasiveness[is.na(df$Invasiveness)] = ss496$Invasiveness[match(rownames(df), ss496$Sample_ID)][is.na(df$Invasiveness)]
# df$Driver[is.na(df$Driver)] = ss496$Driver_Group[match(rownames(df), ss496$Sample_ID)][is.na(df$Driver)]
#
# custom_shapes <- c("THCA" = 16,  # Circle
#                    "Pediatric" = 1)  # Diamond
#
# pdf('~/Documents/HPC_share/thyroid/figures/20250203_merged_tsne_drivers_braf_ras_2.pdf', width=6, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = df, aes(x = tSNE1, y = tSNE2, color = Driver, shape = Study)) +
#     scale_shape_manual(values = custom_shapes) +  # Set custom shapes for each study
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()
#
# png('~/Documents/HPC_share/thyroid/figures/20250112_merged_tsne_drivers_braf_ras_2.png', width=6, height=5, units='in', res=300)
#
# p <- ggplot() +
#     geom_point(data = df, aes(x = tSNE1, y = tSNE2, color = Driver, shape = Study)) +
#     scale_shape_manual(values = custom_shapes) +  # Set custom shapes for each study
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()
#
# ##
#
# pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_invasiveness_brafras_2.pdf', width=6, height=5, onefile=FALSE)
# p <- ggplot() +
#     geom_point(data = df, aes(x = tSNE1, y = tSNE2, color = Invasiveness, shape = Study)) +
#     scale_shape_manual(values = custom_shapes) +  # Set custom shapes for each study
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray")
#     )
# plot(p)
# dev.off()
#
# png('~/Documents/HPC_share/thyroid/figures/20250112_merged_tsne_invasiveness_brafras_2.png', width=6, height=5, units='in', res=300)
#
# # Create the plot
# p <- ggplot(df, aes(x = tSNE1, y = tSNE2, color = Invasiveness, shape = Study)) +
#     geom_point() +
#     scale_shape_manual(values = custom_shapes) +  # Set custom shapes for each study
#     theme(legend.title = element_blank(),
#           panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgray"))
#
# # Draw the plot
# print(p)
#
# # Close the PNG device
# dev.off()


## CHHOOSING COLORS _______________________________________________________________
# pal.bands(brewer.brbg, brewer.piyg, brewer.prgn, brewer.puor, brewer.rdbu,
#           brewer.rdgy, brewer.rdylbu, brewer.rdylgn, brewer.spectral,
#           main="Brewer diverging")
gnuplot()
pal.bands(gnuplot)
pal.bands(brewer.puor)
pal.bands(brewer.rdbu)
pal.bands(brewer.prgn)
pal.bands(parula)
library(pals)
library(grDevices)
plot_palette <- function(pal, name) {
    n <- length(pal)
    plot(1:n, rep(1,n), col=pal, pch=19, cex=3,
         ylim=c(0,1), type="p", main=name,
         xlab="Color Index", ylab="", yaxt="n")
}

pdf('~/Documents/HPC_share/thyroid/figures/20250202_pal_puor.pdf', width=10, height=8, onefile=FALSE)
all_pals <- ls("package:pals")[!grepl("^\\.", ls("package:pals"))]
pal_functions <- all_pals[sapply(all_pals, function(x) is.function(get(x, "package:pals")))]
par(mfrow=c(4,1))
for(pal_name in pal_functions) {
    pal <- try(get(pal_name, "package:pals")(), silent=TRUE)
    if(!inherits(pal, "try-error") && is.character(pal)) {
        plot_palette(pal, pal_name)
    }
}
dev.off()







