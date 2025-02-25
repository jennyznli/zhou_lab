packages <- c("tidyr", "dplyr", "plotly", "readr", "writexl", "readxl",
              "here", "pvclust", "DESeq2", "stringr", "ggplot2", "sesame",
              "Rtsne", "impute", "pheatmap", "limma", "AnnotationHub",
              "tidyverse", "ggrepel", "clusterProfiler", "msigdbr",
              "fgsea", "enrichplot", "CytoMethIC", "RColorBrewer")
invisible(lapply(packages, require, character.only = TRUE))

# PREPROCESSING _________________________________________________________________________
probes = list(
    HM450 = sesameAnno_buildAddressFile("~/Documents/HPC_share/thyroid/prev_data/HM450.hg38.manifest.tsv.gz")$Probe_ID,
    EPICv2 = sesameAnno_buildAddressFile("~/Documents/HPC_share/thyroid/prev_data/EPICv2.hg38.manifest.tsv.gz")$Probe_ID)

probes$EPICv2_prefix = sapply(strsplit(probes$EPICv2,"_"), function(x) x[1])

probesOv = with(probes, intersect(HM450, EPICv2_prefix)) #394442
table(substr(probesOv,1,2))
saveRDS(probesOv, file="~/Documents/HPC_share/thyroid/prev_data/20240917_overlap_HM450_EPICv2.rds")

## OVERLAPPING CGS ________________________________________________________________________
ovs = readRDS("~/Documents/HPC_share/thyroid/prev_data/20240917_Overlap_HM450_EPICv2.rds") #394442
ovs = sort(ovs[substr(ovs,1,2)=="cg"]) #391466
saveRDS(ovs, file="~/Documents/HPC_share/thyroid/prev_data/20240917_overlap_HM450_EPICv2_CG.rds")

## OVERLAPPING BETAS AND CLEAN SAMPLESHEETS________________________________________________________________________
probesOv = readRDS("~/Documents/HPC_share/thyroid/prev_data/20240917_overlap_HM450_EPICv2_CG.rds")  #391466

betas88 = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240320_thyroid136_betas_condensed.rds")[probesOv,]
betas496 = readRDS(file = "~/Documents/HPC_share/thyroid/data/20240915_thca_betas.rds")[probesOv,]
ss88 = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")
ss496 = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")

ss88 = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx") %>%
    dplyr::filter(Include_In_Analysis == "1")
ss411 = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx") %>%
    dplyr::filter(!(Driver_Group == "Indeterminate")) %>%
    dplyr::filter(!(Invasiveness == "NA"))


sel_cols <- c("Source", "Sample_ID", "Cluster_Group", "Invasiveness", "Driver_Group")
ss88 <- ss88[, sel_cols]
ss411 <- ss411[, sel_cols]
ss499 <- rbind(ss88, ss411)
# write_xlsx(ss499, "~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_499.xlsx")

betas88 <- betas88[,match(ss88$IDAT, colnames(betas88))] #391466     88
betas411 <- betas496[,match(ss411$Sample_ID, colnames(betas496))] #391466    411
colnames(betas88) <- ss88$Sample_ID

betas499 <- cbind(betas88, betas411)
# saveRDS(betas499, here("data", "20250203_betas499_unprocessed.rds"))

## PREPROCESSING MERGED BETAS________________________________________________________________________
betas499 <- readRDS(here("data", "20250203_betas499_unprocessed.rds"))
betas499 <- cleanMatrixForClusterW(betas499)
# Filter rows with >0.50 missingness and columns with >0.50 missingness.
# Before:  391466 rows and  499 columns.
# After:  341275 rows and  499 columns.
betas499 <- impute.knn(betas499, k = 10, rng.seed=1234)$data

# saveRDS(betas499, here("data", "20250203_betas499_processed.rds"))

## COMPUTE LEUKOCYTE FRACTION _______________________________________________________________
betas499 <- readRDS(here("data", "20250203_betas499_processed.rds"))
cmi_model <- readRDS(here("prev_data", "LeukoFrac_EPIC_20240614.rds"))
leuko = cmi_predict(betas499, cmi_model)

df <- as.data.frame(sapply(leuko, function(x) x))
rownames(df) <- colnames(betas499)
colnames(df) <- c("Leukocyte_Fraction")
saveRDS(df, here("data", "20250203_thyroid_thca_499_leuko.rds"))

## add to existing sheet
ss499 <- read_xlsx("~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_499.xlsx")
ss499$Leukocyte_Fraction = df$Leukocyte_Fraction
# write_xlsx(ss499, "~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_499.xlsx")

## MAKE TSNE ________________________________
ss499 <- read_xlsx("~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_499.xlsx")
betas499 <- readRDS(here("data", "20250203_betas499_processed.rds"))

mtx = bSubMostVariable(betas499, 3000) #341275    499
mtx = t(mtx) #499 3000

set.seed(12345)
tsne = Rtsne(mtx, dims=2, perplexity=22)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
rownames(df) = rownames(mtx)
saveRDS(df, here("data", "2025_thyroid_thca_499_tsne_coords.rds")) #499   2

# add to spreadsheet
ss499$tSNE1 <- df$tSNE1
ss499$tSNE2 <- df$tSNE2
write_xlsx(ss499, "~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_499.xlsx")

## PLOT TSNE ________________________________
ss499 <- read_xlsx("~/Documents/HPC_share/thyroid/ss/20250203_thyroid_thca_499.xlsx")

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
    "DICER1" = "#b47cffff"
)

pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_invasiveness2.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss499, aes(x = tSNE1, y = tSNE2, color = Cluster_Group, shape = Source)) +
    scale_color_manual(values = cluster_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_drivers2.pdf', width=6.5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss499, aes(x = tSNE1, y = tSNE2, color = Driver_Group, shape = Source)) +
    scale_color_manual(values = driver_colors) +
    scale_shape_manual(values = custom_shapes) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250202_merged_tsne_leuko2.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss499, aes(x = tSNE1, y = tSNE2, color = Leukocyte_Fraction, shape = Source)) +
    scale_shape_manual(values = custom_shapes) +
    scale_color_gradientn(colors=gnuplot(20)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", color = NA),  # Added color = NA
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

# pal.bands(brewer.prgn)

# DIFF METH ANALYSIS FOR JUST THCA _________________________________________________________________________

# PREPROCESSING THCA BETAS ______________________________
# ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx")
# thca = load(file = "/Users/jennyzli/Documents/HPC_share/thyroid/THCA.rda")
# betas <- get("betas")
#
# # rownames are weird, truncate them
# colnames(betas) <- substr(colnames(betas), 1, 15)
# betas = betas[, colnames(betas) %in% ss$Sample_ID] #485577    496
# saveRDS(betas, "~/Documents/HPC_share/thyroid/20240915_thca_betas.rds")

# PREPROCESSING _________________________________________________________________________
betas496 <- readRDS("~/Documents/HPC_share/thyroid/data/20240915_thca_betas.rds") #485577    496
betas496 <- cleanMatrixForClusterW(betas496)
# Filter rows with >0.50 missingness and columns with >0.50 missingness.
# Before:  485577 rows and  496 columns.
# After:  395802 rows and  496 columns.

betas496 <- t(impute.knn(betas496, k = 10, rng.seed=1234)$data)

rs <- grep('rs', colnames(betas496))
betas496 = betas496[,-rs] #496 395737
betas496 = t(betas496)
saveRDS(betas496, here("data", "20250204_thca_betas_processed.rds"))

# DIFF METHYLATION ON HPC_________________________________________________________________________
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx") %>% filter(!(Invasiveness == "NA"))
betas496 <- readRDS(here("data", "20250204_thca_betas_processed.rds"))
betas496 <- betas496[, colnames(betas496) %in% ss$Sample_ID] #395737    449

se = SummarizedExperiment(betas496, colData = ss)

se_ok = checkLevels(assay(se), colData(se)$Invasiveness)
colData(se)$Invasiveness <- relevel(factor(colData(se)$Invasiveness), "Low")
se = se[se_ok,] #395737    496
smry = DML(se, ~Invasiveness)
res = summaryExtractTest(smry)
head(res)

saveRDS(smry, here("diff_meth", "20250204_thca_smry_invasiveness.rds"))
saveRDS(res, here("diff_meth", "20250204_thca_res_invasiveness.rds"))

# VISUALIZATION ___________________________________
res <- readRDS(here("diff_meth", "20250204_thca_res_invasiveness.rds"))
res = testEnrichment(res1$Probe_ID[res1$Est_CLUSTER_GROUPGROUP2 > 0.2], platform="EPIC", universe=res1$Probe_ID)
min(res$Est_InvasivenessHigh)

query <- res$Probe_ID[res$Est_InvasivenessHigh < -0.05]
x<- testEnrichment(query, platform="HM450")
pdf(here("figures", "20250204_thca_invasiveness_all_enrichment_hyper.pdf"), width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(x)
dev.off()

# HIGH V LOW INVASIVE TFBS HYPO
a <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], "TFBS", platform="HM450", universe=res$Probe_ID)
pdf(here("figures", "20250204_thca_invasiveness_TFBS_enrichment_hyper.pdf"),  width = 4, height = 5, onefile = FALSE)
plotDot2(a, n_max = 20, n_min = 20)
dev.off()

# HIGH V LOW INVASIVE TFBS HYPER
b <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], "TFBS", platform="HM450", universe=res$Probe_ID)
pdf(here("figures", "20250204_thca_invasiveness_TFBS_enrichment_hypo.pdf"), width = 4, height = 5,  onefile = FALSE)
plotDot2(b, n_max = 20, n_min = 20)
dev.off()

# HIGH V LOW INVASIVE TISSUE HYPO
c <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], "tissueSignature", platform="HM450", universe=res$Probe_ID)
pdf(here("figures", "20250204_thca_invasiveness_tissue_enrichment_hyper.pdf"), width = 4, height = 5,  onefile = FALSE)
plotDot(c, n_max = 20, n_min = 20)
dev.off()

# HIGH V LOW INVASIVE TISSUE HYPER
d <- testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], "tissueSignature", platform="HM450", universe=res$Probe_ID)
pdf(here("figures", "20250204_thca_invasiveness_tissue_enrichment_hypo.pdf"), width = 4, height = 5,  onefile = FALSE)
plotDot(d, n_max = 20, n_min = 20)
dev.off()


# TSNE PLOTTING_________________________________________________________________________
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20241023_thca_master.xlsx") %>% filter(!(Invasiveness == "NA"))
betas <- readRDS(here("data", "20250204_thca_betas_processed.rds")) # 395737    496
betas449 <- betas[, colnames(betas) %in% ss$Sample_ID] #395737    449
betas449 <- betas449[,match(ss$Sample_ID, colnames(betas449))]

mtx = bSubMostVariable(betas449, 3000)
mtx = t(mtx) #449 3000

set.seed(12345)
tsne = Rtsne(mtx, dims=2, perplexity=22)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
rownames(df) = rownames(mtx)
# saveRDS(df, here("data", "20250204_thca496_tsne_coords.rds")) #449   2

df <- readRDS(here("data", "20250204_thca496_tsne_coords.rds")) #449   2

ss$tSNE1 <- df$tSNE1
ss$tSNE2 <- df$tSNE2

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

pdf('~/Documents/HPC_share/thyroid/figures/20250204_thca449_tsne_invasiveness.pdf', width=5, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = Invasiveness)) +
    scale_color_manual(values = invasiveness_colors) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250204_thca449_tsne_driver.pdf', width=5.5, height=5, onefile=FALSE)
ss_clean <- ss[!(ss$Driver_Group == "Indeterminate"), ]
p <- ggplot() +
    geom_point(data = ss_clean, aes(x = tSNE1, y = tSNE2, color = Driver_Group)) +
    scale_color_manual(values = driver_colors) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250204_thca449_tsne_differentiation.pdf', width=4.9, height=5, onefile=FALSE)
ss$Differentiation_Score <- as.numeric(ss$Differentiation_Score)
ss_clean <- ss[!is.na(ss$Differentiation_Score), ]
p <- ggplot() +
    geom_point(data = ss_clean, aes(x = tSNE1, y = tSNE2, color = Differentiation_Score)) +
    scale_color_gradientn(colors = brewer.bupu(20)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250204_thca449_tsne_brs.pdf', width=5, height=5, onefile=FALSE)
ss$ERK_Score <- as.numeric(ss$ERK_Score)
ss_clean <- ss[!is.na(ss$ERK_Score), ]
p <- ggplot() +
    geom_point(data = ss_clean, aes(x = tSNE1, y = tSNE2, color = ERK_Score)) +
    scale_color_gradientn(colors = brewer.bugn(20)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20250204_thca449_tsne_erk.pdf', width=4.9, height=5, onefile=FALSE)
ss$BRS_Score <- as.numeric(ss$BRS_Score)
ss_clean <- ss[!is.na(ss$BRS_Score), ]
p <- ggplot() +
    geom_point(data = ss_clean, aes(x = tSNE1, y = tSNE2, color = BRS_Score)) +
    scale_color_gradientn(colors = brewer.blues(20)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linewidth = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()






