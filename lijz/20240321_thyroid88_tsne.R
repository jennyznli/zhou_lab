library(sesame)
library(Rtsne)
library(ggplot2)
library(plotly)
ss = read_excel("~/Documents/HPC_share/thyroid/20231102_thyroid_master.xlsx")
ss88 = ss %>% dplyr::filter(INCLUDE_IN_ANALYSIS == "1")
ss1 = ss %>% dplyr::filter(Histology == "PTmC") #LN metastasis

##
betas = readRDS(file = "~/Documents/HPC_share/thyroid/20240320_thyroid136_betas_condensed.rds") #930659
# betas1 = readRDS(file = "~/Documents/HPC_share/thyroid/20240320_thyroid136_betas.rds") #937690

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

rs <- grep('rs', rownames(betas_knn))
betas1 = betas_knn[-b,] # 856541    136

# se = SummarizedExperiment(betas1, colData=ss)

bSubMostVariable <- function(betas, n=3000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

mtx = bSubMostVariable(betas1, 3000)
mtx = t(mtx)
mtx88 = mtx[rownames(mtx) %in% ss88$IDAT,] ## 88 3000

## try most important probes instead!!_______________________________________________________________

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
# betas <- cleanMatrixForClusterW(betas) #855437 rows and  88 columns
#
# # knn imputation
# betas1 <- impute.knn(betas, k = 10, rng.seed=1234)
# betas_knn <- betas1$data
#
# ##
#
# imp = read_tsv(col_names = FALSE, "~/Documents/HPC_share/thyroid/20240610_features.tsv")
# colnames(imp) = c("Model", "CpG", "Accuracy", "Gini")
# imp = imp %>% arrange(desc(Accuracy))
#
# pos = imp %>% filter(Accuracy > 0)
# rs <- grep('rs', pos$CpG)
# # print(length(rs))
# mostImp = pos[1:3000,]$CpG
#
# betas_knn = t(betas_knn)
# mtx88_mi = betas_knn[rownames(betas_knn) %in% ss88$IDAT, colnames(betas_knn) %in% mostImp]

# TSNE _________________________________________________________________________

set.seed(123)
tsne = Rtsne(mtx88, dims=2, perplexity=9.4)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")

ggplotly(ggplot(df) + geom_point(aes(tSNE1, tSNE2)))

# INCLUDE LOW CONFIDENCE  _______________________________________________________________

ss88$tSNE1 <- df88$tSNE1
ss88$tSNE2 = df88$tSNE2

low_conf = subset(ss88, ss88$CONFIDENCE == "Low")

#  invasiveness as pdf
pdf('~/Documents/HPC_share/thyroid/20240703_tsne88_invasiveness.pdf', family="ArialMT", width=7, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss88, aes(x = tSNE1, y = tSNE2, color = Predicted_Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Predicted Invasiveness:", Predicted_Invasiveness))) +
    geom_point(data = low_conf, aes(x = tSNE1, y = tSNE2, color = Predicted_Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Predicted Invasiveness:", Predicted_Invasiveness)), size=2, shape = 18) +
    ggtitle("Invasiveness tSNE") +
    theme_bw()
plot(p)
dev.off()

#  predicted invasiveness as plotly
p <- ggplot() +
    geom_point(data = ss88, aes(x = tSNE1, y = tSNE2, color = Predicted_Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Predicted Invasiveness:", Predicted_Invasiveness))) +
    geom_point(data = low_conf, aes(x = tSNE1, y = tSNE2, color = Predicted_Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Predicted Invasiveness:", Predicted_Invasiveness)), size=2, shape = 18) +
    ggtitle("88 Thyroid Predicted Invasiveness tSNE") +
    theme_bw()
pp <- ggplotly(p, tooltip = c("text"))
pp

# invasiveness
p <- ggplot() +
    geom_point(data = ss88, aes(x = tSNE1, y = tSNE2, color = Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness))) +
    geom_point(data = low_conf, aes(x = tSNE1, y = tSNE2, color = Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness)), size=2, shape = 18) +
    ggtitle("88 Thyroid Invasiveness tSNE") +
    theme_bw()
pp <- ggplotly(p, tooltip = c("text"))
pp

## PTC V FTC

p <- ggplot() +
    geom_point(data = ss88, aes(x = tSNE1, y = tSNE2, color = Histology, text = paste("Sample:", Sample_ID, "\n", "Histology:", Histology))) +
    geom_point(data = low_conf, aes(x = tSNE1, y = tSNE2, color = Histology, text = paste("Sample:", Sample_ID, "\n", "Histology:", Histology)), size=2, shape = 18) +
    ggtitle("88 Thyroid Histology tSNE") +
    theme_bw()
pp <- ggplotly(p, tooltip = c("text"))
pp

p <- ggplot() +
    geom_point(data = ss88, aes(x = tSNE1, y = tSNE2, color = Driver_Group, text = paste("Sample:", Sample_ID, "\n", "Histology:", Driver_Group))) +
    geom_point(data = low_conf, aes(x = tSNE1, y = tSNE2, color = Driver_Group, text = paste("Sample:", Sample_ID, "\n", "Histology:", Driver_Group)), size=2, shape = 18) +
    ggtitle("88 Thyroid Driver Group tSNE") +
    theme_bw()
pp <- ggplotly(p, tooltip = c("text"))
pp








x <- ggplot() +
    geom_point(data = ss88, aes(x = tSNE1, y = tSNE2, color = Invasiveness), text = paste("Sample:", ss88$Sample_ID, "\n", "Invasiveness:", ss88$Invasiveness)) +
    geom_point(data = low_conf, aes(x = tSNE1, y = tSNE2, color = Invasiveness), text = paste("Sample:", ss88$Sample_ID, "\n", "Invasiveness:", ss88$Invasiveness)) +
    ggtitle("88 Thyroid Invasiveness tSNE") +
    theme_bw()
ggplotly(x)
ggplotly(ggplot(df) + geom_point(aes(tSNE1, tSNE2)))


low_conf136 = subset(ss136_10k, ss136_10k$`Low-confidence` == "Low-confidence")

##PLOTLY
x <- ggplot() +
    geom_point(data = ss136_10k, aes(x = tSNE1, y = tSNE2, color = Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness))) +
    geom_point(data = low_conf136, aes(x = tSNE1, y = tSNE2, color = Invasiveness, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness)), size=2, shape = 18) +
    ggtitle("136 Thyroid Invasiveness tSNE") +
    theme_bw()
ggplotly(x)


#THIS DOESN'T WORK??
ggsave(
    filename = paste0("20240321_thyroid90_", "Invasiveness", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

# LYMPH NODE RATIO  _______________________________________________________________

na_LN = subset(ss90_10k, is.na(ss90_10k$LNP_Decimal))
LNP = subset(ss90_10k, !is.na(ss90_10k$LNP_Decimal))

x <- ggplot() +
    geom_point(data = na_LN, aes(x = tSNE1, y = tSNE2, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness)), shape = 1) +
    geom_point(data = LNP, aes(x = tSNE1, y = tSNE2, color = LNP_Decimal, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness)), shape = 16) +
    scale_color_gradient(low = "blue", high = "orange") +
    ggtitle("90 Thyroid LN Positive Ratio tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid90_", "LNP", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

na_LN136 = subset(ss136_10k, is.na(ss136_10k$LNP_Decimal))
LNP136 = subset(ss136_10k, !is.na(ss136_10k$LNP_Decimal))

x <- ggplot() +
    geom_point(data = na_LN136, aes(x = tSNE1, y = tSNE2, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness)), shape = 1) +
    geom_point(data = LNP136, aes(x = tSNE1, y = tSNE2, color = LNP_Decimal, text = paste("Sample:", Sample_ID, "\n", "Invasiveness:", Invasiveness)), shape = 16) +
    scale_color_gradient(low = "blue", high = "orange") +
    ggtitle("136 Thyroid LN Positive Ratio tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid136_", "LNP", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)



# 10K TSNES  _______________________________________________________________

ss136_10k = ss
ss90_10k = ss
ss90_10k = ss90_10k %>% dplyr::filter(IncludeInAnalysis == "1")

ss136_10k$tSNE1 = df10k$tSNE1
ss136_10k$tSNE2 = df10k$tSNE2

ss90_10k$tSNE1 = df90_10k$tSNE1
ss90_10k$tSNE2 = df90_10k$tSNE2

write_csv(ss136_10k,"~/Documents/HPC_share/thyroid/20240326_thyroid136_tSNE_10k.csv")
write_csv(ss90_10k,"~/Documents/HPC_share/thyroid/20240326_thyroid90_tSNE_10k.csv")

ggplotly(ggplot(ss136_10k) + geom_point(aes(tSNE1, tSNE2, color=Sex, text = (paste("Sample:", Sample_ID, "Invasiveness:", Invasiveness)))))
ggplotly(ggplot(ss90_10k) + geom_point(aes(tSNE1, tSNE2, color=Sex, text = (paste("Sample:", Sample_ID, "Invasiveness:", Invasiveness)))))

# TSNE COMPARISONS  _______________________________________________________________

setwd("/Users/jennyzli/Documents/HPC_share/thyroid/thyroid136_plots/")
variables = c(
    "Invasiveness",
    "Sex",
    "Invasiveness",
    "T", "N", "M"
)

cont_variables = c(
    "ProbeSuccessRate",
    "Age at Surgery"
)

# CONSTANT SIZE
for (n in variables) {
    x = ggplot() +
        geom_point(data = ss90_10k, aes(x = tSNE1, y = tSNE2, color = get(n)), size=.75) +
        ggtitle(paste0("90 Thyroid ", n, " tSNE")) +
        theme_bw()
    # print(x)
    ggsave(
        filename = paste0("20240321_thyroid90_", n, ".png"),
        plot = x,
        device = NULL,
        path = NULL,
        scale = 1,
        width = 5,
        height = 4,
        units = "in",
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
    )
}


## PROBE SUCCESS RATE
x <- ggplot(ss90_10k, aes(tSNE1, tSNE2))+
    geom_point(aes(color = ProbeSuccessRate), size=.75) +
    scale_color_gradient(low = "blue", high = "red") +
    ggtitle("90 Thyroid ProbeSuccessRate tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid90_", "ProbeSuccessRate", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## AGE - DIDN'T WORK
# x <- ggplot(ss90_10k, aes(tSNE1, tSNE2))+
#     geom_point(aes(color = get("Age at Surgery")), size=.75) +
#     scale_color_gradient(low = "blue", high = "red") +
#     ggtitle("90 Thyroid Age at Surgery tSNE") +
#     theme_bw()
# ggsave(
#     filename = paste0("20240321_thyroid136_", "Age at Surgery", ".png"),
#     plot = x,
#     device = NULL,
#     path = NULL,
#     scale = 1,
#     width = 5,
#     height = 4,
#     units = "in",
#     dpi = 300,
#     limitsize = TRUE,
#     bg = NULL,
# )

## DRIVER
x = ggplot() +
    geom_point(data = ss90_10k, aes(x = tSNE1, y = tSNE2, color = Driver), size=.75) +
    ggtitle("90 Thyroid Driver tSNE") +
    theme_bw()
ggplotly(x)
ggplotly(ggplot(ss90_10k) + geom_point(aes(tSNE1, tSNE2, color=Driver, text = (paste("Sample:", Sample_ID, "Invasiveness:", Invasiveness)))))

ggsave(
    filename = paste0("20240321_thyroid90_", "Driver", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 5,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## DRIVER
x = ggplot() +
    geom_point(data = ss90_10k, aes(x = tSNE1, y = tSNE2, color = Driver_Group), size=.75) +
    ggtitle("90 Thyroid Driver_Group tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid90_", "Driver Group", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## HISTOLOGY
x = ggplot() +
    geom_point(data = ss90_10k, aes(x = tSNE1, y = tSNE2, color = Histology), size=.75) +
    ggtitle("90 Thyroid Histology tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid90_", "Histology", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## PTC VARIANT
x = ggplot() +
    geom_point(data = ss90_10k, aes(x = tSNE1, y = tSNE2, color = get("PTC variant")), size=.75) +
    ggtitle("90 Thyroid PTC variant tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid90_", "Variant", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)


## REPEAT WITH 136

# CONSTANT SIZE
for (n in variables) {
    x = ggplot() +
        geom_point(data = ss136_10k, aes(x = tSNE1, y = tSNE2, color = get(n)), size=.75) +
        ggtitle(paste0("136 Thyroid ", n, " tSNE")) +
        theme_bw()
    # print(x)
    ggsave(
        filename = paste0("20240321_thyroid136_", n, ".png"),
        plot = x,
        device = NULL,
        path = NULL,
        scale = 1,
        width = 5,
        height = 4,
        units = "in",
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
    )
}

## PROBE SUCCESS RATE
x <- ggplot(ss136_10k, aes(tSNE1, tSNE2))+
    geom_point(aes(color = ProbeSuccessRate), size=.75) +
    scale_color_gradient(low = "blue", high = "red") +
    ggtitle("136 Thyroid ProbeSuccessRate tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid136_", "ProbeSuccessRate", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## AGE - DIDN'T WORK
# x <- ggplot(ss90_10k, aes(tSNE1, tSNE2))+
#     geom_point(aes(color = get("Age at Surgery")), size=.75) +
#     scale_color_gradient(low = "blue", high = "red") +
#     ggtitle("90 Thyroid Age at Surgery tSNE") +
#     theme_bw()
# ggsave(
#     filename = paste0("20240321_thyroid136_", "Age at Surgery", ".png"),
#     plot = x,
#     device = NULL,
#     path = NULL,
#     scale = 1,
#     width = 5,
#     height = 4,
#     units = "in",
#     dpi = 300,
#     limitsize = TRUE,
#     bg = NULL,
# )

## DRIVER
x = ggplot() +
    geom_point(data = ss136_10k, aes(x = tSNE1, y = tSNE2, color = Driver), size=.75) +
    ggtitle("136 Thyroid Driver tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid136_", "Driver", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## DRIVER
x = ggplot() +
    geom_point(data = ss136_10k, aes(x = tSNE1, y = tSNE2, color = Driver_Group), size=.75) +
    ggtitle("136 Thyroid Driver_Group tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid136_", "Driver Group", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## HISTOLOGY
x = ggplot() +
    geom_point(data = ss136_10k, aes(x = tSNE1, y = tSNE2, color = Histology), size=.75) +
    ggtitle("136 Thyroid Histology tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid136_", "Histology", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)

## PTC VARIANT
x = ggplot() +
    geom_point(data = ss136_10k, aes(x = tSNE1, y = tSNE2, color = get("PTC variant")), size=.75) +
    ggtitle("136 Thyroid PTC variant tSNE") +
    theme_bw()
ggsave(
    filename = paste0("20240321_thyroid136_", "Variant", ".png"),
    plot = x,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height = 4,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
)
# SEPARATE CLUSTERS _______________________________________________________________

ggplotly(ggplot(df) + geom_point(aes(tSNE1, tSNE2)))

ss_filter <- filter(ss, tSNE2 > 0)
ss_filter <- filter(ss_filter, tSNE1 > 8)

ss_filter2 <- filter(ss, tSNE2 < 0)
ss_filter2 <- filter(ss_filter2, tSNE1 < 8)
write.xlsx(ss_filter, "~/Downloads/cluster1.xlsx", showNA=FALSE)
write.xlsx(ss_filter2, "~/Downloads/cluster2.xlsx",showNA=FALSE)

for (n in 1:88) {
    if (ss$IDAT[n] %in% ss_filter$IDAT) {
        ss$Cohort[n] = "1"
    } else {
        ss$Cohort[n] = "2"
    }
}

# RESTRICT RS PROBES _______________________________________________________________

b <- grep('rs', rownames(betas))
betas_rs = betas[b,]

set.seed(123)
betas_rs = t(betas_rs)
tsne_rs = Rtsne(betas_rs, dims=2, perplexity=9)

df = as.data.frame(tsne_rs$Y)
colnames(df) = c("tSNE1", "tSNE2")
df$sample = rownames(betas_rs)

ss <- read_excel("~/Downloads/20230509_EPICv2_CHOP_THYROID_Julio_NEW.xlsx")

write_csv(ss,"~/Documents/HPC_share/thyroid/20231006_tSNE_p9.csv")
ss$
ggplotly(ggplot(df) + geom_point(aes(ss$tsne1, ss$tsne2)))

ggplotly(ggplot(ss) + geom_point(aes(ss$tsne1, ss$tsne2, color=Sex, text = (paste("Sample:", ss$Sample_Name...2)))))

# EXAMPLE _______________________________________________________________
se = sesameDataGet("MM285.10.SE.tissue")[1:1000,] # an arbitrary 1000 CpGs
cd = as.data.frame(colData(se)); rownames(cd) = NULL
cd

se_ok = (checkLevels(assay(se), colData(se)$sex) &
             checkLevels(assay(se), colData(se)$tissue))
sum(betas_ok)   # the number of CpGs that passes

se = se[se_ok,]

colData(se)$tissue <- relevel(factor(colData(se)$tissue), "Colon")
colData(se)$sex <- relevel(factor(colData(se)$sex), "Female")

smry = DML(se, ~tissue + sex)
smry
# MODELING _______________________________________________________________

betas_ok = (checkLevels(betas, ss$Cohort))
sum(betas_ok)   # the number of CpGs that passes

betas1 = betas[betas_ok,]

ss$Cohort <- relevel(factor(ss$Cohort), "1")

smry = DML(betas1, ~ss$Cohort)
smry







