# ========================
# Setup and Dependencies
# ========================

required_packages <- c(
    "tidyverse", "randomForest", "SummarizedExperiment", "sesame", "sail",
    "impute", "RSNNS", "e1071", "caret", "ISLR", "pROC", "readxl", "shapviz",
    "future", "future.apply", "logger", "treeshap", "here", "viridis"
)
load_packages <- function(packages) {
    for (package in packages) {
        if (!require(package, character.only = TRUE)) {
            stop(paste("Package", package, "is required but not installed"))
        }
    }
}
load_packages(required_packages)

# ========================
# Helper Functions
# ========================

# ========================
# Data Loading and Preprocessing
# ========================
dir <- file.path(here(), "data")
date <- "20250222"

# Load data
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")

# Load and preprocess betas
betas <- readRDS(here("data", "20250214_thyroid88_processed_condensed_betas.rds"))
betas <- betas[, colnames(betas) %in% ss$IDAT]
betas <- betas[,match(ss$IDAT, colnames(betas))]

# Prepare labels and parameters
labels <- as.factor(ss$Invasiveness)
n_trees <- 500

# ========================
# Make Final Model
# ========================
dir <- file.path(here(), "data", "20250223_subset_final")
# dir.create(dir, recursive = TRUE, showWarnings = FALSE)
date <- "20250223"

# FEATURE SELECTION
batch_size = 10000
n_probes <- nrow(betas)
n_batches <- ceiling(n_probes/batch_size)
fold_models <- list()
fold_importance <- data.frame()

for(s in 0:(n_batches-1)) {
    # Calculate indices for current batch
    start_idx <- s * batch_size + 1
    end_idx <- min((s + 1) * batch_size, n_probes)

    sbetas <- betas[start_idx:end_idx, ]
    # print(paste("Batch dimensions:", dim(sbetas)[1], "x", dim(sbetas)[2]))

    model <- randomForest(
        x = t(sbetas),
        y = as.factor(labels),
        ntree = n_trees,
        importance = TRUE
    )
    fold_models[[s+1]] <- model

    importance <- data.frame(
        Model = s,
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[,"MeanDecreaseGini"]
    )
    fold_importance <- rbind(fold_importance, importance)
}
features <- fold_importance %>%
    arrange(desc(MeanDecreaseAccuracy))

saveRDS(fold_models, file.path(dir, paste0(date, "_fold_models.rds")))
saveRDS(features, file.path(dir, paste0(date, "_fold_importance.rds")))

##
most_imp <- fold_importance %>%
    arrange(desc(MeanDecreaseAccuracy)) %>%
    head(3000)

sel_probes <- rownames(most_imp)
train_betas <- betas[sel_probes, ]

# TRAIN MODEL
set.seed(123)
model <- randomForest(
    x = t(train_betas),
    y = labels,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)

# OOB estimate of  error rate: 10.23%
# Confusion matrix:
#     High Low class.error
# High   53   2  0.03636364
# Low     7  26  0.21212121

# IMPORTANCE SCORES
# importance_scores <- as.data.frame(importance(model))
# imp_results <- data.frame(
#     Feature = rownames(importance_scores),
#     MeanDecreaseAccuracy = importance_scores$MeanDecreaseAccuracy,
#     MeanDecreaseGini = importance_scores$MeanDecreaseGini
# )

# PREDICTIONS
# pred_probs <- predict(model, t(betas[var_probes, ]), type="prob")
# pred_class <- predict(model, t(betas[var_probes, ]))
#
# pred_results <- data.frame(
#     True_Label = labels,
#     Predicted = pred_class,
#     Probability = pred_probs[,2],
#     Sample_ID = ss$Sample_ID
# )

# ROC and metrixs
# roc_obj <- roc(labels, pred_probs[,2])
# conf_matrix <- caret::confusionMatrix(as.factor(pred_class), as.factor(labels))

# Reference
# Prediction High Low
# High   53   0
# Low     2  33

# SUMMARIZED RESULTS
# results <- data.frame(
#     AUC = auc(roc_obj),
#     Accuracy = mean(pred_class == labels),
#     Sensitivity = conf_matrix$byClass["Sensitivity"],
#     Specificity = conf_matrix$byClass["Specificity"],
#     PPV = conf_matrix$byClass["Pos Pred Value"],
#     NPV = conf_matrix$byClass["Neg Pred Value"],
#     F1 = conf_matrix$byClass["F1"]
# )

saveRDS(model, file.path(dir, paste0(date, "_final_subset_model_1k.rds")))
# saveRDS(pred_results, file.path(dir, paste0(date, "_predictions_.rds")))
# write.csv(imp_results, file.path(dir, paste0(date, "_imp_values.csv")))
model <- readRDS(file.path(dir, paste0(date, "_final_subset_model_1k.rds")))

# ========================
# ROC Analysis
# ========================
# roc_df <- data.frame(
#     FPR = 1 - roc_obj$specificities,
#     TPR = roc_obj$sensitivities
# )
#
# p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
#     geom_line(color = "#2c7fb8", size = 1.2) +
#     geom_abline(linetype = "dashed", alpha = 0.5) +  # diagonal reference line
#     theme_minimal() +
#     labs(
#         title = "ROC Curve",
#         subtitle = paste("AUC =", round(auc(roc_obj), 3)),
#         x = "False Positive Rate (1 - Specificity)",
#         y = "True Positive Rate (Sensitivity)"
#     ) +
#     theme(
#         plot.title = element_text(size = 14, face = "bold"),
#         plot.subtitle = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12)
#     ) +
#     coord_equal()  # force aspect ratio 1:1
# pdf(file.path(fig_dir, paste0(date, "_thyroid88_roc_mostvar.pdf")), height = 3.5, width = 4)
# plot(p)
# dev.off()
# ========================
# SHAP Analysis
# ========================
train_data <- as.data.frame(t(train_betas))
train_labels <- factor(ifelse(as.factor(ss$Invasiveness) == "Low", 0, 1),
                       levels = c(0, 1)) # CRUCIAL STEP BRUH
model_shp <- randomForest(
    x = train_data,
    y = train_labels,
    ntree = 500,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
) #it's somehow diff than the original model?? check this

unified_model <- randomForest.unify(model_shp, train_data)
shap <- treeshap(unified_model, train_data, interactions = TRUE)
shp <- shapviz(shap, X = train_data)

# most important features
most_imp <- sort(colMeans(abs(shap$shaps)), decreasing = TRUE)
print(head(most_imp))

## PLOTTING
p <- sv_importance(shp,
                   kind = "beeswarm",
                   max_display = 30)
# p <- sv_importance(shp, kind = "both") # if u want both beeswarm and barplot
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_beeswarm_subset.pdf")), height = 8, width = 6)
plot(p)
dev.off()

# specify specific feature interaction
p <- sv_dependence(shp,
                   v = "cg00768487",
                   color_var = "cg11868461")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence_subset.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

# a little less correlated...
p <- sv_dependence(shp,
                   v = "cg00768487",
                   color_var = "cg10579279")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence2_subset.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

# plots multiple feature interactions
imp_cpgs <- names((most_imp))
p <- sv_dependence(shp, v = imp_cpgs[1], color_var = imp_cpgs[2:10], interactions = TRUE)
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap2_dependence3.pdf")), height = 10, width = 10)
plot(p)
dev.off()

# interaction plot
# p <- sv_interaction(shp, max_display = 2, size = 3) # max displays
p <- sv_interaction(shp, max_display = 20)
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_interaction_subset.pdf")), height = 20, width = 20)
plot(p)
dev.off()

# force plot 1
p <- sv_force(shp, row_id = "207222790148_R07C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_force1.pdf")), height = 4, width = 7)
plot(p)
dev.off()
p <- sv_force(shp, row_id = "207700160043_R03C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_force2.pdf")), height = 4, width = 7)
plot(p)
dev.off()

# waterfall plot
p <- sv_waterfall(shp, row_id = "207222790148_R07C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_waterfall1.pdf")), height = 7, width = 7)
plot(p)
dev.off()
p <- sv_waterfall(shp, row_id = "207700160043_R03C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_waterfall2.pdf")), height = 7, width = 7)
plot(p)
dev.off()

# ========================
# Enrichment Analysis
# ========================
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
se = se[sel_probes,]

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Invasiveness <- relevel(factor(colData(se)$Invasiveness), "Low")

# se = se[1:10,]

smry = DML(se, ~Invasiveness + Sex + Actual_Age)
res = summaryExtractTest(smry)
head(res)

saveRDS(smry, here("diff_meth", "20250223_thyroid88_smry_invasiveness_3k_condensed.rds"))
saveRDS(res, here("diff_meth", "20250223_thyroid88_res_invasiveness_3k_condensed.rds"))


# res = readRDS(here("data", "20250214_thyroid88_res_invasiveness_condensed.rds"))

# Low v High TFBS HYPO
pdf(here("figures", "20250223_thyroid88_hyper_rfmodel_invasiveness_condensed.pdf"), width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh > 0.2], "TFBSconsensus", platform="EPIC", universe=res$Probe_ID), n_min = 40)
dev.off()

pdf(here("figures", "20250223_thyroid88_hypo_rfmodel_invasiveness_condensed.pdf"), width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], "TFBSconsensus", platform="EPIC", universe=res$Probe_ID), n_min = 40)
dev.off()

pdf(here("figures", "20250223_thyroid88_hypo_rfmodel_invasiveness_tissue_condensed.pdf"), width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res$Probe_ID[res$Est_InvasivenessHigh < -0.2], "tissueSignature", platform="EPIC", universe=res$Probe_ID), n_min = 40)
dev.off()







