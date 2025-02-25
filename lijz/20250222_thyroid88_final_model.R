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
bSubMostVariable <- function(betas, n=3000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

# ========================
# Data Loading and Preprocessing
# ========================
dir <- file.path(here(), "data")
date <- "20250223"

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

dir <- file.path(here(), "data", "20250222_mostvar_final")
date <- "20250223"

# FEATURE SELECTION
train_betas <- bSubMostVariable(betas)
var_probes <- rownames(train_betas)

# TRAIN MODEL
model <- randomForest(
    x = t(train_betas),
    y = as.factor(labels),
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)
# High Low class.error
# High   51   4  0.07272727
# Low     8  25  0.24242424

# IMPORTANCE SCORES
# importance_scores <- as.data.frame(importance(model))
# imp_results <- data.frame(
#     Feature = rownames(importance_scores),
#     Importance = importance_scores$MeanDecreaseAccuracy,
#     Gini = importance_scores$MeanDecreaseGini
# )

# PREDICTIONS
# pred_probs <- predict(model, t(betas[var_probes, ]), type="prob")
# pred_class <- predict(model, t(betas[var_probes, ]))
# prob don't wanna predict on trainin gdata??

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
#     # AUC = auc(roc_obj),
#     # Accuracy = mean(pred_class == labels),
#     Sensitivity = conf_matrix$byClass["Sensitivity"],
#     Specificity = conf_matrix$byClass["Specificity"],
#     PPV = conf_matrix$byClass["Pos Pred Value"],
#     NPV = conf_matrix$byClass["Neg Pred Value"],
#     F1 = conf_matrix$byClass["F1"]
# )

# OOB estimate of  error rate: 14.77%
# Confusion matrix:
#     High Low class.error
# High   51   4  0.07272727
# Low     9  24  0.27272727

saveRDS(model, file.path(dir, paste0(date, "_final_model_mostvar_3k.rds")))
# saveRDS(pred_results, file.path(dir, paste0(date, "_predictions_.rds")))
# write.csv(imp_results, file.path(dir, paste0(date, "_imp_values.csv")))

# ========================
# ROC Analysis
# ========================

# plot(roc_obj,
#      main = "ROC Curve",
#      col = "blue",
#      print.auc = TRUE)  # This will print AUC on the plot
#
# # plot(roc_obj,main ="ROC curve -- Logistic Regression ")
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
# pdf(file.path(fig_dir, paste0(date, "_thyroid88_roc_mostvar.pdf")), height = 6, width = 6)
# plot(p)
# dev.off()

# ========================
# SHAP Analysis
# ========================
train_data <- as.data.frame(t(train_betas))
train_labels <- factor(ifelse(as.factor(ss$Invasiveness) == "Low", 0, 1),
                       levels = c(0, 1)) # CRUCIAL STEP BRUH
model <- randomForest(
    x = train_data,
    y = train_labels,
    ntree = 500,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)

unified_model <- randomForest.unify(model, train_data)
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
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_beeswarm.pdf")), height = 8, width = 6)
plot(p)
dev.off()

# specify specific feature interaction
p <- sv_dependence(shp,
                   v = "cg09690283",
                   color_var = "cg01207474")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

# a little less correlated...
# p <- sv_dependence(shp,
#                    v = "cg00768487",
#                    color_var = "cg10579279")
# pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence2.pdf")), height = 3.5, width = 4)
# plot(p)
# dev.off()

# plots multiple feature interactions
imp_cpgs <- names((most_imp))
p <- sv_dependence(shp, v = imp_cpgs[1], color_var = imp_cpgs[2:10], interactions = TRUE)
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence3.pdf")), height = 10, width = 10)
plot(p)
dev.off()

# interaction plot
# p <- sv_interaction(shp, max_display = 2, size = 3) # max displays
p <- sv_interaction(shp, max_display = 20)
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_interaction.pdf")), height = 20, width = 20)
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







