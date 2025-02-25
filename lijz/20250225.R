# Script: Thyroid Cancer Methylation Analysis Pipeline
# Description: Complete analysis pipeline for thyroid cancer methylation data
# Last Updated: 2025-01-30

required_packages <- c(
    "tidyverse", "randomForest", "SummarizedExperiment", "sesame", "sail",
    "impute", "RSNNS", "e1071", "caret", "ISLR", "pROC", "readxl", "shapviz",
    "future", "future.apply", "logger"
)
load_packages <- function(packages) {
    for (package in packages) {
        if (!require(package, character.only = TRUE)) {
            stop(paste("Package", package, "is required but not installed"))
        }
    }
}

# PREPROCESSING ON HPC_______________________________________________________________
packages <- c("tidyverse", "randomForest", "SummarizedExperiment", "sesame", "sail",
              "impute", "RSNNS", "e1071", "caret", "ISLR", "pROC")
lapply(packages, require, character.only = TRUE)

meta = read_excel("/home/lijz/20230818_thyroid_cancer/20231102_thyroid_master.xlsx")
betas = readRDS("/home/lijz/20230818_thyroid_cancer/20240320_thyroid136_betas_condensed.rds")
meta = meta %>% dplyr::filter(Include_In_Analysis == "1")
print(dim(meta))

print(dim(betas))

se = SummarizedExperiment(betas, colData=meta)
print("loaded")

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
imputeRowMean <- function(mtx) {
    k <- which(is.na(mtx), arr.ind=TRUE)
    mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
    mtx
}
mtx <- cleanMatrixForClusterW(assay(se)) %>%
    imputeRowMean(.)

 print("cleaned")

labels <- as.factor(meta$Invasiveness)
print(labels)
saveRDS(mtx,"/home/lijz/20230818_thyroid_cancer/20250125_thyroid88_betas_condensed_processed.rds")

# betas <- readRDS("~/Documents/HPC_share/thyroid/data/20250125_thyroid88_betas_condensed_processed.rds")
# 855437     88

# MAKE FOLDS ON MY COMPUTER_______________________________________________________________
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")
ss = ss %>% dplyr::filter(Include_In_Analysis == "1")

# Create folds
set.seed(42)
folds <- createFolds(ss$Invasiveness, k=10)
saveRDS(folds, "~/Documents/HPC_share/thyroid/data/20250125_thyroid_folds.rds")

# TRY JUST MOST VARIABLE 3K SELECTION_______________________________________________________________
# ========================
# Setup and Dependencies
# ========================
library(tidyverse)
library(randomForest)
library(pROC)
library(caret)
library(readxl)
library(here)
library(shapr)

# ========================
# Helper Functions
# ========================
bSubMostVariable <- function(betas, n=3000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

save_fold_results <- function(model, pred_results, shap_results, fold, model_dir, feat_dir) {
    saveRDS(model, file.path(model_dir, paste0("20250210_final_model_fold", fold, ".rds")))
    saveRDS(pred_results, file.path(model_dir, paste0("20250210_predictions_fold", fold, ".rds")))
    write.csv(shap_results, file.path(feat_dir, paste0("20250210_shap_values_fold", fold, ".csv")))
}

# ========================
# Data Loading and Preprocessing
# ========================
# Set up directories
model_dir <- file.path(here(), "models")
feat_dir <- file.path(model_dir, "features")
dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
folds <- readRDS(file.path(here(), "data", "20250125_thyroid_folds.rds"))
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")

# Load and preprocess betas
betas <- readRDS(file.path(here(), "data", "20250125_thyroid88_betas_condensed_processed.rds"))
betas <- betas[, colnames(betas) %in% ss$IDAT]
betas <- betas[,match(ss$IDAT, colnames(betas))]

# Prepare labels and parameters
labels <- as.factor(ss$Invasiveness)
n_trees <- 500

# Initialize results containers
models_list <- list()
pred_results <- list()
results <- data.frame()
roc_curves <- list()
shap_results <- list()
feature_lists <- list()


# ========================
# Cross-Validation Loop
# ========================

for(i in 1:10) {
    print(paste("Processing fold", i))
    train_idx <- unlist(folds[-i])
    valid_idx <- folds[[i]]

    # Select 3k most variable features
    bSubMostVariable <- function(betas, n=3000) {
        std <- apply(betas, 1, sd, na.rm=TRUE)
        betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
    }
    train_betas <- bSubMostVariable(betas[,train_idx])
    var_probes <- rownames(train_betas)
    feature_lists[[i]] <- var_probes  # Store selected features

    # Train model
    model <- randomForest(
        x = t(train_betas),
        y = as.factor(labels[train_idx]),
        ntree = n_trees,
        importance = TRUE  # Enable importance calculation
    )
    models_list[[i]] <- model

    # Calculate SHAP values for this fold
    importance_scores <- as.data.frame(importance(model))
    top_features <- rownames(importance_scores[order(importance_scores$MeanDecreaseAccuracy, decreasing = TRUE),])[1:20]

    # Store SHAP results for this fold
    shap_results[[i]] <- data.frame(
        Feature = rownames(importance_scores),
        Importance = importance_scores$MeanDecreaseAccuracy,
        Fold = i
    )

    # Predictions
    pred_probs <- predict(model, t(betas[var_probes, valid_idx]), type="prob")
    preds <- predict(model, t(betas[var_probes, valid_idx]))

    pred_results[[i]] <- data.frame(
        Fold = i,
        True_Label = labels[valid_idx],
        Predicted = preds,
        Probability = pred_probs[,2],
        Sample_ID = valid_idx
    )

    # ROC and metrics
    roc_obj <- roc(labels[valid_idx], pred_probs[,2])
    roc_curves[[i]] <- roc_obj
    conf_matrix <- caret::confusionMatrix(as.factor(preds), as.factor(labels[valid_idx]))

    results <- rbind(results, data.frame(
        Fold = i,
        AUC = auc(roc_obj),
        Accuracy = mean(preds == labels[valid_idx]),
        Sensitivity = conf_matrix$byClass["Sensitivity"],
        Specificity = conf_matrix$byClass["Specificity"],
        PPV = conf_matrix$byClass["Pos Pred Value"],
        NPV = conf_matrix$byClass["Neg Pred Value"]
    ))

    # Save fold-specific results
    saveRDS(model, file.path(model_dir, paste0("20250210_final_model_fold", i, ".rds")))
    saveRDS(pred_results[[i]], file.path(model_dir, paste0("20250210_predictions_fold", i, ".rds")))
    write.csv(shap_results[[i]], file.path(feat_dir, paste0("20250210_shap_values_fold", i, ".csv")))

    print(paste("Completed fold", i))
}

# ========================
# Process Feature Importance
# ========================
all_shap <- bind_rows(shap_results)
avg_importance <- all_shap %>%
    group_by(Feature) %>%
    summarize(
        Mean_Importance = mean(Importance),
        SD_Importance = sd(Importance),
        Times_Selected = n()
    ) %>%
    arrange(desc(Mean_Importance))

# Create feature overlap analysis
feature_overlap <- data.frame(
    Feature = unique(unlist(feature_lists)),
    Times_Selected = sapply(unique(unlist(feature_lists)),
                            function(x) sum(sapply(feature_lists, function(y) x %in% y)))
)

# ========================
# Save Final Results
# ========================
# Save model results
saveRDS(models_list, file.path(model_dir, "20250214_final_model_list.rds"))
saveRDS(pred_results, file.path(model_dir, "20250214_pred_results.rds"))
saveRDS(results, file.path(model_dir, "20250214_results_df.rds"))
saveRDS(roc_curves, file.path(model_dir, "20250214_roc_curves.rds"))

# Save feature importance results
write.csv(avg_importance, file.path(feat_dir, "20250214_average_importance.csv"))
write.csv(feature_overlap, file.path(feat_dir, "20250214_feature_overlap.csv"))
saveRDS(shap_results, file.path(feat_dir, "20250214_shap_results_list.rds"))
saveRDS(feature_lists, file.path(feat_dir, "20250214_feature_lists.rds"))

# ========================
# Print Summary Statistics
# ========================
print("SHAP Analysis Summary:")
print(paste("Total unique features:", nrow(feature_overlap)))
print(paste("Features selected in all folds:", sum(feature_overlap$Times_Selected == 10)))
print("\nTop 10 most important features by mean importance:")
print(head(avg_importance, 10))


# MODEL ANALYSIS______________________________________________________________
# Set up directories
model_dir <- file.path(here(), "models")
feat_dir <- file.path(model_dir, "features")
fig_dir <- file.path(here(), "figures")
dir.create(fig_dir, showWarnings = FALSE)

# Load all required data
results <- readRDS(file.path(model_dir, "20250214_results_df.rds"))
pred_results <- readRDS(file.path(model_dir, "20250214_pred_results.rds"))
roc_curves <- readRDS(file.path(model_dir, "20250214_roc_curves.rds"))
models_list <- readRDS(file.path(model_dir, "20250214_final_model_list.rds"))
shap_results <- readRDS(file.path(feat_dir, "20250214_shap_results_list.rds"))

# ========================
# Performance Analysis
# ========================

# Calculate overall performance metrics
performance_summary <- results %>%
    summarize(
        Mean_AUC = mean(AUC),
        SD_AUC = sd(AUC),
        Mean_Accuracy = mean(Accuracy),
        SD_Accuracy = sd(Accuracy),
        Mean_Sensitivity = mean(Sensitivity),
        SD_Sensitivity = sd(Sensitivity),
        Mean_Specificity = mean(Specificity),
        SD_Specificity = sd(Specificity),
        Mean_PPV = mean(PPV),
        SD_PPV = sd(PPV),
        Mean_NPV = mean(NPV),
        SD_NPV = sd(NPV)
    )

# Print performance metrics with confidence intervals
print_performance_metrics <- function(performance_summary) {
    print("Overall Performance Metrics (Mean ± SD):")
    for(metric in names(performance_summary)) {
        if(startsWith(metric, "Mean_")) {
            metric_name <- gsub("Mean_", "", metric)
            mean_val <- performance_summary[[metric]]
            sd_val <- performance_summary[[paste0("SD_", metric_name)]]
            ci_lower <- mean_val - 1.96 * sd_val / sqrt(10)
            ci_upper <- mean_val + 1.96 * sd_val / sqrt(10)
            cat(sprintf("%s: %.3f ± %.3f (95%% CI: %.3f-%.3f)\n",
                        metric_name, mean_val, sd_val, ci_lower, ci_upper))
        }
    }
}

# ========================
# Visualization Functions
# ========================

# Create performance metrics boxplot
create_performance_boxplot <- function(results) {
    results_long <- results %>%
        select(-Fold) %>%
        gather(Metric, Value)

    ggplot(results_long, aes(x = Metric, y = Value)) +
        geom_boxplot(fill = "lightblue") +
        geom_jitter(width = 0.2, alpha = 0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Performance Metrics Across Folds",
             y = "Value",
             x = "Metric")
}

# Create ROC curve plot
plot_combined_roc <- function(roc_curves) {
    plot(0:1, 0:1, type = "n", xlim = c(1,0), ylim = c(0,1),
         xlab = "Specificity", ylab = "Sensitivity",
         main = "ROC Curves for All Folds")

    colors <- rainbow(length(roc_curves))
    for(i in seq_along(roc_curves)) {
        lines(roc_curves[[i]]$specificities, roc_curves[[i]]$sensitivities,
              col = colors[i], lwd = 2)
    }

    legend("bottomright", legend = paste("Fold", 1:length(roc_curves)),
           col = colors, lwd = 2)
}

# Create calibration plot
create_calibration_plot <- function(all_predictions) {
    ggplot(all_predictions, aes(x = Probability, fill = True_Label)) +
        geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
        theme_minimal() +
        labs(title = "Model Calibration Plot",
             x = "Predicted Probability",
             y = "Count")
}

# ========================
# Feature Importance Analysis
# ========================

analyze_feature_importance <- function(models_list, shap_results) {
    # Combine SHAP results
    all_shap <- bind_rows(shap_results)

    # Get importance metrics from each model
    importance_metrics <- lapply(1:length(models_list), function(i) {
        imp <- importance(models_list[[i]])
        data.frame(
            Feature = rownames(imp),
            MeanDecreaseAccuracy = imp[,3],
            MeanDecreaseGini = imp[,4],
            Fold = i
        )
    })

    all_importance <- bind_rows(importance_metrics)

    # Calculate average metrics
    feature_metrics <- all_importance %>%
        group_by(Feature) %>%
        summarize(
            Mean_MDA = mean(MeanDecreaseAccuracy),
            SD_MDA = sd(MeanDecreaseAccuracy),
            Mean_MDG = mean(MeanDecreaseGini),
            SD_MDG = sd(MeanDecreaseGini),
            Times_Selected = n()
        ) %>%
        arrange(desc(Mean_MDA))

    # Add SHAP values
    feature_metrics <- all_shap %>%
        group_by(Feature) %>%
        summarize(
            Mean_SHAP = mean(Importance),
            SD_SHAP = sd(Importance)
        ) %>%
        right_join(feature_metrics, by = "Feature")

    return(feature_metrics)
}

# ========================
# Generate All Plots
# ========================

generate_importance_plots <- function(feature_metrics) {
    # Top features plot
    p3 <- feature_metrics %>%
        head(20) %>%
        gather(Metric, Value, c(Mean_MDA, Mean_MDG, Mean_SHAP)) %>%
        ggplot(aes(x = reorder(Feature, Value), y = Value, fill = Metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() +
        scale_fill_viridis(discrete = TRUE) +
        theme_minimal() +
        labs(title = "Top 20 Features by Different Importance Metrics",
             x = "Feature",
             y = "Importance Score")

    # Metric correlation plot
    p4 <- ggplot(feature_metrics,
                 aes(x = Mean_MDA, y = Mean_SHAP)) +
        geom_point(aes(color = Mean_MDG), alpha = 0.6) +
        scale_color_viridis() +
        theme_minimal() +
        labs(title = "Correlation between Importance Metrics",
             x = "Mean Decrease Accuracy",
             y = "SHAP Importance",
             color = "Mean Decrease\nGini")

    # Feature stability plot
    p5 <- feature_metrics %>%
        head(20) %>%
        ggplot(aes(x = reorder(Feature, Times_Selected),
                   y = Times_Selected,
                   fill = Mean_MDA)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_viridis() +
        theme_minimal() +
        labs(title = "Feature Selection Stability (Top 20)",
             x = "Feature",
             y = "Times Selected",
             fill = "Mean Decrease\nAccuracy")

    return(list(p3 = p3, p4 = p4, p5 = p5))
}

# ========================
# Main Analysis Execution
# ========================

# Run performance analysis
print_performance_metrics(performance_summary)

# Calculate confusion matrix and metrics
all_predictions <- bind_rows(pred_results)
overall_conf_matrix <- table(
    Predicted = all_predictions$Predicted,
    Actual = all_predictions$True_Label
)

# Calculate per-class metrics
overall_metrics <- data.frame(
    Class = levels(all_predictions$True_Label),
    N = table(all_predictions$True_Label),
    Precision = diag(overall_conf_matrix) / colSums(overall_conf_matrix),
    Recall = diag(overall_conf_matrix) / rowSums(overall_conf_matrix)
)
overall_metrics$F1 <- 2 * (overall_metrics$Precision * overall_metrics$Recall) /
    (overall_metrics$Precision + overall_metrics$Recall)

# Generate all plots
p1 <- create_performance_boxplot(results)
p2 <- create_calibration_plot(all_predictions)
importance_plots <- generate_importance_plots(
    analyze_feature_importance(models_list, shap_results)
)

# ========================
# Save Results
# ========================

# Save plots
pdf(file.path(fig_dir, "20250214_model_analysis.pdf"), height = 6, width = 6)
print(p1)
plot_combined_roc(roc_curves)
print(p2)
print(importance_plots$p3)
print(importance_plots$p4)
print(importance_plots$p5)
dev.off()

# Save results
saveRDS(performance_summary,
        file.path(model_dir, "20250214_performance_summary.rds"))
write.csv(overall_metrics,
          file.path(model_dir, "20250214_class_metrics.csv"))

# Print final insights
cat("\nAdditional Insights:\n")
cat("1. Model Stability: CV of AUC =",
    sd(results$AUC)/mean(results$AUC) * 100, "%\n")
cat("2. Class Balance in Test Sets:\n")
print(table(all_predictions$True_Label))






















# FEATURE SELECTION_______________________________________________________________
packages <- c("tidyverse", "randomForest", "SummarizedExperiment", "sesame",
              "impute", "caret", "pROC", "readxl")
lapply(packages, require, character.only = TRUE)

folds <- readRDS("/home/lijz/20230818_thyroid_cancer/20250125_thyroid_folds.rds")
meta <- read_excel("/home/lijz/20230818_thyroid_cancer/20231102_thyroid_master.xlsx")
betas <- readRDS("/home/lijz/20230818_thyroid_cancer/20250125_thyroid88_betas_condensed_processed.rds")
meta <- meta %>% dplyr::filter(Include_In_Analysis == "1")
labels <- meta$Invasiveness
n_trees <- 500

model_dir <- "/home/lijz/20230818_thyroid_cancer/fmodel/"
feat_dir <- paste0(model_dir, "features/")
# dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)

n_probes <- nrow(betas)
batch_size <- 10000
n_batches <- ceiling(n_probes/batch_size)

for(i in 1:10) {
    print(paste("Starting fold", i))
    feat_log <- paste0(feat_dir, "20250125_", i, "features.csv")

    # Split data
    train_idx <- unlist(folds[-i])
    valid_idx <- folds[[i]]
    print(paste("Training samples:", length(train_idx)))

    # Initialize importance collection
    fold_importance <- data.frame()
    fold_models <- list()

    for(s in 0:(n_batches-1)) {
        print(paste("Processing batch", s, "/", n_batches - 1))

        # Calculate indices for current batch
        start_idx <- s * batch_size + 1
        end_idx <- min((s + 1) * batch_size, n_probes)

        sbetas <- betas[start_idx:end_idx, train_idx, drop=FALSE]
        print(paste("Batch dimensions:", dim(sbetas)[1], "x", dim(sbetas)[2]))

        model <- randomForest(
            x = t(sbetas),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )
        fold_models[[s+1]] <- model

        # Collect importance
        importance <- data.frame(
            Fold = i,
            Model = s,
            Feature = rownames(model$importance),
            MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"]
        )
        fold_importance <- rbind(fold_importance, importance)
    }
    saveRDS(fold_models, file.path(model_dir, paste0("20250125_models_fold", i, ".rds")))

    write.csv(fold_importance, feat_log)
    print("finished")
}

# FEATURE OPTIMIZATION ON HPC_______________________________________________________________

packages <- c("tidyverse", "randomForest", "SummarizedExperiment", "sesame",
              "impute", "caret", "pROC", "readxl", "shapviz")
lapply(packages, require, character.only = TRUE)

folds <- readRDS("/home/lijz/20230818_thyroid_cancer/20250125_thyroid_folds.rds")
meta <- read_excel("/home/lijz/20230818_thyroid_cancer/20231102_thyroid_master.xlsx")
betas <- readRDS("/home/lijz/20230818_thyroid_cancer/20250125_thyroid88_betas_condensed_processed.rds")
meta <- meta %>% dplyr::filter(Include_In_Analysis == "1")
labels <- meta$Invasiveness
n_trees <- 500

model_dir <- "/home/lijz/20230818_thyroid_cancer/fmodel/"
feat_dir <- paste0(model_dir, "features/")

# Function for model evaluation
evaluate_features <- function(feature_counts = c(1000, 2000, 3000, 5000), folds, betas, labels) {
    cv_results <- data.frame()

    for(n_feat in feature_counts) {
        print(paste("Testing", n_feat, "features"))
        for(i in 1:10) {
            print(paste("Fold", i))
            tryCatch({
                feat_file <- paste0(feat_dir, "20250125_", i, "features.csv")
                features <- read.csv(feat_file) %>%
                    arrange(desc(MeanDecreaseAccuracy)) %>%
                    head(n_feat) %>%
                    pull(Feature)

                train_idx <- unlist(folds[-i])
                valid_idx <- folds[[i]]

                model <- randomForest(
                    x = t(betas[features, train_idx]),
                    y = as.factor(labels[train_idx]),
                    ntree = 500
                )

                saveRDS(model, paste0(model_dir, "model_", n_feat, "_fold", i, ".rds"))

                pred_probs <- predict(model, t(betas[features, valid_idx]), type="prob")
                preds <- predict(model, t(betas[features, valid_idx]))

                roc_obj <- roc(labels[valid_idx], pred_probs[,2])
                cv_results <- rbind(cv_results, data.frame(
                    Features = n_feat,
                    Fold = i,
                    AUC = auc(roc_obj),
                    Accuracy = mean(preds == labels[valid_idx])
                ))
            }, error = function(e) {
                print(paste("Error in fold", i, ":", e$message))
            })
        }
        saveRDS(cv_results, paste0(model_dir, "20250125_feature_optimization_results.rds"))
    }
    return(cv_results)
}

# Run evaluation
results <- evaluate_features(feature_counts = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000),
                             folds = folds,
                             betas = betas,
                             labels = labels)

# MODEL DEVELOPMENT ON HPC_______________________________________________________________

# Load packages and check success
packages <- c("tidyverse", "randomForest", "SummarizedExperiment", "sesame",
              "impute", "caret", "pROC", "readxl", "shapviz")
success <- lapply(packages, require, character.only = TRUE)
print("Packages loaded")

# Load data
folds <- readRDS("/home/lijz/20230818_thyroid_cancer/20250125_thyroid_folds.rds")
meta <- read_excel("/home/lijz/20230818_thyroid_cancer/20231102_thyroid_master.xlsx")
betas <- readRDS("/home/lijz/20230818_thyroid_cancer/20250125_thyroid88_betas_condensed_processed.rds")
print(paste("Data dimensions - meta:", dim(meta)[1], "betas:", dim(betas)[1]))
n_trees = 500

meta <- meta %>% dplyr::filter(Include_In_Analysis == "1")
labels <- meta$Invasiveness
print(paste("Filtered samples:", length(labels)))

# Setup directories
model_dir <- "/home/lijz/20230818_thyroid_cancer/fmodel/"
feat_dir <- paste0(model_dir, "features/")
dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)

# Load features
features_all <- list()
print("Loading feature importance files...")
for(i in 1:10) {
    feat_file <- paste0(feat_dir, "20250125_", i, "features.csv")
    features_all[[i]] <- read.csv(feat_file) %>%
        arrange(desc(MeanDecreaseAccuracy))
    print(paste("Loaded fold", i, "features"))
}

# Select top features
top_features <- lapply(features_all, function(x) head(x$Feature, 3000))
saveRDS(top_features, paste0(model_dir, "20250125_top_features.rds"))

# Model evaluation
results <- data.frame()
pred_results <- list()  # Store all predictions

for(i in 1:10) {
    print(paste("Processing fold", i))
    train_idx <- unlist(folds[-i])
    valid_idx <- folds[[i]]

    print(paste("Training samples:", length(train_idx), "Validation samples:", length(valid_idx)))

    model <- randomForest(
        x = t(betas[top_features[[i]], train_idx]),
        y = as.factor(labels[train_idx]),
        ntree = n_trees
    )

    pred_probs <- predict(model, t(betas[top_features[[i]], valid_idx]), type="prob")
    preds <- predict(model, t(betas[top_features[[i]], valid_idx]))

    # Store predictions
    pred_results[[i]] <- data.frame(
        Fold = i,
        True_Label = labels[valid_idx],
        Predicted = preds,
        Probability = pred_probs[,2]
    )

    roc_obj <- roc(labels[valid_idx], pred_probs[,2])
    results <- rbind(results, data.frame(
        Fold = i,
        AUC = auc(roc_obj),
        Accuracy = mean(preds == labels[valid_idx])
    ))

    # Save model and predictions
    saveRDS(model, paste0(model_dir, "20250125_final_model_fold", i, ".rds"))
    saveRDS(pred_results[[i]], paste0(model_dir, "20250125_predictions_fold", i, ".rds"))
    print(paste("Completed fold", i))

    # SHAP analysis
    # shap <- shapviz(model, X_pred = t(betas[top_features[[i]], valid_idx]))
    # saveRDS(shap, paste0(model_dir, "20250125_shap_fold", i, ".rds"))
}

# Save final results
saveRDS(results, paste0(model_dir, "20250125_cv_results.rds"))
saveRDS(pred_results, paste0(model_dir, "20250125_all_predictions.rds"))
write.csv(results, paste0(model_dir, "20250125_cv_results.csv"))

print("Final Results:")
print(summary(results))

# MODEL ANALYSIS ON COMPUTER_______________________________________________________________

sum <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/analysis_summary.rds")
res <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/cv_results.rds")
folds <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/cv_folds.rds")


pred <- readRDS(here("data", "20250125_all_predictions.rds"))
# write.csv(results, paste0(model_dir, "20250125_cv_results.csv"))

# SHAP analysis
shap <- shapviz(model, X_pred = t(betas[top_features[[i]], valid_idx]))
saveRDS(shap, paste0(model_dir, "20250125_shap_fold", i, ".rds"))





