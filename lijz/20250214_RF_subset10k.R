# ========================
# Setup and Dependencies
# ========================

required_packages <- c(
    "tidyverse", "randomForest", "SummarizedExperiment", "sesame", "sail",
    "impute", "RSNNS", "e1071", "caret", "ISLR", "pROC", "readxl", "shapviz",
    "future", "future.apply", "logger", "shapr", "here", "viridis"
)
load_packages <- function(packages) {
    for (package in packages) {
        if (!require(package, character.only = TRUE)) {
            stop(paste("Package", package, "is required but not installed"))
        }
    }
}
load_packages(required_packages)

# MAKE FOLDS ON MY COMPUTER_______________________________________________________________
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")
ss = ss %>% dplyr::filter(Include_In_Analysis == "1")

# Create folds
set.seed(42)
folds <- createFolds(ss$Invasiveness, k=10)
saveRDS(folds, "~/Documents/HPC_share/thyroid/data/20250125_thyroid_folds.rds")

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
# Set up directories
model_dir <- file.path(here(), "models")
feat_dir <- file.path(model_dir, "features")
dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
folds <- readRDS(file.path(here(), "data", "20250125_thyroid_folds.rds"))
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")

# Load and preprocess betas
betas <- readRDS(here("data", "20241029_thyroid88_betas_processed.rds"))
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

subset_feature_selection <- function(betas, batch_size, n_feat, train_idx, labels, n_trees = 500, fold = NULL) {
    n_probes <- nrow(betas)
    n_batches <- ceiling(n_probes/batch_size)
    fold_models <- list()
    fold_importance <- data.frame()

    for(s in 0:(n_batches-1)) {
        # Calculate indices for current batch
        start_idx <- s * batch_size + 1
        end_idx <- min((s + 1) * batch_size, n_probes)

        sbetas <- betas[start_idx:end_idx, drop=FALSE]
        print(paste("Batch dimensions:", dim(sbetas)[1], "x", dim(sbetas)[2]))

        model <- randomForest(
            x = t(sbetas[, train_idx]),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )
        fold_models[[s+1]] <- model

        importance <- data.frame(
            Fold = if(is.null(fold)) NA else fold,
            Model = s,
            Feature = rownames(model$importance),
            MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"]
        )
        fold_importance <- rbind(fold_importance, importance)
    }
    features <- fold_importance %>%
        arrange(desc(MeanDecreaseAccuracy)) %>%
        head(n_feat) %>%
        pull(Feature)

    return(features)
}

for(i in 1:10) {
    print(paste("Processing fold", i))
    train_idx <- unlist(folds[-i])
    valid_idx <- folds[[i]]

    # Feature selection subset 10k
    sel_probes <- subset_feature_selection(
        betas = betas,
        batch_size = 10000,
        n_feat = 3000,
        train_idx = train_idx,
        labels = labels,
        n_trees = n_trees,
        fold = i
    )

    train_betas <- betas[sel_probes, train_idx]
    feature_lists[[i]] <- sel_probes

    # Train model
    model <- randomForest(
        x = t(train_betas),
        y = as.factor(labels[train_idx]),
        ntree = n_trees,
        importance = TRUE
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
    pred_probs <- predict(model, t(betas[sel_probes, valid_idx]), type="prob")
    preds <- predict(model, t(betas[sel_probes, valid_idx]))

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
# Summary Statistics
# ========================
# Total unique features: 3857
dim(avg_importance)

# Features selected in all folds: 2235
feature_overlap_10 = feature_overlap %>% filter(Times_Selected == 10)

# Top 10 most important features by mean importance:
print(head(avg_importance, 10))

# ========================
# Model Analysis
# ========================
# Set up directories
model_dir <- file.path(here(), "data", "20250214")
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

print_performance_metrics(performance_summary)
# "Overall Performance Metrics (Mean ± SD):"
# AUC: 0.901 ± 0.175 (95% CI: 0.792-1.009)
# Accuracy: 0.864 ± 0.118 (95% CI: 0.791-0.936)
# Sensitivity: 0.930 ± 0.120 (95% CI: 0.856-1.004)
# Specificity: 0.783 ± 0.264 (95% CI: 0.620-0.947)
# PPV: 0.883 ± 0.139 (95% CI: 0.797-0.969)
# NPV: 0.902 ± 0.162 (95% CI: 0.801-1.002)

# [1] "Overall Performance Metrics (Mean ± SD):"
# AUC: 0.908 ± 0.131 (95% CI: 0.826-0.989)
# Accuracy: 0.867 ± 0.123 (95% CI: 0.791-0.944)
# Sensitivity: 0.936 ± 0.114 (95% CI: 0.865-1.007)
# Specificity: 0.775 ± 0.272 (95% CI: 0.606-0.944)
# PPV: 0.873 ± 0.149 (95% CI: 0.781-0.965)
# NPV: 0.877 ± 0.202 (95% CI: 0.751-1.002)

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
pdf(file.path(fig_dir, "20250215_thyroid88_model_summary_boxplots.pdf"), height = 6, width = 6)
print(p1)
dev.off()

# think there's prob with this one!!
pdf(file.path(fig_dir, "20250215_thyroid88_model_calibration.pdf"), height = 6, width = 6)
print(p2)
dev.off()

pdf(file.path(fig_dir, "20250215_thyroid88_model_roc_curves.pdf"), height = 6, width = 6)
plot_combined_roc(roc_curves)
dev.off()

pdf(file.path(fig_dir, "20250215_thyroid88_model_importance1.pdf"), height = 6, width = 6)
print(importance_plots$p3)
dev.off()

pdf(file.path(fig_dir, "20250215_thyroid88_model_importance2.pdf"), height = 6, width = 6)
print(importance_plots$p4)
dev.off()

pdf(file.path(fig_dir, "20250215_thyroid88_model_importance3.pdf"), height = 6, width = 6)
print(importance_plots$p5)
dev.off()

# Save results
saveRDS(performance_summary,
        file.path(model_dir, "20250215_performance_summary.rds"))
write.csv(overall_metrics,
          file.path(model_dir, "20250215_class_metrics.csv"))

# Print final insights
cat("\nAdditional Insights:\n")
cat("1. Model Stability: CV of AUC =",
    sd(results$AUC)/mean(results$AUC) * 100, "%\n")
cat("2. Class Balance in Test Sets:\n")
print(table(all_predictions$True_Label))

# 1. Model Stability: CV of AUC = 14.48302 %
# High  Low
# 55   33



