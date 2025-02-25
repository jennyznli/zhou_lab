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
dir <- file.path(here(), "data", "20250222_mostvar")
model_dir <- file.path(dir, "models")
pred_dir <- file.path(dir, "predictions")
feat_dir <- file.path(dir, "features")
sum_dir <- file.path(dir, "summary")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)
date <- "20250222"

# Load data
folds <- readRDS(file.path(here(), "data", "20250125_thyroid_folds.rds"))
ss <- read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    filter(Include_In_Analysis == "1")

# Load and preprocess betas
betas <- readRDS(here("data", "20250214_thyroid88_processed_condensed_betas.rds"))
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
imp_results <- list()
feature_lists <- list()

# ========================
# Cross-Validation Loop
# ========================

# 3K MOST VARIABLE
for(i in 1:10) {
    print(paste("Processing fold", i))
    train_idx <- unlist(folds[-i])
    valid_idx <- folds[[i]]

    # FEATURE SELECTION
    bSubMostVariable <- function(betas, n=3000) {
        std <- apply(betas, 1, sd, na.rm=TRUE)
        betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
    }
    train_betas <- bSubMostVariable(betas[,train_idx])
    var_probes <- rownames(train_betas)
    feature_lists[[i]] <- var_probes  # Store selected features

    # TRAIN MODEL
    model <- randomForest(
        x = t(train_betas),
        y = as.factor(labels[train_idx]),
        ntree = n_trees,
        importance = TRUE
    )
    models_list[[i]] <- model

    # IMPORTANCE SCORES
    importance_scores <- as.data.frame(importance(model))
    imp_results[[i]] <- data.frame(
        Feature = rownames(importance_scores),
        Importance = importance_scores$MeanDecreaseAccuracy,
        Gini = importance_scores$MeanDecreaseGini,
        Fold = i
    )

    # PREDICTIONS
    pred_probs <- predict(model, t(betas[var_probes, valid_idx]), type="prob")
    preds <- predict(model, t(betas[var_probes, valid_idx]))

    pred_results[[i]] <- data.frame(
        Fold = i,
        True_Label = labels[valid_idx],
        Predicted = preds,
        Probability = pred_probs[,2],
        Sample_ID = valid_idx
    )

    # ROC and metrixs
    roc_obj <- roc(labels[valid_idx], pred_probs[,2])
    roc_curves[[i]] <- roc_obj
    conf_matrix <- caret::confusionMatrix(as.factor(preds), as.factor(labels[valid_idx]))

    # SUMMARIZED RESULTS
    results <- rbind(results, data.frame(
        Fold = i,
        AUC = auc(roc_obj),
        Accuracy = mean(preds == labels[valid_idx]),
        Sensitivity = conf_matrix$byClass["Sensitivity"],
        Specificity = conf_matrix$byClass["Specificity"],
        PPV = conf_matrix$byClass["Pos Pred Value"],
        NPV = conf_matrix$byClass["Neg Pred Value"],
        F1 = conf_matrix$byClass["F1"]
    ))

    # SAVE
    saveRDS(model, file.path(model_dir, paste0(date, "_final_model_fold", i, ".rds")))
    saveRDS(pred_results[[i]], file.path(pred_dir, paste0(date, "_predictions_fold", i, ".rds")))
    write.csv(imp_results[[i]], file.path(feat_dir, paste0(date, "_imp_values_fold", i, ".csv")))

    print(paste("Completed fold", i))
}

# ========================
# Process Feature Importance
# ========================
all_imp <- bind_rows(imp_results)  #30000     4
avg_importance <- all_imp %>%
    group_by(Feature) %>%
    summarize(
        Mean_Importance = mean(Importance),
        SD_Importance = sd(Importance),
        Times_Selected = n()
    ) %>%
    arrange(desc(Mean_Importance))

# ========================
# Save Final Results
# ========================
saveRDS(models_list, file.path(sum_dir, paste0(date, "_final_model_list.rds")))
saveRDS(pred_results, file.path(sum_dir, paste0(date, "_pred_results.rds")))
saveRDS(results, file.path(sum_dir, paste0(date, "_results_df.rds")))
saveRDS(roc_curves, file.path(sum_dir, paste0(date, "_roc_curves.rds")))

# Save feature importance results
write.csv(avg_importance, file.path(sum_dir, paste0(date, "_average_importance.csv")))
saveRDS(imp_results, file.path(sum_dir, paste0(date, "_imp_results_list.rds")))
saveRDS(feature_lists, file.path(sum_dir, paste0(date, "_feature_lists.rds")))

# ========================
# Summary Statistics
# ========================
# Total unique features: 3824
dim(avg_importance)

# Features selected in all folds: 2239
feature_overlap_10 = avg_importance %>% filter(Times_Selected == 10)
dim(feature_overlap_10)

# Top 10 most important features by mean importance:
print(head(avg_importance, 10))

# ========================
# Model Analysis
# ========================
dir <- file.path(here(), "data", "20250223_subset")
model_dir <- file.path(dir, "models")
pred_dir <- file.path(dir, "predictions")
feat_dir <- file.path(dir, "features")
sum_dir <- file.path(dir, "summary")
date <- "20250223"
fig_dir <- file.path(here(), "figures")

# Load all required data
# feature_lists_rf <- readRDS(file.path(feat_dir, "20250214_feature_lists.rds"))
# feature_lists_mv <- readRDS("~/Documents/HPC_share/thyroid/data/20250222_mostvar/summary/20250222feature_lists.rds")

results <- read.csv(file.path(sum_dir, paste0(date, "_average_metrics.csv")))
pred_results <- readRDS(file.path(model_dir, paste0(date,  "_pred_results.rds")))
roc_curves <- readRDS(file.path(sum_dir, paste0(date, "_roc_curves.rds")))
models_list <- readRDS(file.path(sum_dir, paste0(date,  "_final_model_list.rds")))
imp_results <- readRDS(file.path(sum_dir, paste0(date,  "_imp_results_list.rds")))

#
# ========================
# Feature Stability
# ========================
# Create frequency table
# 13678 unique ones
# feature_overlap <- data.frame(
#     Feature = unique(unlist(feature_lists)),
#     Times_Selected = sapply(
#         unique(unlist(feature_lists)),
#         function(x) sum(sapply(feature_lists, function(y) x %in% y))
#     )
# ) #9159    2
# feature_overlap_mv <- data.frame(
#     Feature = unique(unlist(feature_lists_mv)),
#     Times_Selected = sapply(
#         unique(unlist(feature_lists_mv)),
#         function(x) sum(sapply(feature_lists_mv, function(y) x %in% y))
#     )
# ) #3824    2

# cumulative_data <- bind_rows(
#     # Method 1 (most variable)
#     data.frame(
#         Method = "Most Variable",
#         Folds = 1:10,
#         Proportion = sapply(1:10, function(x) {
#             sum(feature_overlap_mv$Times_Selected >= x) / nrow(feature_overlap_mv)
#         })
#     ),
#     # Method 2 (RF importance)
#     data.frame(
#         Method = "RF Importance",
#         Folds = 1:10,
#         Proportion = sapply(1:10, function(x) {
#             sum(feature_overlap_rf$Times_Selected >= x) / nrow(feature_overlap_rf)
#         })
#     )
# )

# p <- ggplot(cumulative_data, aes(x = Folds, y = Proportion, color = Method)) +
#     geom_line(linewidth = 1.2) +
#     geom_point(size = 3) +
#     scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
#     scale_x_continuous(breaks = 1:10) +
#     scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#     labs(
#         title = "Cumulative Feature Stability",
#         x = "Number of Folds a Feature Appears In (≥X)",
#         y = "Proportion of Features"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#         panel.grid.major = element_line(color = "grey90"),
#         panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5, face = "bold")
#     )

feature_lists <- list.files(feat_dir, pattern = "_importance_fold", full.names = TRUE) |>
    lapply(readRDS)

count_feature_overlap <- function(feature_lists) {
    all_features <- lapply(feature_lists, function(mat) mat[,1])
    combined_features <- unlist(all_features)
    feature_counts <- table(combined_features)
    overlap_df <- data.frame(
        feature = names(feature_counts),
        count = as.numeric(feature_counts),
        percentage = as.numeric(feature_counts) / length(feature_lists) * 100
    )
    overlap_df <- overlap_df[order(-overlap_df$count), ]
    return(overlap_df)
}

overlap_results <- count_feature_overlap(feature_lists)

cumulative_data <- data.frame(
    Folds = 1:10,
    Proportion = sapply(1:10, function(x) {
        sum(overlap_results$count >= x) / nrow(overlap_results)
    })
)

p <- ggplot(cumulative_data, aes(x = Folds, y = Proportion)) +
    geom_line(linewidth = 1.2, color = "#2c7bb6") +
    geom_point(size = 3, color = "#2c7bb6") +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "Cumulative Feature Stability",
        x = "Number of Folds a Feature Appears In (≥X)",
        y = "Proportion of Features"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

pdf(file.path(fig_dir, paste0(date, "_thyroid88_feature_stability_rf2.pdf")), height = 6, width = 6)
print(p)
dev.off()

# p <- ggplot(overlap_results, aes(x = count)) +
#     geom_histogram(binwidth = 1) +
#     theme_minimal() +
#     labs(x = "Number of Lists", y = "Count of Features",
#          title = "Distribution of Feature Overlap")
#
# pdf(file.path(fig_dir, paste0(date, "_thyroid88_feature_stability_rf.pdf")), height = 6, width = 6)
# print(p)
# dev.off()

data <- data.frame(
    Folds = 1:10,
    Proportion = sapply(1:10, function(x) {
        sum(feature_counts_df$Count >= x) / nrow(feature_counts_df)
    })
)

p <- ggplot(cumulative_data, aes(x = Folds, y = Proportion)) +
    geom_line(linewidth = 1.2, color = "#2c7bb6") +
    geom_point(size = 3, color = "#2c7bb6") +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "Cumulative Feature Stability",
        x = "Number of Folds a Feature Appears In (≥X)",
        y = "Proportion of Features"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

pdf(file.path(fig_dir, paste0(date, "_thyroid88_feature_stability_rf.pdf")), height = 6, width = 6)
print(p)
dev.off()





feature_lists <- list.files(feat_dir, pattern = "_importance_fold", full.names = TRUE) |>
    lapply(readRDS)




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
        Mean_F1 = mean(F1),
        SD_F1 = sd(F1),
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

# i don't remember what this is...
# "Overall Performance Metrics (Mean ± SD):"
# AUC: 0.901 ± 0.175 (95% CI: 0.792-1.009)
# Accuracy: 0.864 ± 0.118 (95% CI: 0.791-0.936)
# Sensitivity: 0.930 ± 0.120 (95% CI: 0.856-1.004)
# Specificity: 0.783 ± 0.264 (95% CI: 0.620-0.947)
# PPV: 0.883 ± 0.139 (95% CI: 0.797-0.969)
# NPV: 0.902 ± 0.162 (95% CI: 0.801-1.002)

# subsets of 10k
# [1] "Overall Performance Metrics (Mean ± SD):"
# AUC: 0.908 ± 0.131 (95% CI: 0.826-0.989)
# Accuracy: 0.878 ± 0.140 (95% CI: 0.791-0.965)
# Sensitivity: 0.936 ± 0.114 (95% CI: 0.865-1.007)
# Specificity: 0.783 ± 0.315 (95% CI: 0.588-0.978)
# PPV: 0.893 ± 0.139 (95% CI: 0.807-0.979)
# NPV: 0.860 ± 0.227 (95% CI: 0.719-1.001)

# subsets of 10k (second try)
# AUC: 0.908 ± 0.131 (95% CI: 0.826-0.989)
# Accuracy: 0.878 ± 0.119 (95% CI: 0.804-0.952)
# Sensitivity: 0.775 ± 0.272 (95% CI: 0.606-0.944)
# Specificity: 0.952 ± 0.077 (95% CI: 0.905-1.000)
# F1: 0.808 ± 0.220 (95% CI: 0.672-0.944)

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

analyze_feature_importance <- function(models_list, imp_results) {
    # Combine SHAP results
    all_shap <- bind_rows(imp_results)

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
    analyze_feature_importance(models_list, imp_results)
)

# ========================
# Save Results
# ========================
pdf(file.path(fig_dir, paste0(date, "_thyroid88_model_summary_boxplots.pdf")), height = 6, width = 6)
print(p1)
dev.off()

# think there's prob with this one!!
pdf(file.path(fig_dir, paste0(date, "_thyroid88_model_calibration.pdf")), height = 6, width = 6)
print(p2)
dev.off()

pdf(file.path(fig_dir, paste0(date, "_thyroid88_model_roc_curves.pdf")), height = 6, width = 6)
plot_combined_roc(roc_curves)
dev.off()

pdf(file.path(fig_dir, paste0(date, "_thyroid88_model_importance1.pdf")), height = 6, width = 6)
print(importance_plots$p3)
dev.off()

pdf(file.path(fig_dir, paste0(date, "_thyroid88_model_importance2.pdf")), height = 6, width = 6)
print(importance_plots$p4)
dev.off()

pdf(file.path(fig_dir, paste0(date, "_thyroid88_model_importance3.pdf")), height = 6, width = 6)
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

# ========================
# Final Model Development
# ========================





# ========================
# SHAP Analysis
# ========================

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
                   v = "cg00768487",
                   color_var = "cg11868461")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

# a little less correlated...
p <- sv_dependence(shp,
                   v = "cg00768487",
                   color_var = "cg10579279")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence2.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

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


