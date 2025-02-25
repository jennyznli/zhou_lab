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

                             
