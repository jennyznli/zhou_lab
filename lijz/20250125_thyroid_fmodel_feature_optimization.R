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
                             
                             
                             
