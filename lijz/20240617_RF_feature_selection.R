
# PREPROCESSING ON HPC_______________________________________________________________
packages <- c("tidyverse", "randomForest", "SummarizedExperiment", "sesame", "sail",
              "impute", "RSNNS", "e1071", "caret", "ISLR", "pROC")
lapply(packages, require, character.only = TRUE)

meta = read_excel("/home/lijz/20230818_thyroid_cancer/20231102_thyroid_master.xlsx")
betas = readRDS("/home/lijz/20230818_thyroid_cancer/20240320_thyroid136_betas_condensed.rds")
meta = meta %>% dplyr::filter(INCLUDE_IN_ANALYSIS == 1)
print(dim(meta))

betas = betas[, colnames(betas) %in% meta$IDAT]
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

res <- readRDS(here("data", "20250125_cv_results.rds"))
pred <- readRDS(here("data", "20250125_all_predictions.rds"))
# write.csv(results, paste0(model_dir, "20250125_cv_results.csv"))

# SHAP analysis
shap <- shapviz(model, X_pred = t(betas[top_features[[i]], valid_idx]))
saveRDS(shap, paste0(model_dir, "20250125_shap_fold", i, ".rds"))


