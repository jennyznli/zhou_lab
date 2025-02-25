# Script: Thyroid Cancer Methylation Analysis Pipeline
# Description: Analysis pipeline for thyroid cancer methylation data
# Last Updated: 2025-01-31

# Load required packages
required_packages <- c(
    "tidyverse",
    "randomForest",
    "SummarizedExperiment",
    "sesame",
    "impute",
    "caret",
    "pROC",
    "readxl",
    "shapviz"
)

# Function to load required packages
load_packages <- function(packages) {
    for (package in packages) {
        if (!require(package, character.only = TRUE)) {
            stop(paste("Package", package, "is required but not installed"))
        }
    }
}

# Setup parallel processing through randomForest
setup_parallel <- function() {
    threads <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
    options(rf.cores = threads)
    message("Using ", threads, " threads for randomForest")
}

# Data Loading Functions
load_data <- function(meta_path, betas_path) {
    message("Loading data...")

    tryCatch({
        meta <- readxl::read_excel(meta_path)
        betas <- readRDS(betas_path)

        # Filter included samples
        meta_filtered <- meta %>%
            dplyr::filter(Include_In_Analysis == "1")

        betas_filtered <- betas[, colnames(betas) %in% meta_filtered$IDAT]

        # Validate dimensions
        if (ncol(betas_filtered) != nrow(meta_filtered)) {
            stop("Dimension mismatch between betas and metadata")
        }

        message(sprintf("Loaded data: %d features, %d samples",
                      nrow(betas_filtered), ncol(betas_filtered)))

        return(list(
            betas = betas_filtered,
            meta = meta_filtered
        ))
    }, error = function(e) {
        stop(paste("Error loading data:", e$message))
    })
}

# Data Preprocessing Functions
clean_matrix <- function(mtx, f_row = 0.5, f_col = 0.5) {
    message("Cleaning matrix...")

    namtx <- is.na(mtx)
    good_row <- rowSums(namtx) <= ncol(mtx) * f_row
    good_col <- colSums(namtx) <= nrow(mtx) * f_col

    message(sprintf("Filtering: %d/%d rows and %d/%d columns retained",
                   sum(good_row), nrow(mtx),
                   sum(good_col), ncol(mtx)))

    mtx[good_row, good_col]
}

impute_row_mean <- function(mtx) {
    message("Imputing missing values...")

    k <- which(is.na(mtx), arr.ind = TRUE)
    if (length(k) > 0) {
        mtx[k] <- rowMeans(mtx, na.rm = TRUE)[k[,1]]
        message(sprintf("Imputed %d missing values", length(k)))
    }
    mtx
}

# Create Cross-validation Folds
create_cv_folds <- function(meta, k = 10, seed = 42) {
    message("Creating cross-validation folds...")

    set.seed(seed)
    folds <- createFolds(meta$Invasiveness, k = k)

    message(sprintf("Created %d folds", length(folds)))
    return(folds)
}

# Feature Selection Function
select_features_batch <- function(betas, labels, train_idx, fold_idx, output_dir,
                                n_trees = 500, batch_size = 10000) {
    n_probes <- nrow(betas)
    n_batches <- ceiling(n_probes/batch_size)
    importance_list <- vector("list", n_batches)

    feat_dir <- file.path(output_dir, "feature_importance")
    dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)

    for (b in 1:n_batches) {
        message(sprintf("Processing batch %d/%d for fold %d", b, n_batches, fold_idx))

        start_idx <- (b-1) * batch_size + 1
        end_idx <- min(b * batch_size, n_probes)
        sbetas <- betas[start_idx:end_idx, train_idx, drop = FALSE]

        model <- randomForest(
            x = t(sbetas),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )

        importance_list[[b]] <- data.frame(
            Fold = fold_idx,
            Batch = b,
            Feature = rownames(sbetas),
            MeanDecreaseAccuracy = model$importance[, "MeanDecreaseAccuracy"],
            MeanDecreaseGini = model$importance[, "MeanDecreaseGini"]
        )
    }

    all_importance <- do.call(rbind, importance_list) %>%
        arrange(desc(MeanDecreaseAccuracy))

    # Save all features
    write.csv(all_importance,
              file.path(feat_dir, sprintf("fold%d_all_features.csv", fold_idx)),
              row.names = FALSE)

    # Save top features
    top_features <- all_importance %>%
        head(3000) %>%
        select(Feature, MeanDecreaseAccuracy, MeanDecreaseGini)

    write.csv(top_features,
              file.path(feat_dir, sprintf("fold%d_top_features.csv", fold_idx)),
              row.names = FALSE)

    return(all_importance)
}

# Main Analysis Function
run_analysis <- function(meta_path, betas_path, output_dir,
                        n_trees = 500, n_folds = 10) {
    # Setup
    message("Starting analysis...")
    load_packages(required_packages)
    setup_parallel()

    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Load and preprocess data
    data <- load_data(meta_path, betas_path)
    processed_betas <- clean_matrix(data$betas) %>% impute_row_mean()

    # Create folds
    folds <- create_cv_folds(data$meta, k = n_folds)
    saveRDS(folds, file.path(output_dir, "cv_folds.rds"))

    # Feature selection and model training for each fold
    results <- vector("list", n_folds)

    for (i in 1:n_folds) {
        message(sprintf("Processing fold %d/%d", i, n_folds))

        # Feature selection
        importance <- select_features_batch(
            processed_betas,
            data$meta$Invasiveness,
            unlist(folds[-i]),
            i,
            output_dir,
            n_trees = n_trees
        )

        # Select top features
        top_features <- importance %>%
            arrange(desc(MeanDecreaseAccuracy)) %>%
            head(3000) %>%
            pull(Feature)

        # Train and evaluate model
        train_idx <- unlist(folds[-i])
        valid_idx <- folds[[i]]

        model <- randomForest(
            x = t(processed_betas[top_features, train_idx]),
            y = as.factor(data$meta$Invasiveness[train_idx]),
            ntree = n_trees
        )

        # Make predictions
        pred_probs <- predict(model, t(processed_betas[top_features, valid_idx]), type = "prob")
        preds <- predict(model, t(processed_betas[top_features, valid_idx]))

        # Calculate metrics
        roc_obj <- roc(data$meta$Invasiveness[valid_idx], pred_probs[,2])

        # Store results
        results[[i]] <- list(
            fold = i,
            metrics = data.frame(
                AUC = auc(roc_obj),
                Accuracy = mean(preds == data$meta$Invasiveness[valid_idx])
            ),
            predictions = data.frame(
                True_Label = data$meta$Invasiveness[valid_idx],
                Predicted = preds,
                Probability = pred_probs[,2]
            )
        )

        # Save fold-specific results
        saveRDS(model, file.path(output_dir, sprintf("final_model_fold%d.rds", i)))
    }

    # Aggregate and save final results
    final_metrics <- do.call(rbind, lapply(results, function(x) x$metrics))
    saveRDS(final_metrics, file.path(output_dir, "cv_results.rds"))

    # Create summary
    summary <- list(
        mean_auc = mean(final_metrics$AUC),
        sd_auc = sd(final_metrics$AUC),
        mean_accuracy = mean(final_metrics$Accuracy),
        sd_accuracy = sd(final_metrics$Accuracy)
    )

    saveRDS(summary, file.path(output_dir, "analysis_summary.rds"))

    message("Analysis complete")
    return(summary)
}


# MODEL ANALYSIS ON COMPUTER_______________________________________________________________

sum <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/analysis_summary.rds")
res <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/cv_results.rds")
folds <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/cv_folds.rds")
mod1 <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/final_model_fold1.rds")

feat1_all <- read.csv("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/feature_importance/fold1_all_features.csv")
feat1_3k <- read.csv("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/feature_importance/fold1_top_features.csv")
feat_prev1 <- read_tsv("/Users/jennyzli/Documents/HPC_share/thyroid/prev_data/20240610_RF_850k_features.tsv")
feat_prev2 <- read_tsv("/Users/jennyzli/Documents/HPC_share/thyroid/features/20240610_features.tsv")

colnames(feat_prev1) <- c("CpG", "MeanDecreaseAccuracy", "MeanDecreaseGini")
feat_prev1 <- dplyr::arrange(feat_prev1, desc(MeanDecreaseAccuracy))
feat1_all <- dplyr::arrange(feat1_all, desc(MeanDecreaseAccuracy))

# SHAP analysis
shap <- shapviz(model, X_pred = t(betas[top_features[[i]], valid_idx]))
saveRDS(shap, paste0(model_dir, "20250125_shap_fold", i, ".rds"))


# 1. Check feature selection stability
# analyze_feature_overlap <- function(base_dir) {
#     features_by_fold <- list()
#     for(i in 1:10) {
#         features_by_fold[[i]] <- read.csv(
#             file.path(base_dir, "feature_importance",
#                       sprintf("fold%d_top_features.csv", i))
#         )$Feature
#         print(head(features_by_fold[1]))
#     }
#     # Check overlap between folds
#     common_features <- Reduce(intersect, features_by_fold)
#     return(length(common_features))
# }
# feature_overlap <- analyze_feature_overlap("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/")
#
#
# base_dir <- "/Users/jennyzli/Documents/HPC_share/thyroid/data/final/"
# features_by_fold <- list()
# for(i in 1:10) {
#     features_by_fold[[i]] <- read.csv(
#         file.path(base_dir, "feature_importance",
#                   sprintf("fold%d_top_features.csv", i))
#     )$Feature
#     print(head(features_by_fold[1]))
#
# }
# common_features <- Reduce(intersect, features_by_fold)


## testing this function
ss = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")
ss = ss %>% dplyr::filter(Include_In_Analysis == "1")
folds <- readRDS("/Users/jennyzli/Documents/HPC_share/thyroid/data/final/cv_folds.rds")
betas <- readRDS("~/Documents/HPC_share/thyroid/data/20250125_thyroid88_betas_condensed_processed.rds")
labels <- ss$Invasiveness
unlist(folds[-1])


# Feature selection

output_dir <- "~/Documents/HPC_share/thyroid/data/"

select_features_batch <- function(betas, labels, train_idx, fold_idx, output_dir,
                                  n_trees = 500, batch_size = 10000) {
    n_probes <- nrow(betas)
    n_batches <- ceiling(n_probes/batch_size)
    importance_list <- vector("list", n_batches)

    feat_dir <- file.path(output_dir, "feature_importance")
    dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)

    for (b in 1:n_batches) {
        message(sprintf("Processing batch %d/%d for fold %d", b, n_batches, fold_idx))

        start_idx <- (b-1) * batch_size + 1
        end_idx <- min(b * batch_size, n_probes)
        sbetas <- betas[start_idx:end_idx, train_idx, drop = FALSE]
        print(start_idx)
        print(end_idx)
#         model <- randomForest(
#             x = t(sbetas),
#             y = as.factor(labels[train_idx]),
#             ntree = n_trees,
#             importance = TRUE
#         )
#
#         importance_list[[b]] <- data.frame(
#             Fold = fold_idx,
#             Batch = b,
#             Feature = rownames(sbetas),
#             MeanDecreaseAccuracy = model$importance[, "MeanDecreaseAccuracy"],
#             MeanDecreaseGini = model$importance[, "MeanDecreaseGini"]
#         )
    }

    # all_importance <- do.call(rbind, importance_list) %>%
    #     arrange(desc(MeanDecreaseAccuracy))
    #
    # # Save all features
    # write.csv(all_importance,
    #           file.path(feat_dir, sprintf("fold%d_all_features.csv", fold_idx)),
    #           row.names = FALSE)
    #
    # # Save top features
    # top_features <- all_importance %>%
    #     head(3000) %>%
    #     select(Feature, MeanDecreaseAccuracy, MeanDecreaseGini)
    #
    # write.csv(top_features,
    #           file.path(feat_dir, sprintf("fold%d_top_features.csv", fold_idx)),
    #           row.names = FALSE)
    #
    # return(all_importance)
}
importance <- select_features_batch(
    betas,
    labels,
    unlist(folds[-1]),
    1,
    output_dir,
    n_trees = 500
)
train_idx <- unlist(folds[-1])
start_idx <- 790001
end_idx <- 800000
sbetas <- betas[start_idx:end_idx, train_idx, drop = FALSE]



run_analysis <- function(meta_path, betas_path, output_dir,
                         n_trees = 500, n_folds = 10) {
    # Setup
    message("Starting analysis...")
    load_packages(required_packages)
    setup_parallel()

    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Load and preprocess data
    data <- load_data(meta_path, betas_path)
    processed_betas <- clean_matrix(data$betas) %>% impute_row_mean()

    # Create folds
    folds <- create_cv_folds(data$meta, k = n_folds)
    saveRDS(folds, file.path(output_dir, "cv_folds.rds"))

    # Feature selection and model training for each fold
    results <- vector("list", n_folds)

    for (i in 1:n_folds) {
        message(sprintf("Processing fold %d/%d", i, n_folds))

        # Feature selection
        importance <- select_features_batch(
            processed_betas,
            data$meta$Invasiveness,
            unlist(folds[-i]),
            i,
            output_dir,
            n_trees = n_trees
        )

        # Select top features
        top_features <- importance %>%
            arrange(desc(MeanDecreaseAccuracy)) %>%
            head(3000) %>%
            pull(Feature)

        # Train and evaluate model
        train_idx <- unlist(folds[-i])
        valid_idx <- folds[[i]]

        model <- randomForest(
            x = t(processed_betas[top_features, train_idx]),
            y = as.factor(data$meta$Invasiveness[train_idx]),
            ntree = n_trees
        )

        # Make predictions
        pred_probs <- predict(model, t(processed_betas[top_features, valid_idx]), type = "prob")
        preds <- predict(model, t(processed_betas[top_features, valid_idx]))

        # Calculate metrics
        roc_obj <- roc(data$meta$Invasiveness[valid_idx], pred_probs[,2])

        # Store results
        results[[i]] <- list(
            fold = i,
            metrics = data.frame(
                AUC = auc(roc_obj),
                Accuracy = mean(preds == data$meta$Invasiveness[valid_idx])
            ),
            predictions = data.frame(
                True_Label = data$meta$Invasiveness[valid_idx],
                Predicted = preds,
                Probability = pred_probs[,2]
            )
        )

        # Save fold-specific results
        saveRDS(model, file.path(output_dir, sprintf("final_model_fold%d.rds", i)))
    }

    # Aggregate and save final results
    final_metrics <- do.call(rbind, lapply(results, function(x) x$metrics))
    saveRDS(final_metrics, file.path(output_dir, "cv_results.rds"))

    # Create summary
    summary <- list(
        mean_auc = mean(final_metrics$AUC),
        sd_auc = sd(final_metrics$AUC),
        mean_accuracy = mean(final_metrics$Accuracy),
        sd_accuracy = sd(final_metrics$Accuracy)
    )

    saveRDS(summary, file.path(output_dir, "analysis_summary.rds"))

    message("Analysis complete")
    return(summary)
}






