.libPaths(c("/home/lijz/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

custom_lib <- "/home/lijz/R/x86_64-pc-linux-gnu-library/4.4"
if(dir.exists(custom_lib)) {
    .libPaths(c(custom_lib, .libPaths()))
} else {
    stop("Custom library path not found: ", custom_lib)
}

# ========================
# Load Packages
# ========================
required_packages <- c(
    "randomForest", "pROC", "caret", "mltools", 
    "readxl", "dplyr", "BiocParallel", "batchtools"
)

# Check package availability
missing <- setdiff(required_packages, rownames(installed.packages()))
if(length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse=", "))
}

suppressPackageStartupMessages({
    library(randomForest)
    library(pROC)
    library(caret)
    library(mltools)
    library(readxl)
    library(dplyr)
    library(BiocParallel)
    library(batchtools)
})


library(BiocParallel)
library(batchtools)
library(here)

# ========================
# Configuration
# ========================
dir <- file.path(here(), "20250222")
model_dir <- file.path(dir, "models")
pred_dir <- file.path(dir, "predictions")
feat_dir <- file.path(dir, "features")
sum_dir <- file.path(dir, "summary")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

date <- "20250222"

# Slurm configuration
slurm_config <- BatchtoolsParam(
    workers = 10,
    cluster = "slurm",
    template = "slurm_template.tmpl",
    resources = list(
        walltime = 86400,  # 24 hours
        memory = "64G",
        ncpus = 4
    )
)

# ========================
# Core Functions
# ========================
subset_feature_selection <- function(betas, batch_size, n_feat, train_idx, labels, n_trees = 500, fold = NULL) {
    n_probes <- nrow(betas)
    n_batches <- ceiling(n_probes/batch_size)
    fold_importance <- data.frame()

    for(s in 0:(n_batches-1)) {
        start_idx <- s * batch_size + 1
        end_idx <- min((s + 1) * batch_size, n_probes)
        sbetas <- betas[start_idx:end_idx,]
        
        model <- randomForest(
            x = t(sbetas[, train_idx]),
            y = as.factor(labels[train_idx]),
            ntree = n_trees,
            importance = TRUE
        )
        
        importance <- data.frame(
            Feature = rownames(model$importance),
            MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"],
            Fold = fold,
            Batch = s
        )
        fold_importance <- rbind(fold_importance, importance)
    }
    
    fold_importance %>%
        arrange(desc(MeanDecreaseAccuracy)) %>%
        head(n_feat) %>%
        pull(Feature)
}

process_fold <- function(i) {
    suppressPackageStartupMessages({
        library(randomForest)
        library(pROC)
        library(caret)
        library(mltools)
        library(readxl)
        library(dplyr)
    })
    
    # Load data fresh for each job
    dir <- file.path(here())
    ss <- readxl::read_excel(file.path(dir, "20231102_thyroid_master.xlsx")) |>
        dplyr::filter(Include_In_Analysis == "1")
    
    betas <- readRDS(file.path(dir, "20250214_thyroid88_processed_condensed_betas.rds"))
    betas <- betas[, colnames(betas) %in% ss$IDAT]
    labels <- as.factor(ss$Invasiveness)
    folds <- readRDS(file.path(dir, "20250125_thyroid_folds.rds"))
    
    # Fold processing
    train_idx <- unlist(folds[-i])
    valid_idx <- folds[[i]]
    
    # Feature selection
    sel_probes <- subset_feature_selection(
        betas = betas,
        batch_size = 10000,
        n_feat = 3000,
        train_idx = train_idx,
        labels = labels,
        fold = i
    )
    
    # Model training
    model <- randomForest(
        x = t(betas[sel_probes, train_idx]),
        y = labels[train_idx],
        ntree = 500,
        importance = TRUE
    )
    
    # Save model
    saveRDS(model, file.path(model_dir, paste0(date, "_model_fold", i, ".rds")))
    
    # PREDICTIONS
    pred_probs <- predict(model, t(betas[sel_probes, valid_idx]), type = "prob")
    pred_class <- predict(model, t(betas[sel_probes, valid_idx]))
    
    # ========================
    # METRICS CALCULATION
    # ========================
    # ROC/AUC
    roc_obj <- pROC::roc(labels[valid_idx], pred_probs[,2])
    auc_value <- pROC::auc(roc_obj)
    
    # Confusion Matrix Metrics
    conf_matrix <- caret::confusionMatrix(
        data = as.factor(pred_class),
        reference = as.factor(labels[valid_idx]),
        positive = levels(labels[valid_idx])[2]
    )
    
    # Compile metrics
    metrics <- data.frame(
        Fold = i,
        AUC = as.numeric(auc_value),
        Accuracy = conf_matrix$overall["Accuracy"],
        Sensitivity = conf_matrix$byClass["Sensitivity"],
        Specificity = conf_matrix$byClass["Specificity"],
        PPV = conf_matrix$byClass["Pos Pred Value"],
        NPV = conf_matrix$byClass["Neg Pred Value"],
        F1 = conf_matrix$byClass["F1"],
        MCC = mltools::mcc(pred_class, labels[valid_idx])
    )
    
    # ========================
    # SAVING RESULTS
    # ========================
    # Save model
    saveRDS(model, file.path(model_dir, paste0(date, "_model_fold", i, ".rds")))
    
    # Save predictions
    preds <- data.frame(
        Fold = i,
        Sample = colnames(betas)[valid_idx],
        TrueLabel = labels[valid_idx],
        PredProb = pred_probs[,2],
        PredClass = pred_class
    )
    saveRDS(preds, file.path(pred_dir, paste0(date, "_preds_fold", i, ".rds")))
    
    # Save metrics
    saveRDS(metrics, file.path(sum_dir, paste0(date, "_metrics_fold", i, ".rds")))
    
    # Save importance
    imp <- data.frame(
        Feature = rownames(model$importance),
        MeanDecreaseAccuracy = model$importance[,"MeanDecreaseAccuracy"],
        MeanDecreaseGini = model$importance[,"MeanDecreaseGini"],
        Fold = i
    )
    saveRDS(imp, file.path(feat_dir, paste0(date, "_importance_fold", i, ".rds")))
    
    return(list(fold = i, status = "complete"))
}

# ========================
# Execution & Aggregation
# ========================
# Run parallel jobs
bplapply(1:10, process_fold, BPPARAM = slurm_config)

# After job completion - Run this separately after all jobs finish
# Aggregate results
all_imp <- list.files(feat_dir, pattern = "_importance_fold", full.names = TRUE) |>
    lapply(readRDS) |>
    bind_rows()

all_preds <- list.files(pred_dir, pattern = "_preds_fold", full.names = TRUE) |>
    lapply(readRDS) |>
    bind_rows()

all_metrics <- list.files(sum_dir, pattern = "_metrics_fold", full.names = TRUE) |>
    lapply(readRDS) |>
    bind_rows()

# Calculate average metrics
avg_metrics <- all_metrics |>
    group_by(Fold) |>
    summarise(
        AUC = mean(AUC),
        Accuracy = mean(Accuracy),
        Sensitivity = mean(Sensitivity),
        Specificity = mean(Specificity),
        F1 = mean(F1),
        MCC = mean(MCC)
    )

# Save final summaries
saveRDS(all_imp, file.path(sum_dir, paste0(date, "_all_importance.rds")))
saveRDS(all_preds, file.path(sum_dir, paste0(date, "_all_predictions.rds")))
write.csv(avg_metrics, file.path(sum_dir, paste0(date, "_average_metrics.csv")))

# Print final metrics
cat("\nFinal Average Metrics:\n")
print(avg_metrics)




