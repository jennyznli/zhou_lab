x<-c("sesame", "RSNNS", "readxl", "readr", "dplyr", "caret", "randomForest", "impute", "tidyverse", "caret", "e1071", "pROC")
lapply(x, require, character.only = TRUE)

# PREPROCESSING __________________________________________________

betas = readRDS(file = "~/Documents/HPC_share/thyroid/20240320_thyroid136_betas_condensed.rds") #930659 136
ss = read_excel("~/Documents/HPC_share/thyroid/20231102_thyroid_master.xlsx") #136 34
ss = ss %>% dplyr::filter(INCLUDE_IN_ANALYSIS == "1") #88 34
betas = betas[, colnames(betas) %in% ss$IDAT] #930659   88

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
betas <- cleanMatrixForClusterW(betas) 
# Filter rows with >0.50 missingness and columns with >0.50 missingness.
# Before:  930659 rows and  88 columns.
# After:  855437 rows and  88 columns.

# knn imputation 
betas1 <- impute.knn(betas, k = 10, rng.seed=1234)
betas_knn <- betas1$data #855437 88
betas_knn = t(betas_knn) #88 855437

# MOST IMPORTANT  __________________________________________________

imp = read_tsv(col_names = FALSE, "~/Documents/HPC_share/thyroid/features/20240610_features.tsv") #855437 4
colnames(imp) = c("Model", "CpG", "Accuracy", "Gini")
imp = imp %>% arrange(desc(Accuracy))

pos = imp %>% filter(Accuracy > 0)
rs <- grep('rs', pos$CpG)
# print(length(rs))
mostImp = pos[1:10000,]$CpG

mtx = betas_knn
mtx = mtx[rownames(mtx) %in% ss$IDAT, colnames(mtx) %in% mostImp]

##  MODEL ITERATIONS____________________________________________________________
# Set a global seed for overall reproducibility
set.seed(42)

# Prepare the data
invasiveness <- as.factor(ss$Invasiveness)
data <- data.frame(mtx, invasiveness)

# Create cross-validation folds
folds <- createFolds(invasiveness, k = 10, list = TRUE, returnTrain = FALSE)

# Initialize vectors to store results
accuracy_vec <- precision_vec <- recall_vec <- f1_vec <- auc_vec <- npv_vec <- ppv_vec <- c()
all_predictions <- data.frame()

# Perform cross-validation
for (i in 1:length(folds)) {
    # Set a seed for each fold to ensure reproducibility within the loop
    set.seed(42 + i)
    
    # Split the data
    test_indices <- folds[[i]]
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]
    
    # Train the model
    # Set seed and use other parameters to ensure reproducibility
    rf_model <- randomForest(invasiveness ~ ., data = train_data, 
                             ntree = 500,  # Fix number of trees
                             mtry = sqrt(ncol(train_data) - 1),  # Fix mtry
                             importance = TRUE,
                             seed = 42 + i)  # Set seed for randomForest
    
    # Make predictions
    predictions <- predict(rf_model, test_data)
    probs <- predict(rf_model, test_data, type = "prob")[, 2]
    
    # Calculate metrics
    conf_mat <- confusionMatrix(predictions, test_data$invasiveness, positive = levels(invasiveness)[2])
    roc_curve <- roc(response = test_data$invasiveness, predictor = probs)
    
    # Store metrics
    accuracy_vec <- c(accuracy_vec, conf_mat$overall["Accuracy"])
    precision_vec <- c(precision_vec, conf_mat$byClass["Pos Pred Value"])
    recall_vec <- c(recall_vec, conf_mat$byClass["Sensitivity"])
    f1_vec <- c(f1_vec, conf_mat$byClass["F1"])
    auc_vec <- c(auc_vec, auc(roc_curve))
    npv_vec <- c(npv_vec, conf_mat$byClass["Neg Pred Value"])
    ppv_vec <- c(ppv_vec, conf_mat$byClass["Pos Pred Value"])
    
    # Store predictions
    all_predictions <- rbind(all_predictions, data.frame(predictions, actual = test_data$invasiveness))
}

# Calculate and print results
results <- data.frame(
    Metric = c("Accuracy", "Precision", "Recall", "F1", "AUC", "NPV", "PPV"),
    Mean = c(mean(accuracy_vec), mean(precision_vec), mean(recall_vec), mean(f1_vec), mean(auc_vec), mean(npv_vec), mean(ppv_vec)),
    SD = c(sd(accuracy_vec), sd(precision_vec), sd(recall_vec), sd(f1_vec), sd(auc_vec), sd(npv_vec), sd(ppv_vec))
)
print(results)

# Final confusion matrix
final_cm <- confusionMatrix(as.factor(all_predictions$predictions), as.factor(all_predictions$actual), positive = levels(invasiveness)[2])
print(final_cm)

# Save predictions
write.table(all_predictions, "~/Documents/HPC_share/thyroid/20241010_thyroid88_predictions.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# Set seed before plotting for reproducible jitter if used
set.seed(42)
# Plot ROC curve (using the last fold for illustration)
plot(roc_curve, main = paste("ROC Curve (AUC =", round(auc(roc_curve), 2), ")"))

## VISUALIZATION _______________________________________________________________

# models <- c("randomForest", "SVM", "MLP")
methods <- c("2k most variable", "2k most important", "2k diff meth")
metrics <- c("Accuracy", "Precision", "Recall", "F1 Score", "AUC")
# 2k most variable
values <- matrix(c(
    #rf   #svm   #mlp
    0.8636111, 0.865, 0.8438889,    # Accuracy
    0.9016667, 0.9016667, 0.8516667,      # Precision
    0.7833333, 0.7833333, 0.7583333,    # Recall
    0.7997619, 0.7997619, 0.7583333,      # F1 Score
    0.9391667, 0.9308333, 0.9002778     # AUC
), nrow = 5, byrow = TRUE)

# 2k most imp 
values <- matrix(c(
    #rf   #svm   #mlp
    0.8911111, 0.8886111, 0.8386111,    # Accuracy
    0.91, 0.935, NaN,      # Precision
    0.91, 0.8166667, 0.7083333,    # Recall
    0.8388095, 0.8330952, NaN,      # F1 Score
    0.9616667, 0.945, 0.8627778     # AUC
), nrow = 5, byrow = TRUE)

# 2k diff meth
values <- matrix(c(
    #rf   #svm   #mlp
    0.8525, 0.8761111, 0.855,    # Accuracy
    0.8766667, 0.9016667, 0.8766667,      # Precision
    0.7833333, 0.8166667, 0.7833333,    # Recall
    0.7854762, 0.8197619, 0.7890476,      # F1 Score
    0.9291667, 0.9258333, 0.9327778     # AUC
), nrow = 5, byrow = TRUE)

df <- data.frame(model = rep(models, each = 5),
                 method = rep(methods, each = 5),
                 metric = rep(metrics),
                 value = as.vector(values))

# Plot grouped bar chart
pdf('~/Documents/HPC_share/thyroid/20240617_thyroid88_classifier_2k_mostvar.pdf', family="ArialMT", width=8, height=5, onefile=FALSE)
ggplot(df, aes(x = model, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ metric, nrow=1) +
    coord_cartesian(ylim=c(0.5, 1)) + 
    labs(title = "2k MostVar Performance Metrics Comparison",
         x = "Model", y = "Metric Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf('~/Documents/HPC_share/thyroid/20240617_thyroid88_classifier_2k_diffMeth.pdf', family="ArialMT", width=8, height=5, onefile=FALSE)
ggplot(df, aes(x = model, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ metric, nrow=1) +
    coord_cartesian(ylim=c(0.5, 1)) + 
    labs(title = "2k DiffMeth Performance Metrics Comparison",
         x = "Model", y = "Metric Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf('~/Documents/HPC_share/thyroid/20240617_thyroid88_classifier_2k_hiImp.pdf', family="ArialMT", width=8, height=5, onefile=FALSE)
ggplot(df, aes(x = model, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ metric, nrow=1) +
    coord_cartesian(ylim=c(0.5, 1)) + 
    labs(title = "highImp Performance Metrics Comparison",
         x = "Model", y = "Metric Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## TEST ON GEO DATA _________________________________________________
betas_geo = readRDS("~/Documents/HPC_share/thyroid/20240527_GEO_DICER1_betas.rds")
dim(betas_geo) # 866553     17
epic_probes <- rownames(betas_geo)
betas_knn_epic <- betas_knn[rownames(betas_knn) %in% epic_probes,]
# [1] 689175     88

## feature selection
imp = read_tsv(col_names = FALSE, "~/Documents/HPC_share/thyroid/20240610_features.tsv")
colnames(imp) = c("Model", "CpG", "Accuracy", "Gini")
imp = imp %>% arrange(desc(Accuracy))
pos = imp %>% filter(Accuracy > 0)
rs <- grep('rs', pos$CpG)
print(length(rs))
pos = pos[-rs,]
mostImp = pos[1:4000,]$CpG
sum(colnames(betas_knn_epic) %in% mostImp)

betas_knn_epic = t(betas_knn_epic)
mtx = betas_knn_epic[rownames(betas_knn_epic) %in% ss$IDAT, colnames(betas_knn_epic) %in% mostImp]

set.seed(500)
model = randomForest(mtx, as.factor(ss$Invasiveness))


# cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
#     cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
#                 f_row, f_col))
#     cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
#     namtx = is.na(mtx)
#     good_row = rowSums(namtx) <= ncol(mtx) * f_row
#     good_col = colSums(namtx) <= nrow(mtx) * f_col
#     cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
#     mtx[good_row, good_col]
# }
# betas_geo <- cleanMatrixForClusterW(betas_geo)

# knn imputation 
betas_geo1 <- impute.knn(betas_geo, k = 10, rng.seed=1234)
betas_geo_knn <- betas_geo1$data

betas_geo_test = betas_geo_knn[rownames(betas_geo_knn) %in% colnames(mtx),]
predictions <- predict(model, t(betas_geo_test), type = "response")
probs <- predict(model, t(betas_geo_test), type = "prob")[, 2]

conf_mat <- caret::confusionMatrix(as.factor(predictions), as.factor(invasiveness), positive = levels(invasiveness)[2])
print(conf_mat)

# Calculate metrics
accuracy <- conf_mat$overall["Accuracy"]
precision <- conf_mat$byClass["Pos Pred Value"]
recall <- conf_mat$byClass["Sensitivity"]
f1 <- 2 * (precision * recall) / (precision + recall)
npv <- conf_mat$byClass["Neg Pred Value"]
ppv <- conf_mat$byClass["Pos Pred Value"]

# Calculate AUC
roc_curve <- roc(response = test_data$invasiveness, predictor = as.numeric(probs), levels = rev(levels(test_data$invasiveness)))
auc_value <- auc(roc_curve)

x = list(accuracy = accuracy, precision = precision, recall = recall, f1 = f1, npv = npv, ppv = ppv, auc = auc_value)
print(x)

## ______________________-

# Load necessary libraries
library(randomForest)
library(caret)

# Define training control using caret's trainControl
train_control <- trainControl(
    method = "cv",        # Using cross-validation
    number = 10,          # Number of folds in the cross-validation
    savePredictions = "final", # Save out-of-fold predictions
    classProbs = TRUE    # Store class probabilities for ROC analysis, etc.
)

# Train the model using caret's train function with a Random Forest model
model <- train(
    Species ~ ., 
    data = df,
    method = "rf",       # Random Forest model
    trControl = train_control,
    tuneLength = 3       # Tuning the number of trees in the random forest
)

# Print the model details to check the chosen parameters and performance
print(model)

# Access the out-of-fold predictions
oof_predictions <- model$pred

# Calculate out-of-fold accuracy
accuracy <- sum(oof_predictions$obs == oof_predictions$pred) / nrow(oof_predictions)
print(paste("Out-of-Fold Accuracy:", accuracy))

# You can also perform more detailed performance analysis using the confusion matrix, etc.
confusionMatrix(oof_predictions$pred, oof_predictions$obs)




