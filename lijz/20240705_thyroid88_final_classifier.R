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

invasiveness = as.factor(ss$Invasiveness)
data <- data.frame(mtx, invasiveness)

# Create cross-validation folds
num_folds <- 10
set.seed(42)
folds <- createFolds(invasiveness, k = num_folds, list = TRUE, returnTrain = TRUE)

# Initialize a function to perform cross-validation and calculate metrics
evaluate_model <- function(train_data, test_data, model_func) {
    # Train the model
    model <- model_func(train_data)
    
    predictions <- predict(model, test_data, type = "response")
    probs <- predict(model, test_data, type = "prob")[, 2]
    # 
    # # Make predictions on the test set
    # if (inherits(model, "randomForest")) {
    #     
    # } else if (inherits(model, "nnet")) {
    #     predictions <- predict(model, test_data, type = "class")
    #     probs <- predict(model, test_data, type = "raw")
    # } else {
    #     predictions <- predict(model, test_data)
    #     probs <- attr(predict(model, test_data, probability = TRUE), "probabilities")[, 2]
    # }
    
    # Calculate confusion matrix
    conf_mat <- caret::confusionMatrix(as.factor(predictions), as.factor(test_data$invasiveness), positive = levels(test_data$invasiveness)[2])
    # print(conf_mat)
    
    predictions <- data.frame(predictions, test_data$invasiveness)
    print(dim(predictions))
    
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
    
    return(list(accuracy = accuracy, precision = precision, recall = recall, f1 = f1, npv = npv, ppv = ppv, auc = auc_value, predictions = predictions))
}

# Define model functions
rf_model_func <- function(data) {
    set.seed(1234)
    randomForest(invasiveness ~ ., data = data)
}

models <- list(randomForest = rf_model_func)

results <- list()  
pred <- matrix(ncol=2, nrow=0)
for (model_name in names(models)) {
    model_func <- models[[model_name]]
    
    accuracy_vec <- c()
    precision_vec <- c()
    recall_vec <- c()
    f1_vec <- c()
    auc_vec <- c()
    npv_vec <- c()
    ppv_vec <- c()
    
    for (i in 1:num_folds) {
        # Split the data into training and testing sets
        train_indices <- folds[[i]]
        train_data <- data[train_indices, ]
        test_data <- data[-train_indices, ]
        
        # Evaluate the model
        metrics <- evaluate_model(train_data, test_data, model_func)
        
        # Store the metrics
        accuracy_vec <- c(accuracy_vec, metrics$accuracy)
        precision_vec <- c(precision_vec, metrics$precision)
        recall_vec <- c(recall_vec, metrics$recall)
        f1_vec <- c(f1_vec, metrics$f1)
        auc_vec <- c(auc_vec, metrics$auc)
        npv_vec <- c(npv_vec, metrics$npv)
        ppv_vec <- c(ppv_vec, metrics$ppv)
        pred <- rbind(pred, metrics$predictions)
        print(dim(pred))
    }
    
    # Calculate mean and standard deviation of metrics
    results[[model_name]] <- 
        data.frame(
            Accuracy = c(mean(accuracy_vec), mean(precision_vec), mean(recall_vec), mean(f1_vec), mean(npv_vec), mean(ppv_vec), mean(auc_vec)),
            SD = c(sd(accuracy_vec), sd(precision_vec), sd(recall_vec), sd(f1_vec), sd(npv_vec), sd(ppv_vec), sd(auc_vec))
        )
    rownames(results[[model_name]]) <- c("Accuracy", "Precision", "Recall", "F1", "AUC", "NPV", "PPV")

}

# final confusion matrix
cm <- caret::confusionMatrix(as.factor(pred$predictions), as.factor(pred$test_data.invasiveness), positive = levels(invasiveness)[2])
print(cm)
results[[1]]

pred1 = sort(rownames(pred), pred)
pred1 = pred[order(row.names(pred)), ]

write_tsv(as.data.frame(pred1), "~/Documents/HPC_share/thyroid/20240624_thyroid88_predictions.tsv")

## ROC AUC CURVE _______________________________________________________________

# Train Random Forest Model
rf_model <- randomForest(Species ~ ., data = train_data, ntree = 100)

# Predict probabilities on test set
rf_prob <- predict(rf_model, test_data, type = "prob")[,2]

# Calculate ROC and AUC
roc_obj <- roc(test_data$Species, rf_prob)

# Plot ROC Curve
ggroc(roc_obj) +
    ggtitle(paste("ROC Curve (AUC =", round(auc(roc_obj), 2), ")")) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    theme_minimal()

# Print AUC value
print(paste("AUC:", auc(roc_obj)))



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




