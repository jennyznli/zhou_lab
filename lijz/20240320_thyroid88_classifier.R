x<-c("sesame", "SummarizedExperiment","impute", "RSNNS", "e1071", "readxl", "readr", "dplyr", "randomForest", "caTools", "tidyverse", "caret", "ISLR", "pROC")
lapply(x, require, character.only = TRUE)

# COLLAPSE TO PREFIX_______________________________________________________________
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
betas = readRDS(here("data", "20240320_thyroid136_betas.rds")) #937690

betas2 <- betasCollapseToPfx(betas) #  930659
saveRDS(betas2, "~/Documents/HPC_share/thyroid/20240320_thyroid136_betas_condensed.rds")

# PREPROCESSING IN HPC_______________________________________________________________

meta = read_excel("/home/lijz/20230818_thyroid_cancer/20231102_thyroid_master.xlsx")
betas = readRDS("/home/lijz/20230818_thyroid_cancer/20240320_thyroid136_betas_condensed.rds")
meta = meta %>% dplyr::filter(Include_In_Analysis == 1)
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

# CROSS VALIDATION  _______________________________________________________________

ss = read_excel("~/Documents/HPC_share/thyroid/ss/20231102_thyroid_master.xlsx")
ss = ss %>% dplyr::filter(INCLUDE_IN_ANALYSIS == "1")

# 3k most important, cleaned betas
betas = readRDS(file = "~/Documents/HPC_share/thyroid/mtx.rds")
betas = betas[rownames(betas) %in% ss$IDAT,]

folds = createFolds(ss$Invasiveness, k = 10, list = TRUE, returnTrain = TRUE)
predictions = lapply(folds, function(train) {
    model = randomForest(t(betas[,train]), as.factor(ss$Invasiveness[train]))
    test = setdiff(seq_along(ss$Invasiveness), train)
    predict(model, t(betas[,test]))
})

# predictions = readRDS("~/Documents/HPC_share/thyroid/20240330_predictions.rds")
df = do.call(bind_rows, lapply(names(predictions), function(x) {a = predictions[[x]]; data.frame(Fold=x, prediction=a, IDAT=names(a))}))

df$actual = ss$Invasiveness[match(df$IDAT, ss$IDAT)]
df$prediction = factor(df$prediction, levels=sort(unique(ss$Invasiveness)))
df$actual = factor(df$actual, levels=sort(unique(ss$Invasiveness)))

df1 = df |> group_by(prediction) |> summarize(tp = sum(prediction == actual), fp = sum(prediction != actual))
df2 = df |> group_by(actual) |> summarize(fn = sum(prediction != actual))

df1$fn = df2$fn[match(df1$prediction, df2$actual)]
df1$f1 = with(df1, tp/(tp+(fp+fn)/2))

confusion = t(table(df$actual, df$prediction))
cm[[2]]

cm <- confusionMatrix(df$prediction, df$actual)

# accuracy <- c()
# precision <- c()
# recall <- c()
# f1 <- c()

accuracy <- c(accuracy, cm$overall['Accuracy'])
precision <- c(precision, cm$byClass['Pos Pred Value'])
recall <- c(recall, cm$byClass['Sensitivity'])
f1 <- c(f1, cm$byClass['F1'])

cat("Mean Accuracy: ", mean_accuracy, "\n")
cat("Mean Precision: ", mean_precision, "\n")
cat("Mean Recall: ", mean_recall, "\n")
cat("Mean F1 Score: ", mean_f1, "\n")

# FEATURE SELECTION _______________________________________________________________

imp = read_tsv(col_names = FALSE, "~/Documents/HPC_share/thyroid/20240610_features.tsv")
colnames(imp) = c("Model", "CpG", "Accuracy", "Gini")
imp = imp %>% arrange(desc(Gini))

pos = imp %>% filter(Gini > 0)
rs <- grep('rs', imp$CpG)
top3k = pos[1:3000,]$CpG

###

ss = read_excel("~/Documents/HPC_share/thyroid/20231102_thyroid_master.xlsx")
ss = ss %>% dplyr::filter(INCLUDE_IN_ANALYSIS == "1")

betas = readRDS(file = "~/Documents/HPC_share/thyroid/20240320_thyroid136_betas.rds") #937690
# 3k most important, cleaned betas
betas = readRDS(file = "~/Documents/HPC_share/thyroid/mtx.rds")

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
# imputeRowMean <- function(mtx) {
#     k <- which(is.na(mtx), arr.ind=TRUE)
#     mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
#     mtx
# }
# betas <- cleanMatrixForClusterW(betas) %>%
#     imputeRowMean(.)

betas <- cleanMatrixForClusterW(betas) %>%
    impute.knn(.)

# bSubMostVariable <- function(betas, n=2000) {
#     std <- apply(betas, 1, sd, na.rm=TRUE)
#     betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
# }
# mtx = bSubMostVariable(t(betas), n=3000)

mtx = betas$data

mtx1 = mtx[rownames(mtx) %in% top3k, colnames(mtx) %in% ss$IDAT]

model = randomForest(t(mtx1), as.factor(ss$Invasiveness), xtest = t(mtx1), ytest = as.factor(ss$Invasiveness), keep.forest = TRUE)
# model$test$votes
# model$test$predicted
rf_prob <- predict(model, t(mtx1), type = "prob")
rf_pred <- predict(model, t(mtx1))

# make mdoels
# rf = randomForest(t(mtx1), as.factor(ss$Invasiveness))
#
# mlp = mlp(t(mtx1), as.factor(ss$Invasiveness))
#
# svm = svm(as.factor(ss$Invasiveness) ~ ., data = t(mtx1), kernel = "linear", cost = 10, scale = FALSE)
# print(svm)

# plot(svm, t(mtx1))

## cross validation

model <- nnet((as.factor(ss$Invasiveness)) ~ ., data = t(mtx1), size = 5, decay = 0.1, maxit = 1000, MaxNWts=84581)
# predict(model, betas)
saveRDS(model, "~/Documents/HPC_share/thyroid/20240612_nnet_model.rds")


set.seed(123)
mtx1 = t(mtx1)
mtx1$Invasiveness = ss$Invasiveness
folds = createFolds(ss$Invasiveness, k = 10, list = TRUE, returnTrain = TRUE)
predictions = lapply(folds, function(train) {
    # model = randomForest(t(mtx1[,train]), as.factor(ss$Invasiveness[train]))
    # model = svm(as.factor(ss$Invasiveness[train]) ~ ., data = t(mtx1[,train]), kernel = "linear", cost = 10, scale = FALSE)
    model <- nnet((as.factor(ss$Invasiveness[train])) ~ ., data = mtx1[,train], size = 5, decay = 0.1, maxit = 200)
    test = setdiff(seq_along(ss$Invasiveness), train)
    predict(model, t(mtx1[,test]))
    predict(model, t(mtx1[,test]), type = "prob")[, 2]
})

# saveRDS(predictions, "~/Documents/HPC_share/thyroid/20240612_predictions.rds")
##
df = do.call(bind_rows, lapply(names(predictions), function(x) {a = predictions[[x]]; data.frame(Fold=x, prediction=a, IDAT=names(a))}))

df$actual = ss$Invasiveness[match(df$IDAT, ss$IDAT)]
df$prediction = factor(df$prediction, levels=sort(unique(ss$Invasiveness)))
df$actual = factor(df$actual, levels=sort(unique(ss$Invasiveness)))

df1 = df |> group_by(prediction) |> summarize(tp = sum(prediction == actual), fp = sum(prediction != actual))
df2 = df |> group_by(actual) |> summarize(fn = sum(prediction != actual))

df1$fn = df2$fn[match(df1$prediction, df2$actual)]
df1$f1 = with(df1, tp/(tp+(fp+fn)/2))


# METRICS
cm <- caret::confusionMatrix(as.factor(ss$Invasiveness), rf_pred)

accuracy <- cm$overall["Accuracy"]
precision <- cm$byClass["Pos Pred Value"]
recall <- cm$byClass["Sensitivity"]
f1 <- 2 * (precision * recall) / (precision + recall)

cat("Accuracy: ", accuracy, "\n")
cat("Precision: ", precision, "\n")
cat("Recall: ", recall, "\n")
cat("F1 Score: ", f1, "\n")


rocCurve <- roc(response = as.numeric(testData$Species), predictor = probs, levels = rev(levels(testData$Species)))
plot(rocCurve, main = "ROC Curve")
aucValue <- auc(rocCurve)
cat("AUC: ", aucValue, "\n")



precision <- posPredValue(df$prediction, df$actual, positive="High")
recall <- sensitivity(df$prediction, df$actual, positive="High")
f1 <- (2 * precision * recall) / (precision + recall)

roc_object <- pROC::multiclass.roc(df$actual, as.numeric(df$prediction))
auc(roc_object)

#
# accuracy <- c()
# precision <- c()
# recall <- c()
# f1 <- c()
#
# accuracy <- c(accuracy, cm$overall['Accuracy'])
# precision <- c(precision, cm$byClass['Pos Pred Value'])
# recall <- c(recall, cm$byClass['Sensitivity'])
# f1 <- c(f1, cm$byClass['F1'])
#
# print(" Accuracy: ", accuracy, "\n")
# print(" Precision: ", precision, "\n")
# print(" Recall: ", recall, "\n")
# print(" F1 Score: ", f1, "\n")

## WITH ALL BETAS
# Reference
# Prediction High Low
# High   33  16
# Low     2   9
# performed best.... uhhh

# WITH 3K MOST IMPORTANT BETAS
# Reference
# Prediction High Low
# High   33  16
# Low     2   9

### NEW FINALIZED RF VERSION

# Load necessary libraries
library(caret)
library(randomForest)
library(pROC)

# # Example dataset: Using the Iris dataset for binary classification
# data(iris)
#
# # Convert it to a binary classification problem
# iris <- iris[iris$Species != "setosa", ]
# iris$Species <- factor(iris$Species)

# Define the cross-validation method
train_control <- trainControl(method = "cv", number = 10, savePredictions = TRUE, classProbs = TRUE, summaryFunction = twoClassSummary)

# Train a Random Forest model with cross-validation
set.seed(42)
rfModel <- train(y = as.factor(ss$Invasiveness), x = t(mtx1), method = "rf", trControl = train_control, metric = "ROC")

# Predictions from cross-validation
predictions <- rfModel$pred
probs <- predictions[,"High"]

# Calculate overall metrics from cross-validation results
confMat <- confusionMatrix(predictions$pred, predictions$obs)

# Calculate accuracy
accuracy <- confMat$overall["Accuracy"]

# Calculate precision, recall, and F1 score
precision <- confMat$byClass["Pos Pred Value"]
recall <- confMat$byClass["Sensitivity"]
f1 <- 2 * (precision * recall) / (precision + recall)
npv <- confMat$byClass["Neg Pred Value"]
ppv <- confMat$byClass["Pos Pred Value"]

# Print metrics
cat("Accuracy: ", accuracy, "\n")
cat("Precision: ", precision, "\n")
cat("Recall: ", recall, "\n")
cat("F1 Score: ", f1, "\n")
cat("NPV: ", npv, "\n")
cat("PPV: ", ppv, "\n")

# Plot ROC curve and calculate AUC
rocCurve <- roc(response = predictions$obs, predictor = probs, levels = rev(levels(predictions$obs)))
plot(rocCurve, main = "ROC Curve")
aucValue <- auc(rocCurve)

# Print AUC
cat("AUC: ", aucValue, "\n")
