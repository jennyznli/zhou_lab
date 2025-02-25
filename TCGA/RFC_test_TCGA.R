library(readr)
library(randomForest)

df <- readRDS("~/Documents/HPC_share/TCGA/TCGA_10k_CT33_Betas.rds")
labels <- t(read_tsv("~/Documents/HPC_share/TCGA/TCGA6786_CancerType33_230227.tsv")$label)

set.seed(123)

for (col in colnames(df)) {
  mean_value <- mean(df[, col], na.rm = TRUE)
  
  df[is.na(df[, col]), col] <- mean_value
}

indices <- sample(1:nrow(df))  
train_indices <- indices[1:round(0.8 * nrow(df))] 
test_indices <- indices[(round(0.8 * nrow(df)) + 1):nrow(df)]

train_data <- df[train_indices, ]
train_labels <- as.factor(labels[train_indices])
test_data <- df[test_indices, ]
test_labels <- as.factor(labels[test_indices])

# model_time <- system.time({
#   rf_model <- randomForest(train_data, train_labels, ntree = 500)
# }) 
rf_model <- randomForest(train_data, train_labels, ntree = 500)

print(rf_model)
# OOB estimate of  error rate: 8.95%

# predict_time <- system.time({
#   test_predictions <- predict(rf_model, test_data)
# })
test_predictions <- predict(rf_model, test_data)

prediction_probs <- predict(rf_model, newdata = test_data, type = "prob", na.omit = TRUE)
probabilities <- prediction_probs[cbind(1:length(test_predictions), as.integer(test_predictions))]

test_labels <- as.character(test_labels)
test_accuracy <- sum(test_predictions == test_labels) / length(test_labels)
print(paste("Test Accuracy:", test_accuracy))
#   "Test Accuracy: 0.905674281503316"

IDs <- read_tsv("~/Documents/HPC_share/TCGA/TCGA6786_CancerType33_230227.tsv")

test_GSMs <- IDs[test_indices, "Sample_ID"]
rf_predictions <- data.frame(
  Sample_ID = c(test_GSMs),
  predicted = c(test_predictions),
  actual = c(test_labels),
  probability_score = c(probabilities)
)
saveRDS(rf_model, "~/Documents/HPC_share/TCGA/20230808_rfc_test_ntree=500.rds")
write_tsv(rf_predictions, "~/Documents/HPC_share/TCGA/20230808_rf_predictions.tsv")

# print(model_time["elapsed"])
# print(predict_time["elapsed"])