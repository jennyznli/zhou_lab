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

train_data_ranked <- apply(train_data, 1, rank, ties.method = "average", na.last=TRUE)
train_data_ranked <- t(train_data_ranked)

# model_time_ranked <- system.time({
#     rf_model_ranked <- randomForest(train_data_ranked, train_labels, ntree = 500)
# }) 
rf_model_ranked <- randomForest(train_data_ranked, train_labels, ntree = 500)

print(rf_model_ranked)
# OOB estimate of  error rate: 8.79%

# predict_time <- system.time({
#     test_predictions <- predict(rf_model, test_data)
# })
test_data_ranked <- apply(test_data, 1, rank, ties.method = "average", na.last=TRUE)
test_data_ranked <- t(test_data_ranked)
test_predictions_ranked <- predict(rf_model_ranked, test_data_ranked)

prediction_probs_ranked <- predict(rf_model_ranked, newdata = test_data_ranked, type = "prob", na.omit = TRUE)
probabilities_ranked <- prediction_probs_ranked[cbind(1:length(test_predictions_ranked), as.integer(test_predictions_ranked))]

test_labels <- as.character(test_labels)
test_accuracy_ranked <- sum(test_predictions_ranked == test_labels) / length(test_labels)
print(paste("Test Accuracy:", test_accuracy_ranked))
# "Test Accuracy: 0.905674281503316"

IDs <- read_tsv("~/Documents/HPC_share/TCGA/TCGA6786_CancerType33_230227.tsv")

test_GSMs <- IDs[test_indices, "Sample_ID"]

rf_predictions_ranked <- data.frame(
    Sample_ID = c(test_GSMs),
    predicted = c(test_predictions_ranked),
    actual = c(test_labels),
    probability_score = c(probabilities_ranked)
)

saveRDS(rf_model_ranked, "~/Documents/HPC_share/TCGA/20230808_rfc_test_ntree=500_ranked.rds")
write_tsv(rf_predictions_ranked, "~/Documents/HPC_share/TCGA/20230808_rf_predictions_ranked.tsv")

# print(model_time["elapsed"])
# print(predict_time["elapsed"])

tcga_comparison <- data.frame(
    Sample_ID = c(test_GSMs),
    actual = c(test_labels),
    Ranked_Prediction = c(test_predictions_ranked),
    Normal_Prediction = c(test_predictions)
)
write_tsv(tcga_comparison, "~/Documents/HPC_share/TCGA/20230808_rf_prediction_comparison.tsv")

tcga_class_accuracies <- data.frame(
    Ranked_Class_Accuracy = rf_model_ranked$confusion[,'class.error'],
    Normal_Class_Accuracy = rf_model$confusion[,'class.error']
)
tcga_class_accuracies$Ranked_Order <- rownames(tcga_class_accuracies)[order(tcga_class_accuracies$Ranked_Class_Accuracy, decreasing = TRUE)]
tcga_class_accuracies$Normal_Order <- rownames(tcga_class_accuracies)[order(tcga_class_accuracies$Normal_Class_Accuracy, decreasing = TRUE)]
write_tsv(tcga_class_accuracies, "~/Documents/HPC_share/TCGA/20230808_rf_class_accuracies_comparison.tsv")
