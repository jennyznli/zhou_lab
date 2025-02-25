# ========================
# SHAP Analysis
# ========================
train_idx <- unlist(folds[-1])
test_idx <- folds[[1]]

bSubMostVariable <- function(betas, n=100) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}
betas1 <- t(bSubMostVariable(t(betas)))

train_labels <- factor(ifelse(as.factor(ss$Invasiveness)[train_idx] == "Low", 0, 1),
                       levels = c(0, 1))
train_data <- as.data.frame(betas1[train_idx,])
test_data <- as.data.frame(betas1[test_idx,])

model <- randomForest(
    x = train_data,
    y = train_labels,
    ntree = 500,
    mtry = 2,
    nodesize = 5,
    importance = TRUE
)
unified_model <- randomForest.unify(model, train_data)
shap <- treeshap(unified_model, train_data, interactions = TRUE)
shp <- shapviz(shap, X = train_data)

# most important features
most_imp <- sort(colMeans(abs(shap$shaps)), decreasing = TRUE)
print(head(most_imp))

## PLOTTING
p <- sv_importance(shp,
                   kind = "beeswarm",
                   max_display = 30)
# p <- sv_importance(shp, kind = "both") # if u want both beeswarm and barplot
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_beeswarm.pdf")), height = 8, width = 6)
plot(p)
dev.off()

# specify specific feature interaction
p <- sv_dependence(shp,
                   v = "cg00768487",
                   color_var = "cg11868461")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

# a little less correlated...
p <- sv_dependence(shp,
                   v = "cg00768487",
                   color_var = "cg10579279")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence2.pdf")), height = 3.5, width = 4)
plot(p)
dev.off()

# plots multiple feature interactions
imp_cpgs <- names((most_imp))
p <- sv_dependence(shp, v = imp_cpgs[1], color_var = imp_cpgs[2:10], interactions = TRUE)
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_dependence3.pdf")), height = 10, width = 10)
plot(p)
dev.off()

# interaction plot
# p <- sv_interaction(shp, max_display = 2, size = 3) # max displays
p <- sv_interaction(shp, max_display = 20)
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_interaction.pdf")), height = 20, width = 20)
plot(p)
dev.off()

# force plot 1
p <- sv_force(shp, row_id = "207222790148_R07C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_force1.pdf")), height = 4, width = 7)
plot(p)
dev.off()
p <- sv_force(shp, row_id = "207700160043_R03C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_force2.pdf")), height = 4, width = 7)
plot(p)
dev.off()

# waterfall plot
p <- sv_waterfall(shp, row_id = "207222790148_R07C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_waterfall1.pdf")), height = 7, width = 7)
plot(p)
dev.off()
p <- sv_waterfall(shp, row_id = "207700160043_R03C01")
pdf(file.path(fig_dir, paste0(date, "_thyroid88_shap_waterfall2.pdf")), height = 7, width = 7)
plot(p)
dev.off()


## TEST SHAP CALCULATIONS
# model <- randomForest(
#     x = train_betas,
#     y = as.factor(ss$Invasiveness)[train_idx],
#     ntree = n_trees,
#     importance = TRUE
# )
#
# # IML PACKAGE
# predict_function <- function(model, newdata) {
#     predict(model, newdata, type = "prob")
# }
#
# predictor <- Predictor$new(model = model,
#                            data = as.data.frame(train_betas),
#                            predict.function = predict_function,
#                            y = as.factor(ss$Invasiveness[train_idx])
# )
# shapley <- Shapley$new(predictor, as.data.frame(test_betas))
#
# shap_data <- shapley$results
# shap_data$feature_name <- as.character(shap_data$feature)
# shap_data$feature_value <- as.numeric(sub(".*=", "", shap_data$feature.value))
# feature_values <- as.data.frame(test_betas)  # the data you explained
#
# shap_matrix <- matrix(
#     shap_data$phi,
#     nrow = length(unique(shap_data$observation)),
#     ncol = length(unique(shap_data$feature)),
#     dimnames = list(NULL, unique(shap_data$feature))
# )
#
# sv <- shapviz(
#     shap_values = matrix(shap_data$phi, ncol = ncol(test_betas)),
#     X = feature_values,
#     baseline = mean(predict_function(model, as.data.frame(train_betas)))
# )
# pdf(file.path(fig_dir, paste0(date, "_shap3.pdf")), height = 6, width = 6)
# plot(p)
# dev.off()

# # 1. Bar plot of mean absolute SHAP values (Feature Importance)
# p <- ggplot(shap_data, aes(x = reorder(feature_name, abs(phi)), y = abs(phi))) +
#     geom_bar(stat = "summary", fun = "mean") +
#     coord_flip() +
#     theme_minimal() +
#     labs(x = "CpG Site", y = "Mean |SHAP value|",
#          title = "Feature Importance by Mean |SHAP value|")
#
# # 2. Beeswarm plot (SHAP values vs Feature value)
# p <- ggplot(shap_data, aes(x = phi, y = reorder(feature_name, abs(phi)),
#                       color = feature_value)) +
#     geom_jitter(width = 0, height = 0.2) +
#     scale_color_gradient(low = "blue", high = "red") +
#     theme_minimal() +
#     labs(x = "SHAP value", y = "CpG Site",
#          title = "SHAP Values Distribution by CpG Site",
#          color = "Beta Value")
#
# # 3. Violin plot of SHAP value distributions
# p <- ggplot(shap_data, aes(x = reorder(feature_name, abs(phi)), y = phi)) +
#     geom_violin() +
#     coord_flip() +
#     theme_minimal() +
#     labs(x = "CpG Site", y = "SHAP value",
#          title = "Distribution of SHAP Values by CpG Site")

## FAST SHAP
# shap_values <- fastshap::explain(
#     model,
#     X = as.data.frame(test_betas),
#     pred_wrapper = pred,
#     nsim = 100,
#     shap_only = FALSE
# )
# sv <- shapviz(shap_values)
# p <- sv_importance(sv, show_numbers = TRUE)
# p <- sv_importance(sv, kind = "beeswarm")
# xvars <- c("log_carat", "cut", "color", "clarity")
# p <- sv_dependence(sv, v = xvars)
#
# pdf(file.path(fig_dir, paste0(date, "_treeshap1.pdf")), height = 8, width = 10)
# plot(p)
# dev.off()

## TREE SHAP
# train_labels <- factor(ifelse(as.factor(ss$Invasiveness)[train_idx] == "Low", "0", "1"),
#                        levels = c("0", "1"))
#
# # Combine features and labels into a data frame
# train_data <- cbind(data.frame(target = train_labels), train_betas)
#
# # Train random forest model specifying classification explicitly
# model <- randomForest(
#     target ~ .,
#     data = train_data,
#     ntree = 500,
#     mtry = 2,
#     nodesize = 5,
#     importance = TRUE,
#     type = "classification"  # Explicitly specify classification
# )
#
# unified_model <- randomForest.unify(model, train_data)
# shaps <- treeshap(unified_model, train_data[1:2,])

# # TREE SHAP PACKAGE
# # Install and load required packages
# install.packages("treeshap")
# library(treeshap)
#
# # Your existing data preparation code
# train_idx <- unlist(folds[-1])
# test_idx <- folds[[1]]
#
# # Function to get most variable CpGs
# bSubMostVariable <- function(betas, n=50) {
#     std <- apply(betas, 1, sd, na.rm=TRUE)
#     betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
# }
# betas1 <- t(bSubMostVariable(t(betas)))
# train_betas <- betas1[train_idx,]
# test_betas <- betas1[test_idx,]
#
# train_df <- as.data.frame(train_betas)
# test_df <- as.data.frame(test_betas)
#
# # Set proper column names if they're missing
# colnames(train_df) <- colnames(train_betas)
# colnames(test_df) <- colnames(test_betas)
#
# ranger_model <- ranger(
#     y = as.factor(ss$Invasiveness[train_idx]),
#     x = train_df,
#     num.trees = n_trees,
#     importance = 'impurity'
# )
#
# # Create unified_model object for treeshap
# unified_model <- treeshap::unify(ranger_model, train_df)
#
# unified_model <- treeshap::unify(
#     model = model,
#     X = train_df,
#     target_index = NULL,  # for classification
#     forest = TRUE  # specify it's a random forest
# )
#
#
#
#
# # Calculate SHAP values
# shap_values <- treeshap::treeshap(
#     unified_model,
#     test_df,
#     verbose = FALSE
# )
#
# # Create visualizations
# # 1. Feature Importance Plot
# plot_feature_importance(shap_values)
#
# # 2. SHAP Dependence Plot for specific feature (e.g., first CpG)
# plot_feature_dependence(
#     shap_values,
#     feature = colnames(test_df)[1]
# )
#
# # 3. SHAP Interaction Plot
# plot_interaction(
#     shap_values,
#     feature = colnames(test_df)[1]
# )
#
# # Save plots to PDF
# pdf(file.path(fig_dir, paste0(date, "_treeshap.pdf")), height = 8, width = 10)
# plot_feature_importance(shap_values)
# plot_feature_dependence(shap_values, feature = colnames(test_df)[1])
# dev.off()
#
# # If you want to extract the SHAP values for further analysis:
# shap_df <- as.data.frame(shap_values$shaps)
#
#
#
#
# pred <- function(X.model, newdata) {
#     predict(X.model, newdata, type = "prob")[,2]  # probability of second class
# }
#
# #
# # # Feature importance plot
# # importance_df <- data.frame(
# #     feature = colnames(train_betas),
# #     importance = colMeans(abs(shap_values))
# # )
# #
# # p <- ggplot(importance_df, aes(x = reorder(feature, importance), y = importance)) +
# #     geom_col() +
# #     coord_flip() +
# #     theme_minimal() +
# #     labs(x = "CpG Site", y = "Mean |SHAP value|",
# #          title = "Feature Importance")
# # pdf(file.path(fig_dir, paste0(date, "_shap4.pdf")), height = 6, width = 6)
# # plot(p)
# # dev.off()
#
#
#
#
# predict_function <- function(model, newdata) {
#     predict(model, newdata, type = "prob")
# }
# baseline <- mean(predict_function(model, as.data.frame(train_betas)))
#
#
# unified_model <- randomForest.unify(model, train_betas)
# shap_values <- treeshap(unified_model, t(train_betas[,1:10]))
# #
# # # train_idx <- unlist(folds[-i])
# # sv <- shapviz(shap_values, X = t(train_betas[,1:10]))
# # sv_importance(sv)



