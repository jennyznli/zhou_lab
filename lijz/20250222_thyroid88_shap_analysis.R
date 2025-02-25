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




