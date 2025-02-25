

se <- readRDS("~/Documents/HPC_share/CAP/20221210_Capper_Reference_betas_SE.rds")
# tse <- readRDS("~/Documents/HPC_share/CAP/20221010_Capper_Prospective_betas_SE.rds")

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

samples <- match(colnames(se),colData(se)$SampleID)
labels <- factor(colData(se)[samples,"Capper_Methylation_Class"])


dir <- "~/Documents/HPC_share/chromHMM" #feature dir 
features <- list.files(dir)
# out_dir <- "~/20230812_BTC_featureTesting/" #file output dir
# importance_log <- paste0(out_dir,"20230812_importance.tsv") #log for feature importance for each model
# oob_log <- paste0(out_dir,"/20230812_oob.tsv") #log for oob error rates for each model


# out <- paste0("20230812_RFC_",sub("\\.gz","",f),"_Model.rds") #output model file name
aggregateMeth <- function(feature_path,betas,...) {
    feat <- KYCG_loadDBs(feature_path)
    names(feat) <- vapply(feat,function(x) attr(x,"dbname"),character(1))
    ft_mx <- dbStats(betas=betas,databases=feat,...) # what is this? ...
    na_sum <- colSums(is.na(ft_mx))
    if (sum(na_sum) != 0) {
        remove <- which(na_sum == nrow(ft_mx)) #sometimes no cgs from features in betas mtx, remove
        ft_mx <- ft_mx[,-remove]
    }
    ft_mx
}
# x <- read.delim(path)
# y <- read.delim(path)
# path <- "~/Documents/HPC_share/chromHMM/ChromHMM.20220303"
path <- "~/Documents/HPC_share/TFBSrm.20221005"
feat <- KYCG_loadDBs(path) # separates each chromhmm state with its associated cpgs
names(feat) <- vapply(feat,function(x) attr(x,"dbname"),character(1)) # names them 
ft_mx <- dbStats(betas=mtx,databases=feat)
saveRDS(ft_mx, "~/Documents/HPC_share/chromHMM/20230814_tfbs_ft_mx.rds")
ft_mx <- readRDS("~/Documents/HPC_share/chromHMM/20230713_tfbs_ft_mx.rds")

na_sum <- colSums(is.na(ft_mx))
if (sum(na_sum) != 0) {
    remove <- which(na_sum == nrow(ft_mx)) #sometimes no cgs from features in betas mtx, remove 
    ft_mx <- ft_mx[,-remove]
}

# m <- aggregateMeth(feature_path=path,betas=mtx) #create aggregated meth betas

if (ncol(ft_mx) == 0) next
model <- randomForest(
    x=ft_mx,
    y=labels,
    ntree=500,
    importance=TRUE)
# chromhmm OOB estimate of  error rate: 26.94%

saveRDS(model, "~/Documents/HPC_share/chromHMM/20230814_chromhmm_model.rds")
saveRDS(model, "~/Documents/HPC_share/chromHMM/20230814_tfbs_model.rds")
model <- readRDS("~/Documents/HPC_share/chromHMM/20230814_tfbs_model.rds")

importance <- data.frame(
    Model=rep(path,nrow(model$importance)),
    Feature=rownames(model$importance),
    MeanDecreaseAccuracy=model$importance[,"MeanDecreaseAccuracy"],
    MeanDecreaeGini=model$importance[,"MeanDecreaseGini"]
) %>% arrange(desc(MeanDecreaseAccuracy))
# write_tsv(oob,file=oob_log,append = TRUE)
# write_tsv(importance,file = importance_log,append = TRUE)
print(paste("Done with ",f))

# OOB: 0.4055694
oob <- data.frame(
    model=path,
    oob=model$err.rate[nrow(model$err.rate),"OOB"]
)

# for (f in features) {
#     path <- paste0(dir,f) #feature file path
#     out <- paste0("20230812_RFC_",sub("\\.gz","",f),"_Model.rds") #output model file name
#     m <- aggregateMeth(feature_path=path,betas=mtx) #create aggregated meth betas
#     if (ncol(m) == 0) next
#     model <- randomForest(
#         x=m,
#         y=labels,
#         ntree=500,
#         importance=TRUE
#     )
#     saveRDS(model,paste0(out_dir,out))
#     oob <- data.frame(
#         model=f,
#         oob=model$err.rate[nrow(model$err.rate),"OOB"]
#     )
#     importance <- data.frame(
#         Model=rep(f,nrow(model$importance)),
#         Feature=rownames(model$importance),
#         MeanDecreaseAccuracy=model$importance[,"MeanDecreaseAccuracy"],
#         MeanDecreaeGini=model$importance[,"MeanDecreaseGini"]
#     ) %>% arrange(desc(MeanDecreaseAccuracy))
#     write_tsv(oob,file=oob_log,append = TRUE)
#     write_tsv(importance,file = importance_log,append = TRUE)
#     print(paste("Done with ",f))
# }

tse <- readRDS("~/Documents/HPC_share/CAP/20221010_Capper_Prospective_betas_SE.rds")

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
tmtx <- cleanMatrixForClusterW(assay(tse)) %>%
    imputeRowMean(.)

tft_mx <- dbStats(betas=tmtx,databases=feat)
# na_sum <- colSums(is.na(ft_mx))
# if (sum(na_sum) != 0) {
#     remove <- which(na_sum == nrow(ft_mx)) #sometimes no cgs from features in betas mtx, remove 
#     ft_mx <- ft_mx[,-remove]
# }
samples <- match(colnames(tmtx),colData(tse)$SampleID)
labels <- factor(colData(tse)[samples,"Capper_Methylation_Class"])

predictions <- predict(model, newdata=tft_mx)
write_tsv(as.data.frame(predictions), "~/Documents/HPC_share/chromHMM/20230814_chromhmm_predictions.tsv")
results <- data.frame(
    Sample_ID = c(samples),
    predicted = c(predictions),
    actual = c(labels)
    # probability_score = c(probabilities_ranked)
)
results$predicted66 <- master$Capper_Methylation_Class2[match(results$actual, master$Capper_Methylation_Class)]
results$actual66 <- master$Capper_Methylation_Class2[match(results$predicted, master$Capper_Methylation_Class)]

master <- read_excel("~/Documents/HPC_share/Samplesheets/20220913_CNStumor_master_samplesheet.xlsx")

test_accuracy <- sum(as.character(results$actual) == as.character(results$predicted), na.rm = TRUE) / length(labels)
print(paste("Test Accuracy:", test_accuracy))
# "Test Accuracy: 0.365367180417044"
write_tsv(as.data.frame(results), "~/Documents/HPC_share/chromHMM/20230814_chromhmm_results.tsv")

# master <- read_excel("~/Documents/HPC_share/Samplesheets/20230130_CNSTumorClassLabel.xlsx")



