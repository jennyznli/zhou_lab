x<-c("tidyr", "dplyr", "randomForest","SummarizedExperiment", "readr", "readxl", "sesame")
lapply(x, require, character.only = TRUE)

meta = read_excel("/mnt/isilon/zhou_lab/projects/20220816_BrainTumorClassifier/20220913_CNStumor_master_samplesheet.xlsx")
# meta = read_excel("~/Documents/HPC_share/Samplesheets/20220913_CNStumor_master_samplesheet.xlsx")

meta = meta %>% dplyr::filter(Capper_Cohort!="NA", Capper_Methylation_Class2!="NA")
case_cnt = meta %>% dplyr::filter(!grepl("^Control", Capper_Methylation_Class2)) %>% pivot_wider(id_cols = "Capper_Methylation_Class2", names_from = "Capper_Cohort", values_from = "IDAT", values_fn = length, values_fill=0)
label = case_cnt %>% dplyr::filter(Reference>=5) %>% with(Capper_Methylation_Class2)
length(label) # 82 classes
meta = meta %>% dplyr::filter(Capper_Cohort=="Reference", Capper_Methylation_Class2 %in% label)

se <- readRDS("/mnt/isilon/zhou_lab/projects/20220816_BrainTumorClassifier/20221210_Capper_Reference_betas_SE.rds")
# se <- readRDS("~/Documents/HPC_share/CAP/20221210_Capper_Reference_betas_SE.rds")

se = se[,meta$GSM]

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

mtx <- assay(se)
samples <- match(colnames(mtx),colData(se)$SampleID)
labels <- as.factor(meta$Capper_Methylation_Class2[match(samples, meta$GSM)])
print(labels)

aggregateMeth <- function(feature_path,betas,...) {
    feat <- KYCG_loadDBs(feature_path)
    names(feat) <- vapply(feat,function(x) attr(x,"dbname"),character(1))
    ft_mx <- dbStats(betas=betas,databases=feat,...)
    na_sum <- colSums(is.na(ft_mx))
    if (sum(na_sum) != 0) {
        remove <- which(na_sum == nrow(ft_mx)) #sometimes no cgs from features in betas mtx, remove 
        ft_mx <- ft_mx[,-remove]
    }
    ft_mx
}
dir <- "/home/lijz/20230817_aggregated/set/new/" #feature dir 

features <- list.files(dir)
out_dir <- "~/20230817_aggregated/MODELS2/" #file output dir
importance_log <- paste0(out_dir,"20230820_importance.tsv") #log for feature importance for each model
oob_log <- paste0(out_dir,"/20230820_oob.tsv") #log for oob error rates for each model

for (f in features) {
    path <- paste0(dir,f) #feature file path
    out <- paste0("20230817_RFC_",sub("\\.gz","",f),"_Model.rds") #output model file name
    m <- aggregateMeth(feature_path=path,betas=mtx) #create aggregated meth betas
    if (ncol(m) == 0) next
    model <- randomForest(
        x=m,
        y=labels,
        ntree=500,
        importance=TRUE
    )
    saveRDS(model,paste0(out_dir,out))
    oob <- data.frame(
        model=f,
        oob=model$err.rate[nrow(model$err.rate),"OOB"]
    )
    importance <- data.frame(
        Model=rep(f,nrow(model$importance)),
        Feature=rownames(model$importance),
        MeanDecreaseAccuracy=model$importance[,"MeanDecreaseAccuracy"],
        MeanDecreaeGini=model$importance[,"MeanDecreaseGini"]
    ) %>% arrange(desc(MeanDecreaseAccuracy))
    write_tsv(oob,file=oob_log,append = TRUE)
    write_tsv(importance,file = importance_log,append = TRUE)
    print(paste("Done with ",f))
}

