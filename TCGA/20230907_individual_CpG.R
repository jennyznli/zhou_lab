# find top classes out of all the classifiers (list of top importance ones)
# find all the cpgs associated with that class
# put them all together

importance <- read_tsv(file = "~/Documents/HPC_share/aggregated/20230817_importance.tsv", col_names=FALSE)

colnames(importance) <- c("Database", "Feature", "Accuracy", "Importance")
importance <- importance[-(298:1485),]
write_tsv(importance, "~/Documents/HPC_share/aggregated/20230817_importance.tsv")

importance$Importance <- importance$Importance[order(importance$Importance, decreasing = TRUE)]

# match filename with feature knowledgeset, search for all cpgs within that file
# bind to a file
# list should include knowledgebase set, feature, and individual cpg

#______________________________________________________________________________

x<-c("tidyr", "dplyr", "readr", "stringr")
lapply(x, require, character.only = TRUE)

importance <- read_tsv(file = "/home/lijz/HPC_share/aggregated/20230817_importance.tsv", col_names=FALSE)
dir <- "/home/lijz/20230817_aggregated/set/new/" #feature dir
features <- list.files(dir)
out_dir <- "/home/lijz/20230908_aggregated_ind/" #output dir
probe_log <- paste0(out_dir,"20230908_probe.tsv") #log for top individual probes

limit = 100
subset = importance$Database[1:limit]
databases = unique(subset)

extract1 <- function(string) {
    str_extract(string, '\\b\\w+$')
}

# read in all necessary databases
for (x in databases) {
    assign(x, read_tsv("~/Documents/HPC_share/aggregated/20230817_importance.tsv"))
    x$Knowledgebase = lapply(x$Knowledgebase, extract1)
}

for (n in 1:nrows(importance)) {
    database = importance$Database[n]
    feature = importance$Feature[n]
    focus_db = get(database)

    probes = focus_db$Probe_ID[feature == database$Knowledgebase]
    if (!(probes %in% probe_log)) {
        # IF THE PROBE IS ALREADY IN THE DATASET THEN DONT ADD IT 
    }
    probe_entry = data.frame(
        Database = database, 
        Feature = feature, 
        Probe = probes
    )
    write_tsv(probe_entry,file=probe_log,append = TRUE)
    print(paste("Done with ",n))
}

#______________________________________________________________________________

# ABCompartment.20220911.gz$Knowledgebase = lapply(ABCompartment.20220911.gz$Knowledgebase, extract1)

    database = importance$Database[3]
    feature = importance$Feature[3]
    # focus_db = get(database)                                                         e)
    probes = focus_db$Probe_ID[feature == focus_db$Knowledgebase]
    probe_entry = data.frame(
        Database = database, 
        Feature = feature, 
        Probe = probes
    )
    write_tsv(probe_entry,file=probe_log,append = TRUE)
print(paste("Done with ",n))

# 9/12______________________________________________________________________________
    
# se <- readRDS("~/Documents/HPC_share/CAP/20221210_Capper_Reference_betas_SE.rds")
se <- readRDS("/mnt/isilon/zhou_lab/projects/20220816_BrainTumorClassifier/20221210_Capper_Reference_betas_SE.rds")
features <- read_tsv("/home/lijz/20230908_aggregated_ind/20230817_importance.tsv")
# features <- read_tsv("~/Documents/HPC_share/aggregated/20230817_importance.tsv")
features <- pull(features, Feature)[1486:1585]
# features <- features[1486:1585,2]
out_dir <- "/home/lijz/20230908_aggregated_ind/" # output_dir
db <- read_tsv("/home/lijz/20230817_aggregated/set/new/ChromHMMfullStack.20230515.gz") #chromhmm stuff
# db <- read_tsv("~/Documents/HPC_share/aggregated/ChromHMMfullStack.20230515") #chromhmm stuff

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

samples <- match(colnames(mtx), colData(se)$SampleID)
# match samples with samplesheet?
labels <- factor(colData(se)[samples,"Capper_Methylation_Class"])
mtx = t(mtx)

extract1 <- function(string) {
    str_extract(string, '\\b\\w+$')
}
db$Features = lapply(db$Knowledgebase, extract1)

# name = paste0("", n)
# print(name)
# 
# sub = mtx[1:100,1:100]
# x = subset(db, Features == name)                                    
# mtx1 = subset(mtx, rownames(mtx) %in% x$Probe_ID)
# mtx1 = mtx[which(rownames(mtx) %in% x$Probe_ID),]
# 
# which(rownames(mtx) == "cg00291929")
# mtx[6137,]

for (n in features) {
    name = paste0("", n)
    print(name)
    x = subset(db, Features == name)
    
    mtx1 = subset(mtx, rownames(mtx) %in% x$Probe_ID)
    print(dim(x))
    imputeRowMean <- function(mtx) {
        k <- which(is.na(mtx), arr.ind=TRUE)
        mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
        mtx
    }
    mtx1 = imputeRowMean(mtx1)
    model <- randomForest(
        x=mtx1,
        y=labels,
        ntree=500,
        importance=TRUE
    )
    out <- paste0("20230912_RFC_",sub("\\.gz","",n),"_Model.rds") #output model file name
    saveRDS(model,paste0(out_dir,out))
    print(paste("Done with ",f))
    print("done")
}
# n = "1_GapArtf2"
# name = paste0("", n)
# print(name)
# x = subset(db, Features == name)
# 
# x = subset(db, Features == )

# x = subset(db, db$Features == "1_GapArtf2")
# # these make the same and they're both not right
# mtx2 = mtx1[rownames(mtx1) %in% x$Probe_ID,]
# mtx3 = mtx1[x$Probe_ID,]
# 
# # this works! but i think there are a lot of nas for some reason
# y <- rownames(mtx1)[rownames(mtx) %in% x$Probe_ID]
# mtx2 = subset(mtx1, rownames(mtx1) %in% y)
# 
# mtx2 = t(mtx2)
# model <- randomForest(
#     x=mtx4,
#     y=labels,
#     ntree=500,
#     importance=TRUE
# )

# 9/21______________________________________________________________________________

## TAKE ACCURACIES OF ALL MODELS 
out_dir <- "/home/lijz/20230908_aggregated_ind/" #out dir 
dir <- "/home/lijz/20230908_aggregated_ind/MODELS/" #model dir 
models <- list.files(dir)
importance_log <- paste0(out_dir,"20230921_importance.tsv") #log for feature importance for each model
oob_log <- paste0(out_dir,"/20230921_oob.tsv") #log for oob error rates for each model

for (f in models) {
    path <- paste0(dir,f) #model file path
    model <- readRDS(path)
    oob <- data.frame(
        model=f,
        oob=model$err.rate[nrow(model$err.rate),"OOB Error"]
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

# 9/22______________________________________________________________________________
## Have to make classifier with all features
## collect the top feature sets

library(readr)
library(ExperimentHub)

oob <- as.data.frame(read_tsv("~/Documents/HPC_share/aggregated/20230921_oob.tsv"))
colnames(oob) <- c("Model", "Error")
oob1 = oob[order(oob$Error),]



