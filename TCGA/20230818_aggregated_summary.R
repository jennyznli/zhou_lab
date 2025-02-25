x<-c("tidyr", "dplyr", "readr", "readxl", "ggplot2", "sesame")
lapply(x, require, character.only = TRUE)

#______________________________________________________________________________

importance <- read_tsv(file = "~/Documents/HPC_share/aggregated/20230817_importance.tsv", col_names=FALSE)
error <- read_tsv(file = "~/Documents/HPC_share/aggregated/20230817_oob.tsv", col_names=FALSE)

error <- error[-(16:25),]
colnames(error) <- c("Database", "Error")
error$Accuracy <- 1 - error$Error
write_tsv(error, "~/Documents/HPC_share/aggregated/20230817_error.tsv")

colnames(importance) <- c("Database", "Feature", "Accuracy", "Importance")
importance <- importance[-(298:1485),]
write_tsv(error, "~/Documents/HPC_share/aggregated/20230817_importance2.tsv")

importance$Importance <- importance$Importance[order(importance$Importance, decreasing = TRUE)]
error$Accuracy <- error$Accuracy[order(error$Accuracy, decreasing = TRUE)]

# ACCURACY ______________________________________________________________________

ggplot(error, aes(x = reorder(Database, -Accuracy), y = Accuracy)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(y = "Accuracy", x = "Database") + 
    ylim(0, 1) +
    ggtitle("Accuracy of CpG Aggregation with Tissue Signature Sets") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
    
# FEATURE IMPORTANCE ______________________________________________________________________

ggplot(importance, aes(x = reorder(Feature, -Importance), y = Importance)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(y = "Importance", x = "Database") + 
    ggtitle("Most Important Features of CpG Aggregation with Tissue Signature Sets") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

# CLASS ERRORS ______________________________________________________________________

## COLLECT THEM 
dir <- "~/Documents/HPC_share/MODELS/" #feature dir 
features <- list.files(dir)

out_dir <- "~/Documents/HPC_share/aggregated/" #file output dir
class_error_log <- paste0(out_dir,"20230820_class_error_log.tsv") #log for feature importance for each model

ex <- readRDS("~/Documents/HPC_share/MODELS/20230817_RFC_ABCompartment.20220911_Model.rds")
x <- as.data.frame(t(ex$confusion[,'class.error']))

for (f in features) {
    model <- readRDS(paste0(dir,f))
    class_error <- as.data.frame(t(model$confusion[,'class.error']))
    write_tsv(class_error,file=class_error_log,append = TRUE)
    print(paste("Done with ",f))
}

##########

capr <- readRDS("~/Documents/HPC_share/CAP/20221210_Capper_Reference_betas_SE.rds")

class_errors <- as.data.frame(read_tsv(file = "~/Documents/HPC_share/aggregated/20230820_class_error_log.tsv", col_names=FALSE))
labels <- as.data.frame(t(as.data.frame(strsplit(features, "_"))[3,]))
# mastersheet <- read_excel("~/Documents/HPC_share/Samplesheets/20230130_CNSTumorClassLabel.xlsx")
# length(unique(factor(mastersheet$ClassName2018))) # 91-8 = 82
# length(unique(factor(colData(capr)[,"Capper_Methylation_Class"]))) #74
# unique(factor(colData(capr)[,"WHO_2016_Pathological_Diagnosis"])) #88
# length(unique(factor(colData(capr)[,"Capper_Methylation_Subclass"]))) #25


x <- as.data.frame(t(ex$confusion[,'class.error']))
colnames(class_errors) <- colnames(x)
rownames(class_errors) <- labels$'3'
write_tsv(as.data.frame(class_errors), "~/Documents/HPC_share/aggregated/20230820_aggregated_class_errors.tsv")

class_errors = read_tsv("~/Documents/HPC_share/aggregated/20230820_aggregated_class_errors.tsv")

class_errors = t(class_errors)
col <- colorRampPalette(c("blue","white", "red"))(100)
heatmap(as.matrix(class_errors), Colv = NA, Rowv = NA, scale="row", 
        main = "Class Accuracies of Aggregated CpG Signatures", margins = c(13,18), col = col, 
        xlab="Model", ylab="Class Error")

########## REORDERING

n_class_errors = class_errors
# class_errors = t(class_errors)
Class_Average = apply(n_class_errors, 1, mean)
n_class_errors = as.data.frame(cbind(n_class_errors, Class_Average))
n_class_errors = n_class_errors[order(n_class_errors$Class_Average, decreasing=TRUE),]
n_class_errors <- n_class_errors[,-22]

Model_Average = apply(n_class_errors, 2, mean)
n_class_errors = rbind(n_class_errors, Model_Average)

x = as.data.frame(t(n_class_errors[75,]))
y = order(x$"75", decreasing=TRUE)
n_class_errors = n_class_errors[,y]
n_class_errors <- n_class_errors[-75,]

##FINAL HEATMAP##

class_accuracies = 1-(n_class_errors)
col <- colorRampPalette(c("blue","white", "red"))(300)
heatmap(as.matrix(class_accuracies), Colv = NA, Rowv = NA, scale = "none",
        main = "Class Accuracies of Aggregated CpG Signatures", margins = c(13,18), col = col, 
        xlab="Model", ylab="Class Accuracies")
write_tsv(class_accuracies, "~/Documents/HPC_share/aggregated/20230820_aggregated_class_accuracies.tsv")

gradientLegend(0:1, color=col, nCol=30, pos = 0, side = 4, inside=FALSE)


meta = read_excel("~/Documents/HPC_share/Samplesheets/20220913_CNStumor_master_samplesheet.xlsx")
meta = meta %>% dplyr::filter(Capper_Cohort!="NA", Capper_Methylation_Class2!="NA")
case_cnt = meta %>% dplyr::filter(!grepl("^Control", Capper_Methylation_Class2)) %>% pivot_wider(id_cols = "Capper_Methylation_Class2", names_from = "Capper_Cohort", values_from = "IDAT", values_fn = length, values_fill=0)
labels = case_cnt %>% dplyr::filter(Reference>=5) %>% with(Capper_Methylation_Class2)
length(labels) # 82 classes
meta = meta %>% dplyr::filter(Capper_Cohort=="Reference", Capper_Methylation_Class2 %in% labels)

table <- as.data.frame(table(meta$Capper_Methylation_Class2))
table1 <- as.data.frame(table(meta$Capper_Methylation_Class))

# find size of each class? and take out the small ones

## CLASS ACCURACY VS SIZE
ggplot(error, aes(x = reorder(Database, -Accuracy), y = Accuracy)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(y = "Accuracy", x = "Database") + 
    ylim(0, 1) +
    ggtitle("Accuracy of CpG Aggregation with Tissue Signature Sets") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))





