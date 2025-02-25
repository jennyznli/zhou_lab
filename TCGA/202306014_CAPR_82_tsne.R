x<-c("sesame", "askme", "sesameData", "tidyr", "readxl", "readr", "randomForest")
lapply(x, require, character.only = TRUE)

labels <- read_tsv("/Users/jennyzli/Documents/HPC_share/CAP/CapperRef2682_CNS82.tsv")
DGDdf <- read_tsv("/Users/jennyzli/Documents/HPC_share/DGD/DGDdf.tsv")
CAPdf <- read_tsv("/Users/jennyzli/Documents/HPC_share/DGD/CAPdf.tsv")
CAPdf$Predicted = labels$label[match(CAPdf$sample, labels$Sample_ID)]
labels2 <- read_excel("/Users/jennyzli/Documents/HPC_share/DGD/20230609_DGDdf.xlsx")
DGDdf$Predicted = labels2$CNS82[match(DGDdf$DGD_ID, labels2$DGD)]


correctDF = DGDdf %>% dplyr::filter(Status=="C")
wrongDF = DGDdf %>% dplyr::filter(Status=="W")
sarcomaDF = DGDdf %>% dplyr::filter(Status=="S")
unsuitableDF = DGDdf %>% dplyr::filter(Status=="N")

write.xlsx(CAPdf, "~/Documents/HPC_share/DGD/20230609_CAPdf.xlsx")

CAPdf <- read_excel("~/Documents/HPC_share/DGD/20230609_CAPdf.xlsx")
CAPdf2 <- read_tsv("~/Documents/HPC_share/DGD/CAPdf.tsv")
CAPdf$CNS66 <- CAPdf2$Predicted

# html cohorts
ggplotly(ggplot() + 
             geom_point(data = CAPdf, aes(x = tSNE1, y = tSNE2, text = (paste("Sample:", sample, "\n", "Predicted:", Predicted))), size=.25, shape = 16, color = "gray") + 
             geom_point(data = correctDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis, "\n", "Probe:", Probes), color=Predicted), size=0.5, shape = 16) + 
             geom_point(data = wrongDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis, "\n", "Probe:", Probes), color=Predicted), size = 2, shape = 4) + 
             geom_point(data = sarcomaDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis, "\n", "Probe:", Probes), color=Predicted), size = 1.5, shape = 2) + 
             geom_point(data = unsuitableDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis, "\n", "Probe:", Probes), color=Predicted), size = 1.5, shape = 1) +
             theme_bw())



# html all labels
ggplotly(ggplot() + 
             geom_point(data = CAPdf, aes(x = tSNE1, y = tSNE2, text = (paste("Sample:", sample)), color=Predicted, alpha = .25), size=.25, shape = 16) + 
             geom_point(data = correctDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis), color=Predicted), alpha = correctDF$newScore, size=1, shape = 17) + 
             geom_point(data = wrongDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis), color=Predicted), alpha = wrongDF$newScore, size = 2, shape = 4) + 
             geom_point(data = sarcomaDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis), color=Predicted), alpha = sarcomaDF$newScore, size = 1.5, shape = 2) + 
             geom_point(data = unsuitableDF, aes(tSNE1, tSNE2, text = paste("Sample:", DGD_ID, "\n", "Diagnosis:", Diagnosis), color=Predicted), alpha = unsuitableDF$newScore, size = 1.5, shape = 1) +
             theme_bw())