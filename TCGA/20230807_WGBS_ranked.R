x<-c("tidyr", "dplyr", "randomForest", "readr", "readxl")
lapply(x, require, character.only = TRUE)

coords <- read.delim("~/Documents/HPC_share/WGBS/HM450.hg38.coords", header=FALSE)
betas <- read_tsv(file ="~/Documents/HPC_share/WGBS/wgbs_TCGA_GBM_0128_beta.tsv", col_names = FALSE)
betas1401 <- read_tsv(file ="~/Documents/HPC_share/WGBS/wgbs_TCGA_GBM_1401_beta.tsv", col_names = FALSE)

hiImp <- read_tsv("~/Documents/HPC_share/CAP/EPICv12HM450CG_CapperRefHiImp6636_230117.tsv")
# 482417 original betas, 289249 are NA, 193168 are values

betas <- as.data.frame(betas)
betas$Probe_ID <- coords$V2
betas2 <- betas[!duplicated(betas$Probe_ID), ]
betas3 <- betas2[betas2$Probe_ID %in% hiImp$Probe_ID,]
# 2403/6636 have values 

# 482413 unique cpgs, 4 duplicated
# V1         V2   V3    V4    V5
# 1 217 cg13869341 chr1 15864 15866
# 2 277 cg14008030 chr1 18826 18828

WGBScpgs <- coords$V2[!(duplicated(coords$V2))]
# duplicated cpgs: "cg07434271" "cg12419862" "cg18751682" "cg06335633"
sum(WGBScpgs %in% CAPPcpgs) # all in

compare$a <- (compare$`betas1401$X1` <= 1)
compare$b <- (compare$Sample_0128 <= 1)
sum(!is.na(compare$a == compare$b))

# 118262 cpgs are in common and not NA

