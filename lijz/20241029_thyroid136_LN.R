


## INCLUDING LYMPH NODE METASTASIS

pdf('~/Documents/HPC_share/thyroid/figures/20240916_tsne136_LN.pdf', width=6, height=5, onefile=FALSE)
p <- ggplot() +
    geom_point(data = ss, aes(x = tSNE1, y = tSNE2, color = LYMPH_NODE)) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "#ffffff", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgray")
    )
plot(p)
dev.off()


## LN DIFF METH

res = readRDS(file = "~/Documents/HPC_share/thyroid/diff_meth/20240917_thyroid136_res_LN.rds")

(length(res$Probe_ID[res$Est_LYMPH_NODET > 0.2]) + length(res$Probe_ID[res$Est_LYMPH_NODET < -0.2]))
# [1] 693
length(res$Probe_ID[res$Est_LYMPH_NODET < -0.2])
# 1 v 2
res1 = testEnrichment(res$Probe_ID[res$Est_LYMPH_NODET > 0.2], platform="HM450", universe=res$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240917_thyroid136_LN_enrichment_hypo.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res1, n_label=30)
dev.off()

res1 = testEnrichment(res$Probe_ID[res$Est_LYMPH_NODET < -0.2], platform="HM450", universe=res$Probe_ID)
pdf('~/Documents/HPC_share/thyroid/figures/20240917_thyroid136_LN_enrichment_hyper.pdf', family="ArialMT", width=10, height=5, onefile=FALSE)
KYCG_plotEnrichAll(res1, n_label=30)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20240917_thyroid136_LN_TFBS_enrichment_hypo.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res$Probe_ID[res$Est_LYMPH_NODET > 0.2], "TFBSconsensus", platform="EPIC", universe=res1$Probe_ID), n_min = 40)
dev.off()

pdf('~/Documents/HPC_share/thyroid/figures/20240917_thyroid136_LN_TFBS_enrichment_hyper.pdf', family="ArialMT", width=5, height=7, onefile=FALSE)
plotDot(testEnrichment(res$Probe_ID[res$Est_LYMPH_NODET < -0.2], "TFBSconsensus", platform="EPIC", universe=res1$Probe_ID), n_min = 40)
dev.off()



