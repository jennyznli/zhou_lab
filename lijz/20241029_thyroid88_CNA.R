packages <- c("tidyr", "tidyverse", "stats", "dplyr", "plotly", "readr", "readxl", "here", "pvclust",
              "stringr", "ggplot2", "sesame", "textshape", "impute", "ComplexHeatmap", "cowplot", "circlize")
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

# COPY NUMBER ANALYSIS_______________________________________________________________
sdfs = readRDS(here("data", "20241021_thyroid136_sigdf.rds"))
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>% filter(Include_In_Analysis ==1)

sample_id <- ss$Sample_ID[match(names(segs), ss$IDAT)]

segs <- lapply(sdfs, cnSegmentation)
saveRDS(segs, here("data", "20241021_thyroid136_segs.rds"))

# 136 INDIVIDUAL PLOTS _______________________________________________________________
segs <- read.csv(here("data", "20241021_thyroid136_segs.rds"))

output_folder <- "~/Documents/HPC_share/thyroid/cnv_figures/"
for (i in seq_along(segs)) {
    pdf_filename <- paste0(output_folder, "20241021_thyroid136_cna_", sample_id[i], ".pdf")
    p <- visualizeSegments(segs[[i]])
    ggsave(filename = pdf_filename,
           plot = p,
           width = 9, height = 5,
           units = "in",
           device = "pdf")
}

## PREPROCESSING ## ____________________________________________________________________________________________
df <- read.csv(here("data", "20241030_thyroid88_cn_matrix.csv")) %>%
    filter(chrom != "chrY")
df$Invasiveness = ss$Invasiveness[match(df$ID, ss$Sample_ID)]

chroms <- paste0("chr", c(1:22, "X"))
df$chrom <- factor(df$chrom,
                   levels = paste0("chr", c(1:22, "X")),
                   ordered = FALSE)
chrom_sizes <- df %>%
    group_by(chrom) %>%
    summarize(max_pos = max(loc.end))
chrom_starts <- c(0, chrom_ends[-length(chrom_ends)])
chrom_ends <- cumsum(as.numeric(chrom_sizes$max_pos))
names(chrom_starts) <- chrom_sizes$chrom
chrom_midpoints <- chrom_starts + (as.numeric(chrom_sizes$max_pos) / 2)

df <- df %>%
    group_by(chrom) %>%
    mutate(
        abs_start = loc.start + chrom_starts[as.character(chrom)],
        abs_end = loc.end + chrom_starts[as.character(chrom)]
    )
df$abs_midpoints <- ((df$abs_start + df$abs_end) / 2)

## OVERLAID VERSION _______________________________________________________

p <- ggplot(df, aes(x=abs_starts, y=seg.mean)) +
    geom_point(aes(x=abs_midpoints, y=seg.mean, color=Invasiveness), alpha=0.2, size=0.3) +
    geom_segment(aes(x=abs_start, xend=abs_start + (loc.end-loc.start),
                     y=seg.mean, yend=seg.mean, color=Invasiveness),
                 alpha=0.2, size=0.7) +
    geom_vline(xintercept=chrom_starts[-1], alpha=0.3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1), panel.grid.major.y = element_line(color="grey90")) +
    labs(x="Chromosome Position", y="Copy Number Variation") +
    scale_x_continuous(breaks=chrom_midpoints,
                       labels=gsub("chr", "", names(chrom_starts))) +
    ylim(-2.5, 2.5) +
    scale_color_manual(values=c("High"="red", "Low"="blue"))

ggsave(here("figures", "20241203_thyroid88_cnv_summary.pdf"), p, width=22, height=10)

## FACET HEATMAPS_______________________________________________________
chrom_order <- paste0("chr", c(1:22, "X"))
df$chrom <- factor(df$chrom, levels = chrom_order, ordered = TRUE)

p <- ggplot(df) +
    geom_rect(aes(xmin = loc.start, xmax = loc.end,
                  ymin = as.numeric(factor(ID)) - 0.5,
                  ymax = as.numeric(factor(ID)) + 0.5,
                  fill = seg.mean)) +
    facet_wrap(~chrom,
               scales = "free_x",
               ncol = 4,     # This will create a 4x6 grid (4 columns)
               dir = "h") +  # Fill horizontally
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, limits = c(-1, 3),
        name = "Log2 Ratio"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),  # Make facet labels bold
        panel.spacing = unit(0.5, "lines")  # Adjust space between facets
    ) +
    labs(x = "Position", y = "Sample")

ggsave(here("figures" ,"20241105_thyroid88_faceted_cnv.pdf"), p, width=15, height=15)

## HEATMAP CLUSTERED CHROMOSOME MEANS _______________________________________________________
chrom_means <- df %>%
    group_by(ID, chrom) %>%
    summarize(mean_cn = mean(seg.mean, na.rm = TRUE), .groups = 'drop') %>%
    # Create a matrix suitable for clustering
    pivot_wider(names_from = chrom,
                values_from = mean_cn,
                id_cols = ID) %>%
    column_to_rownames("ID")

# 2. Perform hierarchical clustering on chromosome means
dist_matrix <- dist(chrom_means, method = "euclidean")
hc <- hclust(dist_matrix, method = "complete")  # or try "ward.D2"

# 3. Get the ordered sample names from clustering
ordered_samples <- hc$labels[hc$order]

# 4. Update your original plot with the new sample ordering
# Add the clustering order to your original data
df$ID <- factor(df$ID, levels = ordered_samples)

## HEATMAP CLUSTERED CHROMOSOME MEANS BUT SEGMENTS_______________________________________________________
p <- ggplot(df) +
    geom_segment(aes(x = abs_start, xend = abs_end,
                     y = ID, yend = ID,
                     color = seg.mean),
                 size = 3) +
    scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Log2 Ratio",
        na.value = "gray80"
    ) +
    geom_vline(xintercept = chrom_ends,
               color = "gray80",
               linetype = "dashed") +
    scale_x_continuous(
        breaks = chrom_starts + (chrom_sizes$max_pos/2),
        labels = gsub("chr", "", levels(cnv_data$chrom)),
        name = "Chromosome"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
    ) +
    labs(y = "Sample")

# Save the main plot
ggsave(here("figures", "20241105_thyroid88_cnv_clustered_by_chroms.pdf"),
       p, width=15, height=13)

## HEATMAP CHROMOSOME MEANS CLUSTERED annotated _______________________________________________________
ha_left = rowAnnotation(
    Cluster_Group = ss$Cluster_Group,
    Invasiveness = ss$Invasiveness,
    Driver_Group = ss$Driver_Group,
    col = list(
        Cluster_Group = c(GROUP1 = "#EAB464", GROUP2 = "#62B6CB", GROUP3 = "#E86A92"),
        Invasiveness = c(High = "darkblue", Low = "white", Unknown = "grey"),
        Driver_Group = c(
            "BRAF V600E" = "#F8766D",
            "Kinase Fusion" = "#7CAE00",
            "DICER1" = "#00BFC4",
            "Ras-like" = "#C77CFF",
            "Indeterminate" = "#a3a7ab"
        )
    ),
    show_legend = TRUE,
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
        Cluster_Group = list(title = "Cluster Group", at = c("GROUP1", "GROUP2", "GROUP3")),
        Invasiveness = list(title = "Invasiveness Level"),
        Driver_Group = list(title = "Driver Group")
    )
)

column_ha = HeatmapAnnotation(
    Chromosome_Type = factor(ifelse(colnames(mat_means) %in% c("chrX"), "Sex", "Autosomal")),
    col = list(Chromosome_Type = c("Sex" = "purple", "Autosomal" = "gray")),
    show_legend = TRUE
)

# UNSPLIT VERSION
p <- Heatmap(mat_means,
             name = "Mean Log2 Ratio",
             col = colorRamp2(c(-1, 0, 3), c("blue", "white", "red")),
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             show_row_dend = TRUE,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             left_annotation = ha_left,

             # Additional customization
             row_title = "Samples",
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8)
             ),
             row_dend_width = unit(40, "mm"),
             # SPLIT THE THING
             # split = ss$Cluster_Group,
             # row_gap = unit(5, "mm"),  # Add gap between clusters

             top_annotation = column_ha,
)

pdf(here("figures", "20241105_thyroid88_cnv_means_clustered_annotated.pdf"),
    width = 20,
    height = 15,
    useDingbats = FALSE)  # Better compatibility
draw(p, padding = unit(c(20, 5, 20, 5), "mm"))  # Add padding around the plot
dev.off()

# SPLIT VERSION
p <- Heatmap(mat_means,
             name = "Mean Log2 Ratio",
             col = colorRamp2(c(-1, 0, 3), c("blue", "white", "red")),
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             show_row_dend = TRUE,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             left_annotation = ha_left,

             # Additional customization
             row_title = "Samples",
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8)
             ),
             row_dend_width = unit(40, "mm"),

             split = ss$Invasiveness,
             row_gap = unit(5, "mm"),  # Add gap between clusters

             top_annotation = column_ha,
)

pdf(here("figures", "20241105_thyroid88_cnv_means_clustered_annotated_split.pdf"),
    width = 20,
    height = 15,
    useDingbats = FALSE)  # Better compatibility
draw(p, padding = unit(c(20, 5, 20, 5), "mm"))  # Add padding around the plot
dev.off()


## TRY TO COMBINE SEGMENTS AND STUFF____________________________
# but lowkey doesn't work well ...

# library(ggh4x)
#
# # First ensure your data has the annotations
# # Merge annotations with your data
# df <- df %>%
#     left_join(ss %>%
#                   rownames_to_column("ID"),
#               by = "ID")
#
# # Define color scales for annotations
# annotation_colors <- list(
#     Cluster_Group = c(GROUP1 = "red", GROUP2 = "blue", GROUP3 = "green"),
#     Histology = c(
#         PTC = "red",
#         FTC = "green",
#         "Tall Cell" = "blue",
#         Unknown = "grey",
#         PTmC = "pink",
#         metPTC = "black",
#         Hyalinizing_trabecular = "turquoise",
#         Spindle_Epithelial = "purple",
#         Benign_Adenoma = "yellow"
#     ),
#     Risk = c(High = "darkred", Low = "pink", Unknown = "grey"),
#     Driver_Group = c(
#         "BRAF V600E" = "darkred",
#         "Kinase Fusion" = "red",
#         "DICER1" = "orange",
#         "Ras-like" = "yellow",
#         "Indeterminate" = "grey"
#     )
# )
#
# # Create the plot with annotations
# p <- ggplot(df) +
#     # Add annotation bars on the left
#     geom_tile(aes(x = -1e9, y = ID, fill = Cluster_Group), width = 5e8) +
#     geom_tile(aes(x = -2e9, y = ID, fill = Histology), width = 5e8) +
#     geom_tile(aes(x = -3e9, y = ID, fill = Invasiveness), width = 5e8) +
#     geom_tile(aes(x = -4e9, y = ID, fill = Driver_Group), width = 5e8) +
#     # Add the original CNV segments
#     geom_segment(aes(x = abs_start, xend = abs_end,
#                      y = ID, yend = ID,
#                      color = seg.mean),
#                  size = 3) +
#     # Add color scales for annotations
#     scale_fill_manual(name = "Annotations",
#                       values = unlist(annotation_colors),
#                       guide = guide_legend(override.aes = list(size = 5))) +
#     # Original color scale for CNV data
#     scale_color_gradient2(
#         low = "blue",
#         mid = "white",
#         high = "red",
#         midpoint = 0,
#         limits = c(-1, 1),
#         name = "Log2 Ratio",
#         na.value = "gray80"
#     ) +
#     geom_vline(xintercept = chrom_ends,
#                color = "gray80",
#                linetype = "dashed") +
#     scale_x_continuous(
#         breaks = c(-4e9, -3e9, -2e9, -1e9,
#                    chrom_starts + (chrom_sizes$max_pos/2)),
#         labels = c("Driver", "Risk", "Histology", "Cluster",
#                    gsub("chr", "", levels(cnv_data$chrom))),
#         name = "Chromosome"
#     ) +
#     theme_minimal() +
#     theme(
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         legend.position = "right",
#         legend.box = "vertical",
#         legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
#         legend.box.margin = margin(t = 10, r = 10, b = 10, l = 10)
#     ) +
#     labs(y = "Sample")
#
# # Save the plot
# ggsave(here("figures", "20241105_thyroid88_cnv_clustered_annotated_idek.pdf"),
#        p, width=20, height=13)  # Increased width to accommodate annotations
#
# # Optional: Create a separate legend for the annotations
# legend_plot <- ggplot(df) +
#     geom_tile(aes(fill = Cluster_Group)) +
#     geom_tile(aes(fill = Histology)) +
#     geom_tile(aes(fill = Invasiveness)) +
#     geom_tile(aes(fill = Driver_Group)) +
#     scale_fill_manual(values = unlist(annotation_colors)) +
#     theme_minimal() +
#     theme(legend.position = "right")
#
# ggsave(here("figures", "20241105_thyroid88_cnv_legend.pdf"),
#        legend_plot, width=8, height=10)

## IGV DATA PREPARATION _______________________________________________________
segs <- readRDS(here("data", "20241021_thyroid136_segs.rds"))
seg <- segs[[1]]$seg.signals[,1:6]
write.table(seg,
            file = here("segments", "20241202_thyroid136_1.seg"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
ss = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
# ss = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>% filter(Include_In_Analysis ==1)

for(i in seq_along(segs)) {
    seg <- segs[[i]]$seg.signals[,1:6]

    seg[,1] <- ss$Sample_ID[match(names(segs)[i], ss$IDAT)]

    filename <- sprintf("20241202_thyroid136_%d.seg", i)

    write.table(seg,
                file = here("segments", filename),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)

    message(sprintf("Saved file %d of %d: %s", i, length(segs), filename))
}

ss_annotated = ss %>% select(Driver_Group, Invasiveness, Lymph_Node)

write.table(ss_annotated,
            here("segments", "20241203_thyroid136_annotations.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



