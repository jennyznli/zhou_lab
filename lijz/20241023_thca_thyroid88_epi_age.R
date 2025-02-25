x<-c("tidyr", "dplyr", "plotly", "readr", "readxl","here", "pvclust", "stringr", "ggplot2", "sesame", "Rtsne", "impute", "pheatmap")
lapply(x, require, character.only = TRUE)
here()

## PREDICT EPIGENETIC CLOCK AGES _____________________________________________________________
betas = readRDS(here("data", "20240320_thyroid136_betas_condensed.rds"))
write.csv(betas, here("data", "20241113_thyroid136_betas_condensed.csv"))
p <- read.csv(here("prev_data", "datMiniAnnotation4_fixed.csv"))

## pediatric and thca
model = readRDS(here("prev_data", "AgeClock_HumanHorvathN353_HM450.rds"))
betas88 = readRDS(here("data", "20240320_thyroid136_betas_condensed.rds")) %>% mLiftOver("HM450")
betas496 = readRDS(here("data", "20240915_thca_betas.rds"))

age_predictions88 <- apply(betas88, 2, function(column) predictAge(column, model))
write.csv(age_predictions88, here("data", "20240915_thyroid136_age_predictions_horvath.csv"), row.names = TRUE)

age_predictions496 <- apply(betas496, 2, function(column) predictAge(column, model))
write.csv(age_predictions496, here("data", "20240915_thca496_age_predictions_horvath.csv"), row.names = TRUE)

## COMBINED FOLD CHANGE DOT PLOT ___________________________________________________________________________________________________________________________

## PREPROCESSING
ss88 = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>% select('Source',"Cluster_Group", 'Sample_ID', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex')
ss88$Actual_Age <- as.numeric(ss88$Actual_Age)
ss88$Predicted_Age <- as.numeric(ss88$Predicted_Age)
# ss88$Fold_Change <- ss88$Predicted_Age / ss88$Actual_Age

ss496 = read_excel(here("ss", "20241023_thca_master.xlsx")) %>%  select('Source', "Cluster_Group", 'Sample_ID', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex' )
ss496$Actual_Age <- as.numeric(ss496$Actual_Age)
ss496$Predicted_Age <- as.numeric(ss496$Predicted_Age)
# ss496$Fold_Change <- ss496$Predicted_Age / ss496$Actual_Age

ss = rbind(ss88, ss496)

## PLOTTING ##

# driver group
p <- ggplot(ss, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Driver_Group, shape = Source), size = 1) +
    scale_color_manual(values = c("BRAF V600E" = "#F8766D",
                                  "DICER1" = "#00BFC4",
                                  "Kinase Fusion" = "#7CAE00",
                                  "Ras-like" = "#C77CFF",
                                  "Indeterminate" = "gray",
                                  NA = "gray")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", size = 0.3) +
    scale_shape_manual(values = c("THCA" = 16, "PED" = 17)) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_combined632_actual_predicted_age_driver.pdf"),
       plot = p, width = 7, height = 6)




cluster_colors <- c(
    "HI" = "#f8766dff",
    "HIL" = "#f7b456ff",
    "LI" = "#00bfc4ff"
)
custom_shapes <- c("THCA" = 16,
                   "PED" = 1)

ss_clean <- ss[!(ss$Cluster_Group == "NA"), ]
p <- ggplot(ss_clean, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Cluster_Group, shape = Source), size = 1) +
    scale_color_manual(values = cluster_colors) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", size = 0.3) +
    scale_shape_manual(values = custom_shapes) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20250204_combined632_actual_predicted_age_invasiveness.pdf"),
       plot = p, width = 6, height = 5)

# sex
p <- ggplot(ss, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Sex, shape = Source), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", size = 0.3) +
    scale_shape_manual(values = c("THCA" = 16, "PED" = 17)) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_combined632_actual_predicted_age_sex.pdf"),
       plot = p, width = 7, height = 6)

### PEDIATRIC PLOT 136 ____________________________________________________________
ss136 = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>% select('Source', 'Sample_ID', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex', 'Lymph_Node' )
ss136$Actual_Age <- as.numeric(ss136$Actual_Age)
ss136$Predicted_Age <- as.numeric(ss136$Predicted_Age)

p <- ggplot(ss136, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Sex), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    # scale_shape_manual(values = c("THCA" = 16, "PED" = 17)) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid136_actual_predicted_age_sex.pdf"),
       plot = p, width = 7, height = 6)


p <- ggplot(ss136, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Driver_Group), size = 1) +
    scale_color_manual(values = c("BRAF V600E" = "#F8766D",
                                  "DICER1" = "#00BFC4",
                                  "Kinase Fusion" = "#7CAE00",
                                  "Ras-like" = "#C77CFF",
                                  "Indeterminate" = "#a3a7ab")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    # scale_shape_manual(values = c("THCA" = 16, "PED" = 17)) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid136_actual_predicted_age_driver.pdf"),
       plot = p, width = 7, height = 6)


p <- ggplot(ss136, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Invasiveness), size = 1) +
    scale_color_manual(values = c("High" = "#F8766D",
                                  "Low" = "#00BFC4",
                                  "NA" = "lightgray")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    # scale_shape_manual(values = c("THCA" = 16, "PED" = 17)) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid136_actual_predicted_age_invasiveness.pdf"),
       plot = p, width = 7, height = 6)

### PEDIATRIC PLOT 88 ____________________________________________________________
ss88 = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>%
    dplyr::filter(Include_In_Analysis == 1) %>%
    select('Source', 'Sample_ID', 'Cluster_Group', 'Histology', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex')
ss88$Actual_Age <- as.numeric(ss88$Actual_Age)
ss88$Predicted_Age <- as.numeric(ss88$Predicted_Age)

p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Cluster_Group), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid88_actual_predicted_age_cluster.pdf"),
       plot = p, width = 6, height = 5)

p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Histology), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid88_actual_predicted_age_histology.pdf"),
       plot = p, width = 6, height = 5)

p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Sex), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid88_actual_predicted_age_sex.pdf"),
       plot = p, width = 6, height = 5)

p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Driver_Group), size = 1) +
    scale_color_manual(values = c("BRAF V600E" = "#F8766D",
                                  "DICER1" = "#00BFC4",
                                  "Kinase Fusion" = "#7CAE00",
                                  "Ras-like" = "#C77CFF",
                                  "Indeterminate" = "#a3a7ab")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid88_actual_predicted_age_driver.pdf"),
       plot = p, width = 6, height = 5)

p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Invasiveness), size = 1) +
    scale_color_manual(values = c("High" = "#F8766D",
                                  "Low" = "#00BFC4",
                                  "NA" = "lightgray")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thyroid88_actual_predicted_age_invasiveness.pdf"),
       plot = p, width = 6, height = 5)

### THCA PLOT 496 ____________________________________________________________
ss496 = read_excel(here("ss", "20241023_thca_master.xlsx")) %>%
    select('Source', 'Sample_ID', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex')
ss496$Actual_Age <- as.numeric(ss496$Actual_Age)
ss496$Predicted_Age <- as.numeric(ss496$Predicted_Age)
ss496$Fold_Change <- ss496$Predicted_Age / ss496$Actual_Age

# Pearson’s correlation coefficient (r) was used to measure the strength of the linear association between epigenetic age and chronological age among survivors or noncancer controls.

p <- ggplot(ss496, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Sex), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thca496_actual_predicted_age_sex.pdf"),
       plot = p, width = 7, height = 6)

p <- ggplot(ss496, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Driver_Group), size = 1) +
    scale_color_manual(values = c("BRAF V600E" = "#F8766D",
                                  "DICER1" = "#00BFC4",
                                  "Kinase Fusion" = "#7CAE00",
                                  "Ras-like" = "#C77CFF",
                                  "Indeterminate" = "#a3a7ab")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thca496_actual_predicted_age_driver.pdf"),
       plot = p, width = 7, height = 6)

p <- ggplot(ss496, aes(x = Actual_Age, y = Predicted_Age, size = 0.15)) +
    geom_point(aes(color = Invasiveness), size = 1) +
    scale_color_manual(values = c("High" = "#F8766D",
                                  "Low" = "#00BFC4",
                                  "NA" = "lightgray")) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey", linewidth = 0.3) +
    labs(
        x = "Actual Age",
        y = "Predicted Age",
    ) +
    theme_minimal()
ggsave(here("figures", "20241114_thca496_actual_predicted_age_invasiveness.pdf"),
       plot = p, width = 7, height = 6)

### INTERACTIVE PLOT___________________________________________________________

p <- plot_ly() %>%
    # invasiveness
    add_trace(
        data = ss88,
        x = ~Actual_Age,
        y = ~Predicted_Age,
        type = 'scatter',
        mode = 'markers',
        color = ~Invasiveness,
        name = ~Invasiveness,
        legendgroup = "Invasiveness",
        visible = "legendonly",
        showlegend = TRUE,
        text = ~paste("Sample:", Sample_ID,
                      "<br>Invasiveness:", Invasiveness,
                      "<br>Age:", Actual_Age,
                      "<br>Predicted Age:", round(Predicted_Age, 1))
    ) %>%
    # Driver Groups
    add_trace(
        data = ss88,
        x = ~Actual_Age,
        y = ~Predicted_Age,
        type = 'scatter',
        mode = 'markers',
        color = ~Driver_Group,
        colors = c("BRAF V600E" = "#F8766D",
                   "DICER1" = "#00BFC4",
                   "Kinase Fusion" = "#7CAE00",
                   "Ras-like" = "#C77CFF",
                   "Indeterminate" = "#a3a7ab"),
        name = ~Driver_Group,
        legendgroup = "Driver_Group",
        visible = TRUE,
        showlegend = TRUE,
        text = ~paste("Sample:", Sample_ID,
                      "<br>Driver:", Driver_Group,
                      "<br>Age:", Actual_Age,
                      "<br>Predicted Age:", round(Predicted_Age, 1))
    ) %>%
    # Add diagonal line y=x
    add_trace(
        x = c(0, max(ss88$Actual_Age)),
        y = c(0, max(ss88$Actual_Age)),
        type = 'scatter',
        mode = 'lines',
        line = list(color = 'grey', width = 0.5),
        showlegend = FALSE,
        hoverinfo = 'none'
    ) %>%
    # Layout settings
    layout(
        title = "Age Prediction Interactive Plot",
        xaxis = list(
            title = "Actual Age",
            zeroline = FALSE,
            range = c(0, max(ss88$Actual_Age))
        ),
        yaxis = list(
            title = "Predicted Age",
            zeroline = FALSE,
            range = c(0, max(ss88$Predicted_Age))
        ),
        legend = list(
            title = list(text = "Variables"),
            itemsizing = "constant",
            itemwidth = 30
        ),
        hovermode = "closest"
    )

htmlwidgets::saveWidget(p, here("figures", "20241114_ss88_interactive_age_prediction.html"))

## REGRESSION ANALYSIS EPIGENETIC AGE ACCELERATION __________________________________

## PEDIATRIC COHORT
ss88 = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
ss88$Actual_Age <- as.numeric(ss88$Actual_Age)
ss88$Predicted_Age <- as.numeric(ss88$Predicted_Age)

model <- lm(Predicted_Age ~ Actual_Age, data = ss88)
r2 <- round(summary(model)$r.squared, 3)
pval <- format.pval(summary(model)$coefficients[2,4], digits = 3)

# not specific
p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    annotate("text", x = min(ss88$Actual_Age), y = max(ss88$Predicted_Age),
             label = paste0("R² = ", r2, "\np = ", pval),
             hjust = 0, vjust = 1) +
    theme_minimal() +
    labs(x = "Chronological Age", y = "Predicted Epigenetic Age")
ggsave(here("figures", "20241114_thyroid88_age_acceleration_regression.pdf"),
       plot = p, width = 7, height = 6)

# invasiveness
p <- ggplot(ss88, aes(x = Actual_Age, y = Predicted_Age, color = Invasiveness)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +  # Add reference line
    scale_color_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4")) +
    theme_minimal() +
    labs(x = "Chronological Age", y = "Epigenetic Age (years)") +
    stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~")))
ggsave(here("figures", "20241114_thyroid88_age_acceleration_regression_invasiveness.pdf"),
       plot = p, width = 7, height = 6)


## THCA COHORT ##
ss496 = read_excel(here("ss", "20241023_thca_master.xlsx"))
ss496$Actual_Age <- as.numeric(ss496$Actual_Age)
ss496$Predicted_Age <- as.numeric(ss496$Predicted_Age)

ss496 = read_excel(here("ss", "20241023_thca_master.xlsx")) %>% filter(Invasiveness == "High")
ss496 = read_excel(here("ss", "20241023_thca_master.xlsx"))%>% filter(Invasiveness == "Low")
ss496 = read_excel(here("ss", "20241023_thca_master.xlsx")) %>%
    filter(Invasiveness %in% c("High", "Low"))

# make model
model <- lm(Predicted_Age ~ Actual_Age, data = ss496)
r2 <- round(summary(model)$r.squared, 3)
pval <- format.pval(summary(model)$coefficients[2,4], digits = 3)

# print out stats
intercept <- coef(model)[1]
slope <- coef(model)[2]
r2 <- summary(model)$r.squared
equation <- sprintf("y = %.2f + %.2fx, R² = %.3f", intercept, slope, r2)
print(equation)

# high - "y = 26.40 + 0.70x, R² = 0.450"
# low - "y = 25.25 + 0.71x, R² = 0.377"

# general
p <- ggplot(ss496, aes(x = Actual_Age, y = Predicted_Age)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    annotate("text", x = min(ss496$Actual_Age), y = max(ss496$Predicted_Age),
             label = paste0("R² = ", r2, "\np = ", pval),
             hjust = 0, vjust = 1) +
    theme_minimal() +
    labs(x = "Chronological Age", y = "Predicted Epigenetic Age")
ggsave(here("figures", "20241114_thca496_age_acceleration_regression.pdf"),
       plot = p, width = 7, height = 6)

# INVASIVENESS
p <- ggplot(ss496, aes(x = Actual_Age, y = Predicted_Age, color = Invasiveness)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +  # Reference line
    scale_color_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4")) +
    theme_minimal() +
    labs(x = "Chronological Age", y = "Epigenetic Age (years)") +
    stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~")))
ggsave(here("figures", "20241114_thca496_age_acceleration_regression_invasiveness.pdf"),
       plot = p, width = 8, height = 7)

## AGE ACCELERATION BOXPLOTS __________________________________

# PEDIATRIC 88
ss88 = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
ss88$Actual_Age <- as.numeric(ss88$Actual_Age)
ss88$Predicted_Age <- as.numeric(ss88$Predicted_Age)
ss88 <- ss88 %>%
    filter(Include_In_Analysis == 1) %>%
    filter(Driver_Group %in% c("BRAF V600E", "DICER1", "Kinase Fusion", "Ras-like")) %>%
    mutate(
        age_acceleration = Predicted_Age - Actual_Age,
        age_acceleration_residual = residuals(lm(Predicted_Age ~ Actual_Age))
    )
# by driver group and invasiveness
p <- ggplot(ss88, aes(x = Driver_Group, y = age_acceleration, fill = Invasiveness)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_point(position = position_dodge(width = 0.8), alpha = 0.5) +
    scale_fill_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4")) +
    theme_minimal() +
    labs(x = "Driver Group",
         y = "Age Acceleration (Years)") +
    theme(axis.text.x = element_text(hjust = 1))
ggsave(here("figures", "20241114_thyroid88_age_acceleration_invasiveness_driver.pdf"),
       plot = p, width = 7, height = 6)

# TOTAL PEDIATRIC TO TEST BENIGN...
ss136 = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
ss136$Actual_Age <- as.numeric(ss136$Actual_Age)
ss136$Predicted_Age <- as.numeric(ss136$Predicted_Age)
ss136 <- ss136 %>%
    mutate(
        age_acceleration = Predicted_Age - Actual_Age,
        age_acceleration_residual = residuals(lm(Predicted_Age ~ Actual_Age))
    )

p <- ggplot(ss136, aes(x = Histology, y = age_acceleration_residual)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(y = "Age Acceleration (Years)")
ggsave(here("figures", "20241114_thyroid136_age_acceleration_histology.pdf"),
       plot = p, width = 7, height = 6)


## THCA COHORT
ss496 = read_excel(here("ss", "20241023_thca_master.xlsx")) %>%
    filter(Invasiveness %in% c("Low", "High"))

ss496$Actual_Age <- as.numeric(ss496$Actual_Age)
ss496$Predicted_Age <- as.numeric(ss496$Predicted_Age)
ss496 <- ss496 %>%
    mutate(
        age_acceleration = Predicted_Age - Actual_Age,
        age_acceleration_residual = residuals(lm(Predicted_Age ~ Actual_Age))
    )

p <- ggplot(ss496, aes(x = Driver_Group, y = age_acceleration, fill = Invasiveness)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_point(position = position_dodge(width = 0.8), alpha = 0.5) +
    scale_fill_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4")) +
    theme_minimal() +
    labs(x = "Driver Group",
         y = "Age Acceleration (Years)") +
    theme(axis.text.x = element_text(hjust = 1))
ggsave(here("figures", "20241114_thca496_age_acceleration_invasiveness_driver.pdf"),
       plot = p, width = 7, height = 6)

## BINNED ANALYSIS ___________

# Create age bins
ss88 = read_excel(here("ss", "20231102_thyroid_master.xlsx"))
ss88$Actual_Age <- as.numeric(ss88$Actual_Age)
ss88$Predicted_Age <- as.numeric(ss88$Predicted_Age)
ss88 <- ss88 %>%
    filter(Include_In_Analysis == 1) %>%
    mutate(
        age_acceleration = Predicted_Age - Actual_Age,
        age_acceleration_residual = residuals(lm(Predicted_Age ~ Actual_Age))
    )

ss496 = read_excel(here("ss", "20241023_thca_master.xlsx"))
ss496$Actual_Age <- as.numeric(ss496$Actual_Age)
ss496$Predicted_Age <- as.numeric(ss496$Predicted_Age)
ss496 <- ss496 %>%
    mutate(
        age_acceleration = Predicted_Age - Actual_Age,
        age_acceleration_residual = residuals(lm(Predicted_Age ~ Actual_Age))
    )

ss88 <- ss88 %>%
    mutate(age_bin = case_when(
        Actual_Age < 10 ~ "<10",
        Actual_Age >= 10 & Actual_Age <= 14 ~ "10-14",
        Actual_Age > 14 ~ ">14"
    ))
ss496 <- ss496 %>%
    mutate(age_bin = ifelse(Actual_Age < 45, "<45", "≥45"))

# Create boxplots
p1 <- ggplot(ss88, aes(x = age_bin, y = age_acceleration)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(x = "Age Group", y = "Age Acceleration (Years)",
         title = "Pediatric Cohort") +
    scale_x_discrete(limits = c("<10", "10-14", ">14"))
ggsave(here("figures", "20241114_thyroid88_age_acceleration_binned.pdf"),
       plot = p1, width = 7, height = 6)

p2 <- ggplot(ss496, aes(x = age_bin, y = age_acceleration)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.3) +
    theme_minimal() +
    labs(x = "Age Group", y = "Age Acceleration (Years)",
         title = "Adult Cohort") +
    scale_x_discrete(limits = c("<45", "≥45"))
ggsave(here("figures", "20241114_thca496_age_acceleration_binned.pdf"),
       plot = p2, width = 7, height = 6)

### COMPARISON WITH ACTUAL AGE ____________________________________________________________
# by invasiveness
median_ages <- ss88 %>%
    group_by(ss88$Invasiveness) %>%
    summarize(median_age = median(as.numeric(ss88$`Age at Surgery`), na.rm = TRUE))
print(median_ages)

median_ages <- ss496 %>%
    group_by(ss496$Invasiveness) %>%
    summarize(median_age = median(as.numeric(ss496$Diagnosis_Age), na.rm = TRUE))
print(median_ages)

median_ages <- ss %>%
    group_by(Invasiveness) %>%
    summarize(median_age = median(Predicted_Age, na.rm = TRUE))

# Print the means to check
print(median_ages)
# Invasiveness median_age
# <chr>             <dbl>
# 1 High               63.2
# 2 Low                64.7
# 3 NA                 64.8

mean_age = mean(ss496_age, na.rm = TRUE)
sd_age = sd(ss496_age, na.rm = TRUE)
ss496$Predicted_Age <- as.factor(ss496$Predicted_Age)

x <- ss496 %>%
    select(Predicted_Age) %>%
    mutate(Group = "Adult")
y <- ss88 %>%
    select(Predicted_Age) %>%
    mutate(Group = "Pediatric")
z = rbind(x, y)
colnames(z) = c( "Predicted_Age", "Group")

pdf('~/Documents/HPC_share/thyroid/figures/20240926_thyroid88_age_boxplot.pdf', width=3, height=5, onefile=FALSE)
p <- ggplot(z, aes(x = Group, y = Predicted_Age, fill = Group)) +
    geom_boxplot() +
    # geom_jitter(width = 0.2, alpha = 0.1) +  # Add jittered points for individual samples
    labs(x = "Cohort",
         y = "Mean Epigenetic Age (Years)") +
    theme_minimal() +
    guides(fill = "none")
plot(p)
dev.off()

# could cluster to high and low invasive?

## BINNED PEDIATRIC 136____________________________________________________________________________________________________________________________

ss = read_excel("~/Documents/HPC_share/thyroid_trash/20231102_thyroid_master.xlsx")

ss_age <- ss %>% select('IDAT', 'Actual_Age', 'Predicted_Age')
ss_age$Actual_Age <- as.numeric(ss_age$Actual_Age)
ss_age$Predicted_Age <- as.numeric(ss_age$Predicted_Age)

bin1 <- ss_age %>% filter(Actual_Age < 10)
bin2 <- ss_age %>% filter(Actual_Age >= 10 & Actual_Age <= 14)
bin3 <- ss_age %>% filter(Actual_Age > 14)
#
# median(bin1$Actual_Age) #9.552778
# median(bin2$Actual_Age) #12.37222
# median(bin3$Actual_Age) #16.21528
#
# Calculate medians for actual and predicted ages in each bin
median_bin1_actual <- median(bin1$Actual_Age, na.rm = TRUE)
median_bin2_actual <- median(bin2$Actual_Age, na.rm = TRUE)
median_bin3_actual <- median(bin3$Actual_Age, na.rm = TRUE)

median_bin1_predicted <- median(bin1$Predicted_Age, na.rm = TRUE)
median_bin2_predicted <- median(bin2$Predicted_Age, na.rm = TRUE)
median_bin3_predicted <- median(bin3$Predicted_Age, na.rm = TRUE)

# Create a combined dataset for the barplot
median_data <- data.frame(
    Bin = rep(c('< 10', '10 - 14', '> 14'), each = 2),
    Type = rep(c("Actual Age", "Predicted Age"), times = 3),
    Median_Age = c(median_bin1_actual, median_bin1_predicted,
                   median_bin2_actual, median_bin2_predicted,
                   median_bin3_actual, median_bin3_predicted)
)

pdf('~/Documents/HPC_share/thyroid/figures/20241023_thyroid137_actual_predicted_age_binned.pdf', width=7, height=6, onefile=FALSE)
p <- ggplot(median_data, aes(x = fct_inorder(Bin), y = Median_Age, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +

    # Add labels and title
    labs(
        x = "Age Bin",
        y = "Median Age") +

    # Customize the colors and theme
    scale_fill_manual(values = c("Actual Age" = "blue", "Predicted Age" = "red")) +
    theme_minimal() +
    theme(legend.title = element_blank())  # Remove legend title
plot(p)
dev.off()


## BINNED VERSION THCA 496__________________________________________________________________

ss = read_excel("~/Documents/HPC_share/thyroid/thca_tcga_clinical_data_JRF_07-08-2024-1.xlsx")

ss_age <- ss %>% select(Sample_ID, Diagnosis_Age, Predicted_Age)
ss_age$Predicted_Age <- as.numeric(ss_age$Predicted_Age)
ss_age$Diagnosis_Age <- as.numeric(ss_age$Diagnosis_Age)
ss_age$Fold_Change <- ss_age$Predicted_Age / ss_age$Diagnosis_Age

bin1 <- ss_age %>% filter(Diagnosis_Age < 45)
bin2 <- ss_age %>% filter(Diagnosis_Age >= 45)

# Calculate medians for actual and predicted ages in each bin
median_bin1_actual <- median(bin1$Diagnosis_Age, na.rm = TRUE)
median_bin2_actual <- median(bin2$Diagnosis_Age, na.rm = TRUE)

median_bin1_predicted <- median(bin1$Predicted_Age, na.rm = TRUE)
median_bin2_predicted <- median(bin2$Predicted_Age, na.rm = TRUE)

# Create a combined dataset for the barplot
median_data <- data.frame(
    Bin = rep(c('< 45', '≥ 45'), each = 2),
    Type = rep(c("Actual Age", "Predicted Age"), times = 2),
    Median_Age = c(median_bin1_actual, median_bin1_predicted,
                   median_bin2_actual, median_bin2_predicted)
)

pdf('~/Documents/HPC_share/thyroid/figures/20241023_thca496_actual_predicted_age_binned.pdf', width=5, height=6, onefile=FALSE)
p <- ggplot(median_data, aes(x = Bin, y = Median_Age, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +

    # Add labels and title
    labs(
        x = "Age Bin",
        y = "Median Age") +

    # Customize the colors and theme
    scale_fill_manual(values = c("Actual Age" = "blue", "Predicted Age" = "red")) +
    theme_minimal() +
    theme(legend.title = element_blank())  # Remove legend title
plot(p)
dev.off()

## FOLD CHANGE PLOT ADULT___________________________

pdf('~/Documents/HPC_share/thyroid/figures/20241023_thca496_actual_predicted_age.pdf', width=7, height=6, onefile=FALSE)
p <- ggplot(ss_age, aes(x = Diagnosis_Age, y = Predicted_Age)) +
    # Scatter plot with size and color based on fold change
    geom_point(aes(size = Fold_Change, color = Fold_Change), alpha = 0.7) +

    # Add a line for perfect prediction (Actual = Predicted)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +

    # Customize size scale for better visibility
    # scale_size_continuous(range = c(2, 8)) +

    # Color gradient to highlight fold change
    scale_color_gradient(low = "blue", high = "red") +

    # Titles and labels
    labs(
        x = "Actual Age",
        y = "Predicted Age",
        size = "Fold Change",
        color = "Fold Change") +

    # Theme customization for cleaner look
    theme_minimal()
plot(p)
dev.off()


## BOX PLOT BY INVASIVENESS AND AGE COHORT____________________________________________________________________
ss88 = read_excel(here("ss", "20231102_thyroid_master.xlsx")) %>% filter(Include_In_Analysis == 1) %>% select('Source', 'Sample_ID', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex')
ss88$Actual_Age <- as.numeric(ss88$Actual_Age)
ss88$Predicted_Age <- as.numeric(ss88$Predicted_Age)

ss496 = read_excel(here("ss", "20241023_thca_master.xlsx")) %>%
    select('Source', 'Sample_ID', 'Actual_Age', 'Predicted_Age', 'Driver_Group', 'Invasiveness', 'Sex') %>%
    filter(!(Invasiveness == "NA"))
ss496$Actual_Age <- as.numeric(ss496$Actual_Age)
ss496$Predicted_Age <- as.numeric(ss496$Predicted_Age)

ss = rbind(ss88, ss496)

ss <- ss %>%
    mutate(
        Age_Acceleration = Predicted_Age - Actual_Age,
        Age_Acceleration_residual = residuals(lm(Predicted_Age ~ Actual_Age))
    )

p <- ggplot(ss, aes(x = Source, y = age_acceleration, fill = Invasiveness)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_point(position = position_dodge(width = 0.8), alpha = 0.2) +
    scale_fill_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4")) +
    theme_minimal() +
    labs(x = "Cohort",
         y = "Age Acceleration (Years)") +
    theme(axis.text.x = element_text(hjust = 1))
ggsave(here("figures", "20250204_thyroid88_age_acceleration_invasiveness_cohort.pdf"),
       plot = p, width=4.5, height=4.5)



## BOX PLOT BY INVASIVENESS AND AGE COHORT WITH WILCOX____________________________________________________________________

format_pval <- function(p) {
    if(p < 0.001) return("p < 0.001")
    return(sprintf("p = %.3f", p))
}

cluster_colors <- c(
    "HI" = "#f8766dff",
    "HIL" = "#f7b456ff",
    "LI" = "#00bfc4ff"
)

invasiveness_colors <- c(
    "High" = "#f8766dff",
    "Low" = "#00bfc4ff"
)

stats_list <- lapply(unique(ss$Source), function(src) {
    data_subset <- ss[ss$Source == src, ]
    test <- wilcox.test(Age_Acceleration ~ Invasiveness, data = data_subset)
    data.frame(
        Source = src,
        p.value = test$p.value,
        label = format_pval(test$p.value),
        # Calculate y position for annotation (above the boxplots)
        y.position = max(data_subset$Age_Acceleration, na.rm = TRUE) + 2
    )
})
stats_df <- do.call(rbind, stats_list)

# maube fix the p value?
p <- ggplot(ss, aes(x = Source, y = Age_Acceleration, fill = Invasiveness)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                alpha = 0.3, size = 0.6) +
    scale_fill_manual(values = invasiveness_colors) +
    theme_minimal() +
    labs(x = "Cohort",
         y = "Age Acceleration (Years)") +
    theme(axis.text.x = element_text(hjust = 1)) +
    geom_text(data = stats_df,
              aes(x = c(1, 2), y = 75, label = label),
              inherit.aes = FALSE,
              size = 3)
ggsave(here("figures", "20250204_thyroid88_age_acceleration_invasiveness_cohort.pdf"),
       plot = p, width=5, height=4.5)

# Print the statistical results
print("Wilcoxon test results by cohort:")
print(stats_df)

# Source   p.value     label y.position
# 1    PED 0.0385872 p = 0.039   43.19281
# 2   THCA 0.2089089 p = 0.209   71.14655
