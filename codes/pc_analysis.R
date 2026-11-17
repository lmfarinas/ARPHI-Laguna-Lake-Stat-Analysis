library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)
library(multcompView)
library(patchwork)
library(lubridate)
library(mice)
library(ggmice)
library(pan)
library(psych)
library(vegan)

library(lme4)
library(lmerTest)
library(emmeans)

set.seed(2026)

# Physicochemical parameters, exclude nitrate
pc_vars_all <- c(
  "Nitrate", "Ammonia", "BOD", "DO", "FC", "pH", "Phosphate"
)

pc_vars_raw <- c(
  "Ammonia", "BOD", "DO", "FC", "pH", "Phosphate"
)
# Parameters to log10-transform for multivariate analyses
pc_vars_log <- c(
  "Ammonia", "BOD", "FC", "Phosphate"
)

# Primary multivariate analysis (All six parameters)
pc_vars_primary <- c(
  "log10_Ammonia", "log10_BOD", "DO",
  "log10_FC", "pH", "log10_Phosphate"
)

# five-parameter sensitivity analysis excludes DO because of substantial
# regionally patterned censoring.
pc_vars_five <- c(
  "log10_Ammonia", "log10_BOD",
  "log10_FC", "pH", "log10_Phosphate"
)

n_imputations <- 20
n_iterations <- 20
n_permutations <- 9999

# ----- Load Data -----

sample_metadata_complete <- read.delim(
  "sample_metadata.tsv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    sample_collection_date = as.Date(sample_collection_date),
    Year = year(sample_collection_date)
  )


# sample and metadata checks
cat("Number of water samples:", nrow(sample_metadata_complete), "\n")
cat("Number of sites:", n_distinct(sample_metadata_complete$Site), "\n")
cat("Number of sampling dates:",
    n_distinct(sample_metadata_complete$sample_collection_date), "\n")


# ----- CLASSIFY QUANTIFIED, BDL, AND MISSING VALUES -----
# includes Nitrate
pc_all_long <- sample_metadata_complete %>%
  select(
    campaign_id, water_sample_id, Site, Region, Season,
    sample_collection_date, Year,
    all_of(pc_vars_all)
  ) %>%
  mutate(
    across(all_of(pc_vars_all), as.character)
  ) %>%
  pivot_longer(
    cols = all_of(pc_vars_all),
    names_to = "Parameter",
    values_to = "Raw_Value"
  ) %>%
  mutate(
    Raw_Value = str_trim(Raw_Value),
    Status = case_when(
      is.na(Raw_Value) | Raw_Value == "***" ~ "Missing",
      str_starts(Raw_Value, "<") ~ "BDL",
      TRUE ~ "Quantified"
    ),
    Censoring_Limit = if_else(
      Status == "BDL",
      as.numeric(str_remove(Raw_Value, "^<")),
      NA_real_
    ),
    Value = if_else(
      Status == "Quantified",
      as.numeric(Raw_Value),
      NA_real_
    )
  ) %>%
  mutate(pc_row_id = row_number())

# excludes Nitrate
pc_long <- sample_metadata_complete %>%
  select(
    campaign_id, water_sample_id, Site, Region, Season,
    sample_collection_date, Year,
    all_of(pc_vars_raw)
  ) %>%
  mutate(
    across(all_of(pc_vars_raw), as.character)
  ) %>%
  pivot_longer(
    cols = all_of(pc_vars_raw),
    names_to = "Parameter",
    values_to = "Raw_Value"
  ) %>%
  mutate(
    Raw_Value = str_trim(Raw_Value),
    Status = case_when(
      is.na(Raw_Value) | Raw_Value == "***" ~ "Missing",
      str_starts(Raw_Value, "<") ~ "BDL",
      TRUE ~ "Quantified"
    ),
    Censoring_Limit = if_else(
      Status == "BDL",
      as.numeric(str_remove(Raw_Value, "^<")),
      NA_real_
    ),
    Value = if_else(
      Status == "Quantified",
      as.numeric(Raw_Value),
      NA_real_
    )
  ) %>%
  mutate(pc_row_id = row_number())

# BDL counts by parameter and reporting limit
bdl_summary <- pc_all_long %>%
  filter(Status == "BDL") %>%
  count(Parameter, Censoring_Limit, name = "n") %>%
  arrange(Parameter, Censoring_Limit) %>%
  rename(reporting_limit = Censoring_Limit)

bdl_summary_wide <- pc_all_long %>%
  count(Region, Parameter, name = "n_total") %>%
  left_join(
    pc_all_long %>%
      filter(Status == "BDL") %>%
      count(Region, Parameter, name = "n_BDL"),
    by = c("Region", "Parameter")
  ) %>%
  replace_na(list(n_BDL = 0)) %>%
  mutate("n_BDL/n_total" = paste(n_BDL, n_total, sep = "/")) %>%
  select(-n_BDL, -n_total) %>%
  pivot_wider(
    names_from = Region, 
    values_from = "n_BDL/n_total"
  ) %>%
  left_join(
    pc_all_long %>%
      filter(Status == "BDL") %>%
      count(Parameter, name = "n_BDL"),
    by = "Parameter"
  ) %>%
  replace_na(list(n_BDL = 0)) %>%
  left_join(
    pc_all_long %>%
      count(Season, Parameter, name = "n_total") %>%
      left_join(
        pc_all_long %>%
          filter(Status == "BDL") %>%
          count(Season, Parameter, name = "n_BDL"),
        by = c("Season", "Parameter")
      ) %>%
      replace_na(list(n_BDL = 0)) %>%
      mutate("n_BDL/n_total" = paste(n_BDL, n_total, sep = "/")) %>%
      select(-n_BDL, -n_total) %>%
      pivot_wider(
        names_from = Season, 
        values_from = "n_BDL/n_total"
      ),
    by = "Parameter"
  ) %>%
  relocate(n_BDL, "East Lake", "West Lake", .after = Parameter)

print(bdl_summary_wide)

# Data-quality summary
data_quality_summary <- pc_long %>%
  group_by(Parameter) %>%
  summarise(
    n_total = n(),
    n_quantified = sum(Status == "Quantified"),
    n_BDL = sum(Status == "BDL"),
    n_missing = sum(Status == "Missing"),
    .groups = "drop"
  ) %>%
  mutate(
    quantified_pct = round(n_quantified/n_total*100,2),
    BDL_pct = round(n_BDL/n_total*100,2),
    missing_pct = round(n_missing/n_total*100,2)
  )
  

print(data_quality_summary)

# Check distributions of quantified values
pc_long %>%
  filter(Status == "Quantified") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 25, na.rm = TRUE) +
  facet_wrap(~ Parameter, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Distributions of quantified physicochemical measurements",
    x = "Reported measurement",
    y = "Number of observations"
  )

# ----- LIMIT/SQRT(2) SUBSTITUTION -----
pc_limit_sqrt2_long <- pc_long %>%
  mutate(
    Value_processed = case_when(
      Status == "BDL" ~ Censoring_Limit / sqrt(2),
      Status == "Quantified" ~ Value,
      TRUE ~ NA_real_
    )
  )

# Verify BDL substitution and remaining missing values
pc_limit_sqrt2_long %>%
  group_by(Parameter) %>%
  summarise(
    n_total = n(),
    n_quantified = sum(Status == "Quantified"),
    n_BDL = sum(Status == "BDL"),
    n_missing = sum(Status == "Missing"),
    n_processed = sum(!is.na(Value_processed)),
    .groups = "drop"
  )

# Convert to one row per water sample
pc_limit_sqrt2_wide <- pc_limit_sqrt2_long %>%
  select(
    campaign_id, water_sample_id, Site, Region, Season,
    sample_collection_date, Year,
    Parameter, Value_processed
  ) %>%
  pivot_wider(
    names_from = Parameter,
    values_from = Value_processed
  ) %>%
  arrange(sample_collection_date, Site)

# Confirm missingness before imputation
print(colSums(is.na(pc_limit_sqrt2_wide[, pc_vars_raw])))


# ----- MISSING-DATA IMPUTATION -----

# Missing values occur in Ammonia and BOD in the current dataset.
# Imputation is performed on log10-transformed values to keep
# imputed concentrations positive after back-transformation.

mice_data <- pc_limit_sqrt2_wide %>%
  mutate(
    Site_cluster = as.integer(factor(Site)),
    Season_WET = as.integer(Season == "WET"),
    Year_c = Year - mean(Year, na.rm = TRUE)
  ) %>%
  mutate(
    log10_Ammonia = log10(Ammonia),
    log10_BOD = log10(BOD),
    log10_FC = log10(FC),
    log10_Phosphate = log10(Phosphate)
  ) %>%
  select(-all_of(pc_vars_log))

imputation_vars <- c(
  "log10_Ammonia", "log10_BOD",
  "DO", "log10_FC", "pH", "log10_Phosphate"
)

# Standardize continuous imputation variables for numerical stability
imputation_center <- sapply(
  mice_data[imputation_vars],
  mean,
  na.rm = TRUE
)

imputation_scale <- sapply(
  mice_data[imputation_vars],
  sd,
  na.rm = TRUE
)

mice_data_scaled <- mice_data %>%
  mutate(
    across(
      all_of(imputation_vars),
      ~ (.x - imputation_center[cur_column()]) /
        imputation_scale[cur_column()]
    )
  ) %>%
  select(
    all_of(imputation_vars),
    Season_WET,
    Year_c,
    Site_cluster
  )

## Primary MICE ####
# PRIMARY: TWO-LEVEL MICE WITH SITE CLUSTERING

# 2l.pan requires the cluster variable to be coded -2 in the
# predictor matrix. Fixed-effect predictors are coded 1.
ini_cluster <- mice(
  mice_data_scaled,
  maxit = 0,
  printFlag = FALSE
)

meth_cluster <- ini_cluster$method
pred_cluster <- ini_cluster$predictorMatrix

meth_cluster[] <- ""
meth_cluster[c("log10_Ammonia", "log10_BOD")] <- "2l.pan"

pred_cluster[,] <- 0

# Site random-intercept grouping variable
pred_cluster[
  c("log10_Ammonia", "log10_BOD"),
  "Site_cluster"
] <- -2

# Fixed-effect predictors
cluster_predictors <- c(
  "log10_Ammonia", "log10_BOD",
  "DO", "log10_FC", "pH", "log10_Phosphate",
  "Season_WET", "Year_c"
)

for (target in c("log10_Ammonia", "log10_BOD")) {
  pred_cluster[
    target,
    setdiff(cluster_predictors, target)
  ] <- 1
}

set.seed(2026)

imp_cluster <- mice(
  mice_data_scaled,
  m = n_imputations,
  maxit = n_iterations,
  method = meth_cluster,
  predictorMatrix = pred_cluster,
  printFlag = TRUE
)

## Sensitivity MICE ####
# SENSITIVITY 1: SIMPLE PMM MICE, WITHOUT CLUSTERING

# Same log-transformed and standardized imputation variables,
# but no Site, Season, or Year predictors.
mice_data_simple <- mice_data_scaled %>%
  select(all_of(imputation_vars))

ini_simple <- mice(
  mice_data_simple,
  maxit = 0,
  printFlag = FALSE
)

meth_simple <- ini_simple$method
pred_simple <- ini_simple$predictorMatrix

meth_simple[] <- ""
meth_simple[c("log10_Ammonia", "log10_BOD")] <- "pmm"
pred_simple[,] <- 0

for (target in c("log10_Ammonia", "log10_BOD")) {
  pred_simple[
    target,
    setdiff(imputation_vars, target)
  ] <- 1
}

set.seed(2026)

imp_simple <- mice(
  mice_data_simple,
  m = n_imputations,
  maxit = n_iterations,
  method = meth_simple,
  predictorMatrix = pred_simple,
  printFlag = TRUE,
  ridge = 1e-05
)

## Diagnostics for MICE ####

# Logged events
print(imp_cluster$loggedEvents)
print(imp_simple$loggedEvents)

### Trace / Convergence plots ####

# Save file for imp_cluster
png(filename = "supplementary_files/MICE_cluster.png", width = 7, height = 4, units = "in", res = 600)
plot(imp_cluster)
dev.off()

# Save file for imp_simple
png(filename = "supplementary_files/MICE_simple.png", width = 7, height = 4, units = "in", res = 600)
plot(imp_simple)
dev.off()

### Density plots ####
plot_density <- function(data) {
  # Create the Ammonia plot
  p_ammonia <- ggmice(data, aes(x = .data[["log10_Ammonia"]], group = .imp)) +
    geom_density() +
    scale_color_manual(
      name = "Data",
      values = c("observed" = "blue", "imputed" = "indianred3")
    ) +
    theme(legend.position = "bottom") +
    labs(
      x = "log10_Ammonia (standardized)",
      y = "Density"
    )
  
  # Create the BOD plot
  p_bod <- ggmice(data, aes(x = .data[["log10_BOD"]], group = .imp)) +
    geom_density() +
    scale_color_manual(
      name = "Data",
      values = c("observed" = "blue", "imputed" = "indianred3")
    ) +
    theme(legend.position = "bottom") +
    labs(
      x = "log10_BOD (standardized)",
      y = "Density"
    )
  
  # Combine plots with patchwork and collect legends
  combined_plot <- p_ammonia + p_bod +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  return(combined_plot)
}

# Save file for imp_cluster
p_cluster_density_bod_ammonia <- plot_density(imp_cluster)
ggsave("supplementary_files/MICE_cluster_densityplot.png",
       plot = p_cluster_density_bod_ammonia,
       width = 6, height = 3.5, dpi = 600)

# Save file for imp_simple
p_simple_density_bod_ammonia <- plot_density(imp_simple)
ggsave("supplementary_files/MICE_simple_densityplot.png",
       plot = p_simple_density_bod_ammonia,
       width = 6, height = 3.5, dpi = 600)

### Strip plots ####
plot_stripplot <- function(data) {
  # Create the Ammonia plot
  p_ammonia <- ggmice(data, aes(x = .imp, y = log10_Ammonia)) + 
    geom_jitter(width = 0.25, height = 0) +
    labs(x = "Iteration", y = "log10_Ammonia (standardized)") +
    scale_color_manual(
      name = "Data",
      values = c("observed" = "blue", "imputed" = "indianred3")
    )
  
  # Create the BOD plot
  p_bod <- ggmice(data, aes(x = .imp, y = log10_BOD)) + 
    geom_jitter(width = 0.25, height = 0) +
    labs(x = "Iteration", y = "log10_BOD (standardized)") +
    scale_color_manual(
      name = "Data",
      values = c("observed" = "blue", "imputed" = "indianred3")
    )
  
  # Combine plots with patchwork and collect legends
  combined_plot <- p_ammonia + p_bod +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  return(combined_plot)
}

p_cluster_strip_bod_ammonia <- plot_stripplot(imp_cluster)
ggsave("supplementary_files/MICE_cluster_stripplot.png",
       plot = p_cluster_strip_bod_ammonia,
       width = 7, height = 3.5, dpi = 600)

p_simple_strip_bod_ammonia <- plot_stripplot(imp_simple)
ggsave("supplementary_files/MICE_simple_stripplot.png",
       plot = p_simple_strip_bod_ammonia,
       width = 7, height = 3.5, dpi = 600)

## Back transform imputations ####
# Helper function: reverse standardization and log transformation
back_transform_imputations <- function(imp_object) {

  completed <- complete(imp_object, action = "all")

  lapply(completed, function(dat) {

    dat_original_scale <- dat %>%
      mutate(
        across(
          all_of(imputation_vars),
          ~ .x * imputation_scale[cur_column()] +
            imputation_center[cur_column()]
        )
      )

    dat_original_scale %>%
      transmute(
        Ammonia = 10^log10_Ammonia,
        BOD = 10^log10_BOD,
        DO = DO,
        FC = 10^log10_FC,
        pH = pH,
        Phosphate = 10^log10_Phosphate
      )
  })
}

# Back-transform completed datasets
cluster_completed <- back_transform_imputations(imp_cluster)
simple_completed <- back_transform_imputations(imp_simple)

# Attach sample metadata
attach_metadata <- function(completed_list) {

  lapply(completed_list, function(dat) {
    bind_cols(
      pc_limit_sqrt2_wide %>%
        select(
          campaign_id, water_sample_id, Site, Region, Season,
          sample_collection_date, Year
        ),
      dat
    )
  })
}

pc_mice_cluster <- attach_metadata(cluster_completed)
pc_mice_simple <- attach_metadata(simple_completed)

# Check imputation outputs
tibble(
  Dataset = c("Cluster-aware MICE", "Simple MICE"),
  n_imputations = c(
    length(pc_mice_cluster),
    length(pc_mice_simple)
  ),
  n_samples = c(
    nrow(pc_mice_cluster[[1]]),
    nrow(pc_mice_simple[[1]])
  ),
  n_missing = c(
    sum(is.na(pc_mice_cluster[[1]][pc_vars_raw])),
    sum(is.na(pc_mice_simple[[1]][pc_vars_raw]))
  )
)

# ----- Define Analysis Datasets -----
# Primary: Limit/sqrt(2) + Site-clustered MICE
pc_analysis_primary <- pc_mice_cluster

# Sensitivity 1: Limit/sqrt(2) + simple MICE
pc_analysis_simple <- pc_mice_simple

# Sensitivity 2: Limit/sqrt(2) + complete-case
pc_complete <- pc_limit_sqrt2_wide %>%
  drop_na(all_of(pc_vars_raw))

cat("Primary MICE sample size:", nrow(pc_analysis_primary[[1]]), "\n")
cat("Simple MICE sample size:", nrow(pc_analysis_simple[[1]]), "\n")
cat("Complete-case sample size:", nrow(pc_complete), "\n")


## Transform Data for Multivariate Analysis ####

transform_pc_data <- function(dat) {

  dat %>%
    mutate(
      log10_Ammonia = log10(Ammonia),
      log10_BOD = log10(BOD),
      log10_FC = log10(FC),
      log10_Phosphate = log10(Phosphate)
    )
}

pc_primary_transformed <- lapply(
  pc_analysis_primary,
  transform_pc_data
)

pc_simple_transformed <- lapply(
  pc_analysis_simple,
  transform_pc_data
)

pc_complete_transformed <- transform_pc_data(pc_complete)


# ----- Exploratory PCA -----
# PCA is performed on the six parameters.
# pca_primary_five = Sensitivity analysis with DO excluded

## Run PCA ####
pca_primary <- lapply(
  pc_primary_transformed,
  function(dat) {
    prcomp(
      dat %>% select(all_of(pc_vars_primary)),
      center = TRUE,
      scale. = TRUE
    )
  }
)

pca_primary_five <- lapply(
  pc_primary_transformed,
  function(dat) {
    prcomp(
      dat %>% select(all_of(pc_vars_five)),
      center = TRUE,
      scale. = TRUE
    )
  }
)

pca_simple <- lapply(
  pc_simple_transformed,
  function(dat) {
    prcomp(
      dat %>% select(all_of(pc_vars_primary)),
      center = TRUE,
      scale. = TRUE
    )
  }
)

pca_complete <- prcomp(
  pc_complete_transformed %>%
    select(all_of(pc_vars_primary)),
  center = TRUE,
  scale. = TRUE
)

## Check explained variance across imputations #### 
summarise_pca_variance <- function(pca_list, dataset_name) {

  bind_rows(
    lapply(seq_along(pca_list), function(i) {
      tibble(
        Imputation = i,
        PC = paste0("PC", seq_along(pca_list[[i]]$sdev)),
        Variance = pca_list[[i]]$sdev^2 /
          sum(pca_list[[i]]$sdev^2)
      )
    })
  ) %>%
    group_by(PC) %>%
    summarise(
      Mean = mean(Variance),
      SD = sd(Variance),
      Min = min(Variance),
      Max = max(Variance),
      .groups = "drop"
    ) %>%
    mutate(Dataset = dataset_name)
}

pca_variance_primary <- summarise_pca_variance(
  pca_primary, "Cluster-aware MICE"
)

pca_variance_primary_five <- summarise_pca_variance(
  pca_primary_five, "Cluster-aware MICE"
)

pca_variance_simple <- summarise_pca_variance(
  pca_simple, "Simple MICE"
)

pca_variance_complete <- tibble(
  PC = paste0("PC", seq_along(pca_complete$sdev)),
  Variance = pca_complete$sdev^2 / sum(pca_complete$sdev^2),
  Dataset = "Complete-case"
)

print(pca_variance_primary)
print(pca_variance_primary_five)
print(pca_variance_simple)
print(pca_variance_complete)

## Check explained variance and loadings #### 
# Extract variance explained across imputations
extract_var <- function(pca_list, model_name) {
  map_dfr(seq_along(pca_list), function(i) {
    var_pct <- (pca_list[[i]]$sdev^2) / sum(pca_list[[i]]$sdev^2) * 100
    tibble(
      model = model_name,
      imp = i,
      PC1 = var_pct[1],
      PC2 = var_pct[2]
    )
  })
}

var_df <- bind_rows(
  extract_var(pca_primary, "Primary MICE (w/ DO)"),
  extract_var(pca_primary_five, "Primary MICE (w/o DO)"),
  extract_var(pca_simple, "Simple MICE (w/ DO)")
) %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "var_explained") %>%
  mutate(
    model = factor(
      model, levels = c("Primary MICE (w/o DO)",
                        "Primary MICE (w/ DO)",
                        "Simple MICE (w/ DO)")
      ))

# Complete case variance reference
comp_var_pct <- (pca_complete$sdev^2) / sum(pca_complete$sdev^2) * 100
var_complete_df <- tibble(
  PC = c("PC1", "PC2"),
  var_explained = c(comp_var_pct[1], comp_var_pct[2])
)

# Extract loadings across imputations
extract_loadings <- function(pca_list, model_name) {
  map_dfr(seq_along(pca_list), function(i) {
    as.data.frame(pca_list[[i]]$rotation[, 1:2]) %>%
      rownames_to_column("variable") %>%
      mutate(imp = i, model = model_name)
  })
}

loadings_df <- bind_rows(
  extract_loadings(pca_primary, "Primary MICE (w/ DO)"),
  extract_loadings(pca_primary_five, "Primary MICE (w/o DO)"),
  extract_loadings(pca_simple, "Simple MICE (w/ DO)")
) %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "loading") %>%
  mutate(
    model = factor(
      model, levels = c("Primary MICE (w/o DO)",
                        "Primary MICE (w/ DO)",
                        "Simple MICE (w/ DO)")
    ))

# Complete case loadings reference
loadings_complete_df <- as.data.frame(pca_complete$rotation[, 1:2]) %>%
  rownames_to_column("variable") %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "loading")

plot_variance <- ggplot(var_df, aes(x = model, y = var_explained, fill = model)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7, aes(color = model)) +
  geom_hline(
    data = var_complete_df, 
    aes(yintercept = var_explained), 
    linetype = "dashed", 
    color = "black", 
    linewidth = 0.7
  ) +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    #title = "Variance Explained Stability Across 20 Imputations",
    #subtitle = "Dashed black line indicates Non-imputed complete case dataset",
    x = "Imputation model",
    y = "% Variance Explained"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "gray80", fill = NA)
  )

print(plot_variance)

ggsave("supplementary_files/PCA_var_explained.png",
       plot = plot_variance,
       width = 7, height = 3.5, dpi = 600)

plot_loadings <- ggplot() +
  # Imputed models: Mean loading with Min-Max range bars across 20 imputations
  stat_summary(
    data = loadings_df,
    aes(x = loading, y = variable, color = model),
    fun.data = median_hilow,
    fun.args = list(conf.int = 1),
    geom = "errorbar",
    width = 0.6,        # Height of the end-cap ticks
    linewidth = 1,    # Line thickness
    position = position_dodge(width = 0.6)
  ) +
  stat_summary(
    data = loadings_df,
    aes(x = loading, y = variable, color = model),
    fun.y = median,
    geom = "point",
    size = 2,
    position = position_dodge(width = 0.6)
  ) +
  # Complete case reference points (red diamonds)
  geom_point(
    data = loadings_complete_df,
    aes(x = loading, y = variable, shape = "Complete Case"),
    color = "firebrick",
    size = 2,
    stroke = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  facet_grid(model ~ PC, scales = "free_y") +
  scale_color_brewer(palette = "Set2", name = "Imputed Model") +
  scale_shape_manual(name = "Non-imputed", values = c("Complete Case" = 18)) +
  labs(
    #title = "Factor Loadings Stability Across 20 Imputations",
    #subtitle = "Points & range bars show 20-imputation distribution; Diamonds show Complete Case loadings",
    x = "Factor Loading",
    y = "Variable"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    #legend.position = "bottom",
    panel.border = element_rect(color = "gray80", fill = NA)
  )

print(plot_loadings)

ggsave("supplementary_files/PCA_factor_loadings.png",
       plot = plot_loadings,
       width = 7, height = 6, dpi = 600)

## KMO and Bartlett tests ####

calculate_PCA_adequacy <- function(data, vars) {
  
  pc_data <- data %>%
    select(all_of(vars))
  
  # Correlation matrix
  R <- cor(pc_data, use = "complete.obs")
  
  # KMO test
  kmo <- psych::KMO(R)
  
  # Bartlett's test of sphericity
  bartlett <- psych::cortest.bartlett(
    R,
    n = nrow(pc_data)
  )
  
  tibble(
    n = nrow(pc_data),
    Overall_MSA = kmo$MSA,
    !!!setNames(
      as.list(kmo$MSAi),
      paste0("MSA_", names(kmo$MSAi))
    ),
    Bartlett_chisq = bartlett$chisq,
    Bartlett_df = bartlett$df,
    Bartlett_p = bartlett$p.value
  )
}

# Primary MICE (with DO)
primary_adequacy <- map_dfr(
  seq_along(pc_primary_transformed),
  ~ calculate_PCA_adequacy(
    pc_primary_transformed[[.x]],
    pc_vars_primary
  ) %>%
    mutate(imp = .x)
)

# Primary MICE (without DO)
primary_five_adequacy <- map_dfr(
  seq_along(pc_primary_transformed),
  ~ calculate_PCA_adequacy(
    pc_primary_transformed[[.x]],
    pc_vars_five
  ) %>%
    mutate(imp = .x)
)

# Simple MICE (with DO)
simple_adequacy <- map_dfr(
  seq_along(pc_simple_transformed),
  ~ calculate_PCA_adequacy(
    pc_simple_transformed[[.x]],
    pc_vars_primary
  ) %>%
    mutate(imp = .x)
)

# Complete-case dataset (with DO)
complete_adequacy <- calculate_PCA_adequacy(
  pc_complete_transformed,
  pc_vars_primary
)

# Format median and range
format_median_range <- function(x, digits = 3) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) return(NA_character_)
  
  if (length(x) == 1) {
    return(sprintf(
      paste0("%.", digits, "f"),
      x
    ))
  }
  
  sprintf(
    paste0("%.", digits, "f [%.", digits,
           "f–%.", digits, "f]"),
    median(x),
    min(x),
    max(x)
  )
}

format_p_summary <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 1) {
    return(sprintf("%.2e", x))
  }
  
  sprintf(
    "%.2e [%.2e–%.2e]",
    median(x),
    min(x),
    max(x)
  )
}

# Summarize an imputed dataset
summarize_adequacy_table <- function(data, model_name) {
  
  numeric_cols <- data %>%
    select(-imp) %>%
    select(where(is.numeric)) %>%
    names()
  
  map_dfc(numeric_cols, function(v) {
    
    x <- data[[v]]
    
    if (v == "Bartlett_p") {
      formatted <- format_p_summary(x)
      
    } else {
      digits <- if (v %in% c("n", "Bartlett_df")) 0 else 3
      formatted <- format_median_range(x, digits)
    }
    
    tibble(!!v := formatted)
    
  }) %>%
    mutate(Model = model_name, .before = 1)
}

# Summarize all approaches
adequacy_summary <- bind_rows(
  
  summarize_adequacy_table(
    primary_adequacy,
    "Primary MICE (with DO)"
  ),
  
  summarize_adequacy_table(
    primary_five_adequacy,
    "Primary MICE (without DO)"
  ),
  
  summarize_adequacy_table(
    simple_adequacy,
    "Simple MICE (with DO)"
  ),
  
  summarize_adequacy_table(
    complete_adequacy %>% mutate(imp = 1),
    "Complete-case (with DO)"
  )
  
)

# Arrange columns
adequacy_summary <- adequacy_summary %>%
  select(
    Model,
    n,
    Overall_MSA,
    starts_with("MSA_"),
    Bartlett_chisq,
    Bartlett_df,
    Bartlett_p
  )

# Export files
write_tsv(
  adequacy_summary,
  "supplementary_files/Supplementary_KMO_Bartlett_summary.tsv")

## Plot PCA ####
# PCA scores for the first primary imputation.
# This plot describes the primary imputed dataset only; it does not
# compare independently fitted PCA coordinate systems.

plot_PCA <- function(
    model = c("primary", "primary_five", "simple", "complete"),
    iteration = 1,
    scale_factor = 2,
    show_ellipses = TRUE
) {
  
  model <- match.arg(model)
  
  # Select PCA object and corresponding transformed dataset
  pca_name <- switch(
    model,
    primary = "pca_primary",
    primary_five = "pca_primary_five",
    simple = "pca_simple",
    complete = "pca_complete"
  )
  
  data_name <- switch(
    model,
    primary = "pc_primary_transformed",
    primary_five = "pc_primary_transformed",
    simple = "pc_simple_transformed",
    complete = "pc_complete_transformed"
  )
  
  pca_obj <- get(pca_name, envir = parent.frame())
  data_obj <- get(data_name, envir = parent.frame())
  
  # Handle MICE lists versus complete-case PCA
  if (model == "complete") {
    
    if (iteration != 1) {
      stop("Complete-case PCA only has one dataset. Use iteration = 1.")
    }
    
    pca_fit <- pca_obj
    dat <- data_obj
    
  } else {
    
    if (iteration < 1 || iteration > length(pca_obj)) {
      stop(
        paste0(
          "Invalid iteration. Choose between 1 and ",
          length(pca_obj), "."
        )
      )
    }
    
    pca_fit <- pca_obj[[iteration]]
    dat <- data_obj[[iteration]]
  }
  
  # Prepare PCA scores
  scores <- as.data.frame(pca_fit$x) %>%
    bind_cols(
      dat %>%
        select(water_sample_id, Site, Region, Season)
    ) %>%
    mutate(
      Region = factor(
        Region,
        levels = c(
          "East Lake",
          "West Lake",
          "Quadrat 1",
          "Quadrat 2",
          "Quadrat 3",
          "Quadrat 4"
        )
      ),
      Season = factor(Season)
    )
  
  # Prepare PCA loadings
  loadings <- as.data.frame(pca_fit$rotation) %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(
      PC1 = PC1 * scale_factor,
      PC2 = PC2 * scale_factor
    )
  
  # Calculate variance explained for selected PCA
  variance <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  
  # Model descriptions
  model_labels <- c(
    primary = "Primary analysis \n(Cluster-aware MICE)",
    primary_five = "Sensitivity analysis \n(DO excluded)",
    simple = "Sensitivity analysis \n(Simple MICE)",
    complete = "Sensitivity analysis \n(Complete-case)"
  )
  
  subtitle <- if (model == "complete") {
    model_labels[model]
  } else {
    paste0(
      model_labels[model],
      " - Imputation ", iteration
    )
  }
  
  # Generate plot
  p <- ggplot(
    scores,
    aes(x = PC1, y = PC2, color = Region, shape = Season)
  )
  
  # Add regional ellipses
  if (show_ellipses) {
    p <- p +
      stat_ellipse(
        data = scores,
        aes(
          x = PC1,
          y = PC2,
          group = Region,
          color = Region,
          fill = Region
        ),
        inherit.aes = FALSE,
        geom = "polygon",
        alpha = 0.10,
        linewidth = 0.8,
        level = 0.95,
        show.legend = FALSE
      )
  }
  
  # Add points and loading vectors
  p <- p +
    geom_point(
      size = 2.5,
      alpha = 0.8
    ) +
    geom_segment(
      data = loadings,
      aes(
        x = 0, y = 0,
        xend = PC1, yend = PC2
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(
        length = grid::unit(0.25, "cm")
      ),
      color = "black",
      linewidth = 0.6
    ) +
    geom_text(
      data = loadings,
      aes(
        x = PC1,
        y = PC2,
        label = Parameter
      ),
      inherit.aes = FALSE,
      color = "red",
      size = 3.5,
      vjust = -0.5
    ) +
    theme_minimal() +
    labs(
      #title = "PCA of physicochemical water quality parameters",
      subtitle = subtitle,
      x = paste0(
        "PC1 (",
        round(variance[1] * 100, 1),
        "% variance)"
      ),
      y = paste0(
        "PC2 (",
        round(variance[2] * 100, 1),
        "% variance)"
      ),
      color = "Region",
      shape = "Season"
    )
  
  return(p)
}
### run function ####
# Generate individual plots
p_primary <- plot_PCA(model = "primary", iteration = 1)

p_primary_five <- plot_PCA(model = "primary_five", iteration = 1)

p_simple <- plot_PCA(model = "simple", iteration = 1)

p_complete <- plot_PCA(model = "complete")

# Combine plots
p_PCA_combined <- (
  (p_primary | p_primary_five) /
    (p_simple | p_complete)
) +
  plot_layout(
    guides = "collect"
  )

p_PCA_combined

ggsave("supplementary_files/PCA_plots.png",
       plot = p_PCA_combined,
       width = 9, height = 7, dpi = 600)

ggsave("supplementary_files/PCA_primary_plot.png",
       plot = p_primary,
       width = 7, height = 5, dpi = 600)

# ----- LMM on PC1 -----
## Extract PC1 Scores ####

region_levels <- c(
  "East Lake",
  "West Lake",
  "Quadrat 1",
  "Quadrat 2",
  "Quadrat 3",
  "Quadrat 4"
)

pca_dataset_config <- list(
  primary = list(
    pca = pca_primary,
    data = pc_primary_transformed,
    multiple_imputations = TRUE
  ),
  primary_five = list(
    pca = pca_primary_five,
    data = pc_primary_transformed,
    multiple_imputations = TRUE
  ),
  simple = list(
    pca = pca_simple,
    data = pc_simple_transformed,
    multiple_imputations = TRUE
  ),
  complete = list(
    pca = pca_complete,
    data = pc_complete_transformed,
    multiple_imputations = FALSE
  )
)

extract_pc1_scores <- function(pca_fit, transformed_data) {
  
  as.data.frame(pca_fit$x) %>%
    bind_cols(
      transformed_data %>%
        select(
          water_sample_id,
          campaign_id,
          Site,
          Region,
          Season
        )
    ) %>%
    select(
      water_sample_id,
      campaign_id,
      Site,
      Region,
      Season,
      PC1
    ) %>%
    mutate(
      Region = factor(Region, levels = region_levels),
      Season = factor(Season),
      Site = factor(Site),
      campaign_id = factor(campaign_id)
    )
}


# Extract scores for all datasets
extract_dataset_scores <- function(dataset_config) {
  
  if (dataset_config$multiple_imputations) {
    
    map2(
      dataset_config$pca,
      dataset_config$data,
      extract_pc1_scores
    )
    
  } else {
    
    list(
      extract_pc1_scores(
        dataset_config$pca,
        dataset_config$data
      )
    )
  }
}


# Apply extraction to all four datasets
pc1_scores <- map(
  pca_dataset_config,
  extract_dataset_scores
)

# Access individual datasets:
# pc1_scores$primary
# pc1_scores$primary_five
# pc1_scores$simple
# pc1_scores$complete

### Fit LMMs ####
pc1_formula <- PC1 ~ Season * Region + (1 | Site)


fit_lmm_across_imputations <- function(
    scores_list,
    formula = pc1_formula
) {
  
  model_fits <- map(
    scores_list,
    function(df) {
      fit <- lmer(
        formula,
        data = df,
        REML = FALSE
      )
      # Embed the data object into the call so emmeans can evaluate it anywhere
      fit@call$data <- df
      fit
    }
  )
  
  if (length(model_fits) > 1) {
    mice::as.mira(model_fits)
  } else {
    model_fits[[1]]
  }
}


# Fit the full model for all datasets
pc1_models <- map(
  pc1_scores,
  fit_lmm_across_imputations
)

### Pool fixed-effect estimates ####
pc1_fixed_effects <- imap(
  pc1_models,
  function(model, dataset_name) {
    
    if (dataset_name == "complete") {
      
      broom.mixed::tidy(
        model,
        effects = "fixed",
        conf.int = TRUE
      )
      
    } else {
      
      mice::pool(model) %>%
        summary(conf.int = TRUE)
    }
  }
)

# Example:
pc1_fixed_effects$primary
pc1_fixed_effects$primary_five
pc1_fixed_effects$simple
pc1_fixed_effects$complete

# Write out as supplementary
write_tsv(
  as.data.frame(pc1_fixed_effects$primary),
  "supplementary_files/Supplementary_pc1_fixed_effects_primary.tsv",
)
write_tsv(
  as.data.frame(pc1_fixed_effects$primary_five),
  "supplementary_files/Supplementary_pc1_fixed_effects_primary_five.tsv",
)
write_tsv(
  as.data.frame(pc1_fixed_effects$simple),
  "supplementary_files/Supplementary_pc1_fixed_effects_simple.tsv",
)
write_tsv(
  as.data.frame(pc1_fixed_effects$complete),
  "supplementary_files/Supplementary_pc1_fixed_effects_complete.tsv",
)


### Omnibus tests for interaction and main effects ####
# Will compare
# Full vs. no interaction to check interaction
# No interaction vs. no Region to check region main effect
# No interaction vs. no Season to check season main effect

pc1_reduced_formulas <- list(
  no_interaction = PC1 ~ Season + Region + (1 | Site),
  no_region = PC1 ~ Season + (1 | Site),
  no_season = PC1 ~ Region + (1 | Site)
)


fit_reduced_lmm <- function(
    scores_list,
    formula
) {
  
  model_fits <- map(
    scores_list,
    function(df) {
      fit <- lmer(
        formula,
        data = df,
        REML = FALSE
      )
      # Embed the data object into the call so emmeans can evaluate it anywhere
      fit@call$data <- df
      fit
    }
  )
  
  if (length(model_fits) > 1) {
    mice::as.mira(model_fits)
  } else {
    model_fits[[1]]
  }
}


# Fit reduced models for all datasets
pc1_reduced_models <- purrr::map(
  pc1_reduced_formulas,
  ~ map(
    pc1_scores,
    fit_reduced_lmm,
    formula = .x
  )
)

compare_pc1_models <- function(
    full_model,
    reduced_model,
    multiply_imputed
) {
  
  if (multiply_imputed) {
    
    mice::D1(
      full_model,
      reduced_model
    )
    
  } else {
    
    anova(
      reduced_model,
      full_model
    )
  }
}


# Compare models for each dataset
pc1_omnibus_tests <- imap(
  pc1_scores,
  function(scores, dataset_name) {
    
    is_mi <- length(scores) > 1
    
    full_model <- pc1_models[[dataset_name]]
    
    no_interaction <- pc1_reduced_models$no_interaction[[dataset_name]]
    no_region <- pc1_reduced_models$no_region[[dataset_name]]
    no_season <- pc1_reduced_models$no_season[[dataset_name]]
    
    list(
      interaction = compare_pc1_models(
        full_model,
        no_interaction,
        is_mi
      ),
      
      region = compare_pc1_models(
        no_interaction,
        no_region,
        is_mi
      ),
      
      season = compare_pc1_models(
        no_interaction,
        no_season,
        is_mi
      )
    )
  }
)

# Example:
pc1_omnibus_tests$primary$region$result
pc1_omnibus_tests$primary_five
pc1_omnibus_tests$simple
pc1_omnibus_tests$complete$interaction

models <- c("primary", "primary_five", "simple")

# Bind all nested $result data frames into one master tibble
pc1_omnibus_tests_df <- map_dfr(models, function(m) {
  pc1_omnibus_tests[[m]] %>% 
    map_dfr(~ as_tibble(.x$result), .id = "effect") %>% 
    mutate(model = m, .before = 1)
})

pc1_omnibus_test_complete <- map_dfr(pc1_omnibus_tests$complete, function(eff) {
  # Extract object (handles whether it has a $result sub-element or is direct)
  res <- if (is.list(eff) && "result" %in% names(eff)) eff$result else eff
  
  # Take the 2nd row (full model LRT test comparison)
  as_tibble(res)[2, ]
}, .id = "effect") %>%
  mutate(model = "complete") %>%
  select(model, effect, Chisq, Df, "Pr(>Chisq)")


# Export supplementary
write_tsv(
  pc1_omnibus_tests_df,
  "supplementary_files/Supplementary_pc1_omnibus_tests.tsv",
)
write_tsv(
  pc1_omnibus_test_complete,
  "supplementary_files/Supplementary_pc1_omnibus_test_complete.tsv",
)


### Estimated Marginal Means ####
get_model_fits <- function(model) {
  
  if (inherits(model, "mira")) {
    model$analyses
  } else {
    list(model)
  }
}


calculate_pc1_emmeans <- function(model) {
  
  model_fits <- get_model_fits(model)
  
  region_emmeans <- map(
    model_fits,
    ~ emmeans(.x, ~ Region | Season)
  )
  
  season_emmeans <- map(
    model_fits,
    ~ emmeans(.x, ~ Season | Region)
  )
  
  list(
    region = region_emmeans,
    season = season_emmeans
  )
}


# Calculate EMMs for all datasets
pc1_emmeans <- map(
  pc1_models,
  calculate_pc1_emmeans
)

# Calculate pairwise comparisons
calculate_pc1_pairwise <- function(emmeans_results) {
  
  list(
    region = map(
      emmeans_results$region,
      ~ pairs(.x, adjust = "none")
    ),
    
    season = map(
      emmeans_results$season,
      ~ pairs(.x, adjust = "none")
    )
  )
}

pc1_pairwise <- map(
  pc1_emmeans,
  calculate_pc1_pairwise
)

### Pool EMMs and Pairwise Contrasts ####

pool_mi_estimates <- function(
    results_list,
    grouping_variables,
    estimate_column,
    df_column = "df"
) {
  
  all_results <- map_dfr(
    seq_along(results_list),
    function(iteration) {
      
      as.data.frame(results_list[[iteration]]) %>%
        mutate(imputation = iteration)
    }
  )
  
  pooled_results <- all_results %>%
    group_by(across(all_of(grouping_variables))) %>%
    group_modify(~ {
      
      estimates <- .x[[estimate_column]]
      variances <- .x$SE^2
      
      # Single complete-case model: no MI pooling required
      if (length(estimates) == 1) {
        
        estimate <- estimates
        standard_error <- sqrt(variances)
        degrees_freedom <- .x[[df_column]]
        
        return(
          tibble(
            m = 1,
            estimate = estimate,
            SE = standard_error,
            df = degrees_freedom,
            conf_low = estimate -
              qt(0.975, degrees_freedom) * standard_error,
            conf_high = estimate +
              qt(0.975, degrees_freedom) * standard_error
          )
        )
      }
      
      complete_data_df <- min(.x[[df_column]], na.rm = TRUE)
      
      pooled <- mice::pool.scalar(
        Q = estimates,
        U = variances,
        n = complete_data_df + 1,
        k = 1
      )
      
      pooled_df <- pooled$df
      pooled_se <- sqrt(pooled$t)
      
      tibble(
        m = pooled$m,
        estimate = pooled$qbar,
        SE = pooled_se,
        df = pooled_df,
        conf_low = pooled$qbar -
          qt(0.975, pooled_df) * pooled_se,
        conf_high = pooled$qbar +
          qt(0.975, pooled_df) * pooled_se,
        FMI = pooled$fmi
      )
      
    }) %>%
    ungroup()
  
  pooled_results
}

#### Pool EMMs ####
pool_dataset_emmeans <- function(emmeans_list) {
  
  pool_mi_estimates(
    results_list = emmeans_list,
    grouping_variables = c("Region", "Season"),
    estimate_column = "emmean"
  )
}


pc1_pooled_emmeans <- map(
  pc1_emmeans,
  ~ list(
    region = pool_dataset_emmeans(.x$region),
    season = pool_mi_estimates(
      results_list = .x$season,
      grouping_variables = c("Season", "Region"),
      estimate_column = "emmean"
    )
  )
)

#### Pool pairwise contrasts ####
pool_dataset_contrasts <- function(
    contrast_list,
    grouping_variables
) {
  
  pool_mi_estimates(
    results_list = contrast_list,
    grouping_variables = grouping_variables,
    estimate_column = "estimate"
  ) %>%
    mutate(
      t_ratio = estimate / SE,
      p_value = 2 * pt(
        abs(t_ratio),
        df = df,
        lower.tail = FALSE
      )
    )
}


adjust_pairwise_pvalues <- function(
    contrast_results,
    adjustment_group
) {
  
  contrast_results %>%
    group_by(across(all_of(adjustment_group))) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    ungroup()
}


pc1_pooled_contrasts <- imap(
  pc1_pairwise,
  function(pairwise_results, dataset_name) {
    
    region_contrasts <- pool_dataset_contrasts(
      pairwise_results$region,
      grouping_variables = c("contrast", "Season")
    ) %>%
      adjust_pairwise_pvalues("Season")
    
    season_contrasts <- pool_dataset_contrasts(
      pairwise_results$season,
      grouping_variables = c("contrast", "Region")
    ) %>%
      adjust_pairwise_pvalues("Region")
    
    list(
      region = region_contrasts,
      season = season_contrasts
    )
  }
)


# Example:
pc1_pooled_contrasts$primary$region
pc1_pooled_contrasts$primary$season

pc1_pooled_contrasts$complete$region

pc1_pooled_contrasts$primary_five$region
pc1_pooled_contrasts$simple$region
pc1_pooled_contrasts$complete$region

pc1_contrast_primary <- bind_rows(
  pc1_pooled_contrasts$primary$region,
  pc1_pooled_contrasts$primary$season
)

# Export supplementary
write_tsv(
  pc1_pooled_contrasts$primary$region %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_primary_region.tsv",
)
write_tsv(
  pc1_pooled_contrasts$primary_five$region %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_primary_five_region.tsv",
)
write_tsv(
  pc1_pooled_contrasts$simple$region %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_simple_region.tsv",
)
write_tsv(
  pc1_pooled_contrasts$complete$region %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_complete_region.tsv",
)



write_tsv(
  pc1_pooled_contrasts$primary$season %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_primary_season.tsv",
)
write_tsv(
  pc1_pooled_contrasts$primary_five$season %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_primary_five_season.tsv",
)
write_tsv(
  pc1_pooled_contrasts$simple$season %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_simple_season.tsv",
)
write_tsv(
  pc1_pooled_contrasts$complete$season %>%
    mutate(
      across(
        c(p_value, p_adj),
        ~ if_else(.x < 0.001, sprintf("%.2e", .x), sprintf("%.3f", .x))
      )
    ) %>%
    select(-significance, -m),
  "supplementary_files/Supplementary_pc1_pooled_contrasts_complete_season.tsv",
)

### Plots ####

#### Significance letters ####
get_significance_letters <- function(
    contrast_data,
    region_levels
) {
  
  p_values <- setNames(
    contrast_data$p_adj,
    gsub("\\s*-\\s*", "-", contrast_data$contrast)
  )
  
  compact_letters <- multcompView::multcompLetters(p_values)
  
  letters_data <- tibble(
    Region = factor(
      names(compact_letters$Letters),
      levels = region_levels
    ),
    raw_letter = compact_letters$Letters
  ) %>%
    arrange(Region)
  
  # Reassign letters according to left-to-right region order
  unique_combinations <- unique(letters_data$raw_letter)
  
  letter_mapping <- setNames(
    letters[seq_along(unique_combinations)],
    unique_combinations
  )
  
  letters_data %>%
    mutate(letter = letter_mapping[raw_letter])
}

#### Plot EMMs ####

plot_pc1_emmeans <- function(
    dataset = c("primary", "primary_five", "simple", "complete"),
    show_ns = FALSE
) {
  
  dataset <- match.arg(dataset)
  
  # Extract corresponding EMMs and pairwise contrasts
  emmeans_data <- pc1_pooled_emmeans[[dataset]]$region
  
  region_contrasts <- pc1_pooled_contrasts[[dataset]]$region
  season_contrasts <- pc1_pooled_contrasts[[dataset]]$season
  
  # Calculate overall y-axis range
  y_range <- diff(
    range(
      c(emmeans_data$conf_low, emmeans_data$conf_high),
      na.rm = TRUE
    )
  )
  
  # Plot 1: Region comparisons within each season
  
  season_annotations <- emmeans_data %>%
    group_by(Region) %>%
    summarise(
      y_pos = max(conf_high, na.rm = TRUE) + y_range * 0.05,
      .groups = "drop"
    ) %>%
    left_join(
      season_contrasts,
      by = "Region"
    )
  
  if (!show_ns) {
    season_annotations <- season_annotations %>%
      filter(significance != "ns")
  }
  
  plot_region <- ggplot(
    emmeans_data,
    aes(
      x = Region,
      y = estimate,
      color = Season,
      group = Season
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = conf_low,
        ymax = conf_high
      ),
      position = position_dodge(width = 0.5),
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_text(
      data = season_annotations,
      aes(
        x = Region,
        y = y_pos,
        label = significance
      ),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4.5,
      fontface = "bold"
    ) +
    labs(
      x = "Region",
      y = "Estimated marginal mean of PC1",
      color = "Season"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(
        color = "gray90",
        linewidth = 0.5
      ),
      panel.grid.minor.y = element_line(
        color = "gray90",
        linewidth = 0.3
      ),
      legend.position = "bottom"
    )
  
  # Plot 2: Region comparisons within each season (faceted)
  
  letters_data <- region_contrasts %>%
    group_by(Season) %>%
    group_modify(
      ~ get_significance_letters(
        contrast_data = .x,
        region_levels = levels(emmeans_data$Region)
      )
    ) %>%
    ungroup()
  
  letter_annotations <- emmeans_data %>%
    group_by(Season, Region) %>%
    summarise(
      y_pos = max(conf_high, na.rm = TRUE) + y_range * 0.06,
      .groups = "drop"
    ) %>%
    left_join(
      letters_data,
      by = c("Season", "Region")
    )
  
  plot_season <- ggplot(
    emmeans_data,
    aes(
      x = Region,
      y = estimate,
      color = Season
    )
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = conf_low,
        ymax = conf_high
      ),
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_text(
      data = letter_annotations,
      aes(
        x = Region,
        y = y_pos,
        label = letter
      ),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4
    ) +
    facet_wrap(~ Season) +
    labs(
      x = "Region",
      y = "Estimated marginal mean of PC1"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_line(
        color = "gray90",
        linewidth = 0.5
      ),
      panel.grid.minor.y = element_line(
        color = "gray90",
        linewidth = 0.3
      ),
      legend.position = "bottom"
    )
  
  # Combine plots
  
  combined_plot <- plot_region + plot_season +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  
  # Return plots and annotations
  return(
    combined_plot
  )
}

pc1_plots <- purrr::map(
  c("primary", "primary_five", "simple", "complete"),
  plot_pc1_emmeans
)

names(pc1_plots) <- c(
  "primary",
  "primary_five",
  "simple",
  "complete"
)

#### Export plots ####
ggsave("supplementary_files/emmeans_plot_primary.png",
       plot = pc1_plots$primary,
       width = 7, height = 4, dpi = 600)
ggsave("supplementary_files/emmeans_plot_primary_five.png",
       plot = pc1_plots$primary_five,
       width = 7, height = 4, dpi = 600)
ggsave("supplementary_files/emmeans_plot_simple.png",
       plot = pc1_plots$simple,
       width = 7, height = 4, dpi = 600)
ggsave("supplementary_files/emmeans_plot_complete.png",
       plot = pc1_plots$complete,
       width = 7, height = 4, dpi = 600)

## Random-Effects Analysis ####

### Fit models ####
primary_random_effect_formulas <- list(
  
  site_only =
    PC1 ~ Season * Region + (1 | Site),
  
  campaign_only =
    PC1 ~ Season * Region + (1 | campaign_id),
  
  site_campaign =
    PC1 ~ Season * Region +
    (1 | Site) +
    (1 | campaign_id)
)


# Fit each random-effects structure across primary imputations
primary_random_effect_models <- map(
  primary_random_effect_formulas,
  ~ map(
    pc1_scores$primary,
    function(dataset) {
      
      lmer(
        formula = .x,
        data = dataset,
        REML = FALSE
      )
    }
  )
)

# Convert to mira objects
primary_random_effect_models <- map(
  primary_random_effect_models,
  mice::as.mira
)

### Compare random-effects structures ####

primary_random_effect_comparison <- purrr::imap_dfr(
  primary_random_effect_models,
  function(model, model_name) {
    
    purrr::map_dfr(
      seq_along(model$analyses),
      function(iteration) {
        
        fit <- model$analyses[[iteration]]
        
        # Extract random-effect variances
        vc <- as.data.frame(VarCorr(fit))
        
        get_variance <- function(group) {
          value <- vc$vcov[vc$grp == group]
          
          if (length(value) == 0) {
            return(NA_real_)
          }
          
          value[1]
        }
        
        # Check convergence
        optimizer_code <- fit@optinfo$conv$opt
        
        convergence_messages <-
          fit@optinfo$conv$lme4$messages
        
        convergence_warnings <-
          fit@optinfo$warnings
        
        all_msgs <- c(
          unlist(convergence_messages),
          unlist(convergence_warnings)
        )
        
        conv_msgs <-
          all_msgs[!grepl("boundary \\(singular\\) fit", all_msgs)]
        
        converged <- 
          (is.null(optimizer_code) || optimizer_code == 0) &&
          (length(conv_msgs) == 0)
        
        # Check singularity
        singular <- lme4::isSingular(
          fit,
          tol = 1e-4
        )
        
        tibble(
          imputation = iteration,
          model = model_name,
          
          # Model fit
          AIC = AIC(fit),
          logLik = as.numeric(logLik(fit)),
          
          # Diagnostics
          converged = converged,
          optimizer_code = if (is.null(optimizer_code)) {
            NA_real_
          } else {
            as.numeric(optimizer_code)
          },
          
          convergence_messages = if (
            is.null(convergence_messages)
          ) {
            NA_character_
          } else {
            paste(convergence_messages, collapse = "; ")
          },
          
          convergence_warnings = if (
            is.null(convergence_warnings)
          ) {
            NA_character_
          } else {
            paste(convergence_warnings, collapse = "; ")
          },
          
          singular = singular,
          
          # Random-effect variance components
          variance_site = get_variance("Site"),
          variance_campaign = get_variance("campaign_id"),
          variance_residual = get_variance("Residual")
        )
      }
    )
  }
)


### Calculate delta AIC ####

primary_random_effect_comparison <-
  primary_random_effect_comparison %>%
  group_by(imputation) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    best_model = delta_AIC == min(delta_AIC)
  ) %>%
  ungroup() %>%
  arrange(imputation, delta_AIC)

# Summarize
primary_random_effect_summary <-
  primary_random_effect_comparison %>%
  group_by(model) %>%
  summarise(
    
    # Model fit
    mean_AIC = mean(AIC),
    sd_AIC = sd(AIC),
    
    mean_delta_AIC = mean(delta_AIC),
    sd_delta_AIC = sd(delta_AIC),
    
    proportion_best_model = mean(best_model),
    
    # Diagnostics
    n_converged = sum(converged),
    n_nonconverged = sum(!converged),
    
    n_singular = sum(singular),
    
    # Random-effect variances
    mean_variance_site = mean(
      variance_site,
      na.rm = TRUE
    ),
    sd_variance_site = sd(
      variance_site,
      na.rm = TRUE
    ),
    
    mean_variance_campaign = mean(
      variance_campaign,
      na.rm = TRUE
    ),
    sd_variance_campaign = sd(
      variance_campaign,
      na.rm = TRUE
    ),
    
    mean_variance_residual = mean(
      variance_residual,
      na.rm = TRUE
    ),
    sd_variance_residual = sd(
      variance_residual,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(mean_delta_AIC)


primary_random_effect_summary

### DHARMa diagnostics: Site-only models ####

set.seed(123)

primary_site_only_dharma <- purrr::map(
  primary_random_effect_models$site_only$analyses,
  function(fit) {
    
    DHARMa::simulateResiduals(
      fittedModel = fit,
      n = 1000,
      refit = FALSE
    )
  }
)

names(primary_site_only_dharma) <- paste0(
  "imputation_",
  seq_along(primary_site_only_dharma)
)

# Run tests
primary_site_only_dharma_tests <- purrr::imap_dfr(
  primary_site_only_dharma,
  function(simulation, iteration) {
    
    uniformity <- DHARMa::testUniformity(
      simulation,
      plot = FALSE
    )
    
    dispersion <- DHARMa::testDispersion(
      simulation,
      plot = FALSE
    )
    
    outliers <- DHARMa::testOutliers(
      simulation,
      plot = FALSE
    )
    
    tibble(
      imputation = iteration,
      
      uniformity_statistic = uniformity$statistic,
      uniformity_p = uniformity$p.value,
      
      dispersion_statistic = dispersion$statistic,
      dispersion_p = dispersion$p.value,
      
      outlier_statistic = outliers$statistic,
      outlier_p = outliers$p.value
    )
  }
)

# Summarize DHARMa diagnostics

primary_site_only_dharma_summary <-
  primary_site_only_dharma_tests %>%
  summarise(
    mean_uniformity_statistic = mean(
      uniformity_statistic,
      na.rm = TRUE
    ),
    mean_uniformity_p = mean(
      uniformity_p,
      na.rm = TRUE
    ),
    mean_dispersion_statistic = mean(
      dispersion_statistic,
      na.rm = TRUE
    ),
    mean_dispersion_p = mean(
      dispersion_p,
      na.rm = TRUE
    ),
    mean_outlier_statistic = mean(
      outlier_statistic,
      na.rm = TRUE
    ),
    mean_outlier_p = mean(
      outlier_p,
      na.rm = TRUE
    ),
    n_passed = sum(
      uniformity_p >= 0.05 &
        dispersion_p >= 0.05 &
        outlier_p >= 0.05,
      na.rm = TRUE
    )
  )

primary_site_only_dharma_summary

### Export supplementary ####
write_tsv(
  primary_random_effect_summary %>%
    select(-sd_AIC, -sd_delta_AIC, -proportion_best_model, -n_nonconverged),
  "supplementary_files/Supplementary_pc1_random_effect_comparison.tsv",
)

write_tsv(
  primary_site_only_dharma_summary,
  "supplementary_files/Supplementary_pc1_primary_dharma.tsv",
)


# ----- PERMANOVA ANALYSES -----
## Analysis Functions ####

# Seasonal PERMANOVA: sample-level observations,
# permutations restricted within Site.
run_season_permanova <- function(dat, variables) {

  X <- dat %>%
    select(all_of(variables)) %>%
    scale()

  distance <- dist(X, method = "euclidean")

  set.seed(123)

  adonis2(
    distance ~ Season,
    data = dat,
    permutations = n_permutations,
    strata = dat$Site,
    by = "terms"
  )
}

# Seasonal PERMDISP
run_season_permdisp <- function(dat, variables) {

  X <- dat %>%
    select(all_of(variables)) %>%
    scale()

  distance <- dist(X, method = "euclidean")

  dispersion <- betadisper(
    distance,
    group = dat$Season
  )

  set.seed(123)

  permutest(
    dispersion,
    permutations = n_permutations,
    strata = dat$Site
  )
}

# Regional PERMANOVA and PERMDISP:
# first calculate one mean profile per Site.
run_region_tests <- function(dat, variables) {

  site_profiles <- dat %>%
    group_by(Site, Region) %>%
    summarise(
      across(all_of(variables), mean, na.rm = TRUE),
      .groups = "drop"
    )

  X <- site_profiles %>%
    select(all_of(variables)) %>%
    scale()

  distance <- dist(X, method = "euclidean")

  set.seed(123)

  permanova <- adonis2(
    distance ~ Region,
    data = site_profiles,
    permutations = n_permutations
  )

  dispersion <- betadisper(
    distance,
    group = site_profiles$Region
  )

  set.seed(123)

  permdisp <- permutest(
    dispersion,
    permutations = n_permutations
  )

  list(
    profiles = site_profiles,
    permanova = permanova,
    permdisp = permdisp
  )
}

# Extract PERMANOVA statistics for one term
extract_permanova <- function(result, term, imputation) {

  result_df <- as.data.frame(result)

  tibble(
    Imputation = imputation,
    R2 = result_df[term, "R2"],
    F_statistic = result_df[term, "F"],
    p_value = result_df[term, "Pr(>F)"]
  )
}

# Extract PERMDISP statistics
extract_permdisp <- function(result, imputation) {

  tab <- result$tab

  tibble(
    Imputation = imputation,
    F_statistic = tab[1, "F"],
    p_value = tab[1, "Pr(>F)"]
  )
}

# Summarize statistics across imputations
summarise_permutation_results <- function(results) {

  results %>%
    select(-Imputation) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          median = ~ median(.x, na.rm = TRUE),
          min = ~ min(.x, na.rm = TRUE),
          max = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    )
}

## Primary Multivariate Analysis ####
# Six parameters (including DO)

### Seasonal PERMANOVA and PERMDISP ####
season_permanova_primary <- map(
  pc_primary_transformed,
  ~ run_season_permanova(.x, pc_vars_primary)
)

season_permdisp_primary <- map(
  pc_primary_transformed,
  ~ run_season_permdisp(.x, pc_vars_primary)
)

season_permanova_primary_results <- bind_rows(
  lapply(seq_along(season_permanova_primary), function(i) {
    extract_permanova(
      season_permanova_primary[[i]],
      term = "Season",
      imputation = i
    )
  })
)

season_permdisp_primary_results <- bind_rows(
  lapply(seq_along(season_permdisp_primary), function(i) {
    extract_permdisp(season_permdisp_primary[[i]], i)
  })
)

### Regional PERMANOVA and PERMDISP ####
region_tests_primary <- map(
  pc_primary_transformed,
  ~ run_region_tests(.x, pc_vars_primary)
)

region_permanova_primary_results <- bind_rows(
  lapply(seq_along(region_tests_primary), function(i) {
    extract_permanova(
      region_tests_primary[[i]]$permanova,
      term = "Model",
      imputation = i
    )
  })
)

region_permdisp_primary_results <- bind_rows(
  lapply(seq_along(region_tests_primary), function(i) {
    extract_permdisp(
      region_tests_primary[[i]]$permdisp,
      i
    )
  })
)

# Summarize
primary_permanova_results <- rbind(
  summarise_permutation_results(season_permanova_primary_results) %>%
    mutate(Test = "Season") %>%
    relocate(Test),
  summarise_permutation_results(region_permanova_primary_results) %>%
    mutate(Test = "Region") %>%
    relocate(Test)
)

primary_permdisp_results <- rbind(
  summarise_permutation_results(season_permdisp_primary_results) %>%
    mutate(Test = "Season") %>%
    relocate(Test),
  summarise_permutation_results(region_permdisp_primary_results) %>%
    mutate(Test = "Region") %>%
    relocate(Test)
)

write_tsv(
  primary_permanova_results,
  "supplementary_files/Supplementary_primary_permanova_results.tsv",
)
write_tsv(
  primary_permdisp_results,
  "supplementary_files/Supplementary_primary_permdisp_results.tsv",
)



# ----- SENSITIVITY ANALYSES -----
## SA1: Simple MICE ####

# Repeat multivariate analyses using simple MICE.
season_permanova_simple <- map(
  pc_simple_transformed,
  ~ run_season_permanova(.x, pc_vars_primary)
)

season_permdisp_simple <- map(
  pc_simple_transformed,
  ~ run_season_permdisp(.x, pc_vars_primary)
)

region_tests_simple <- map(
  pc_simple_transformed,
  ~ run_region_tests(.x, pc_vars_primary)
)

season_permanova_simple_results <- bind_rows(
  lapply(seq_along(season_permanova_simple), function(i) {
    extract_permanova(
      season_permanova_simple[[i]],
      "Season",
      i
    )
  })
)

season_permdisp_simple_results <- bind_rows(
  lapply(seq_along(season_permdisp_simple), function(i) {
    extract_permdisp(season_permdisp_simple[[i]], i)
  })
)

region_permanova_simple_results <- bind_rows(
  lapply(seq_along(region_tests_simple), function(i) {
    extract_permanova(
      region_tests_simple[[i]]$permanova,
      "Model",
      i
    )
  })
)

region_permdisp_simple_results <- bind_rows(
  lapply(seq_along(region_tests_simple), function(i) {
    extract_permdisp(
      region_tests_simple[[i]]$permdisp,
      i
    )
  })
)

# Summaries
print(summarise_permutation_results(
  season_permanova_simple_results
))
print(summarise_permutation_results(
  season_permdisp_simple_results
))
print(summarise_permutation_results(
  region_permanova_simple_results
))
print(summarise_permutation_results(
  region_permdisp_simple_results
))

## SA2: Complete-Case ####

# Seasonal tests
season_permanova_complete <- run_season_permanova(
  pc_complete_transformed,
  pc_vars_primary
)

season_permdisp_complete <- run_season_permdisp(
  pc_complete_transformed,
  pc_vars_primary
)

# Regional tests
region_tests_complete <- run_region_tests(
  pc_complete_transformed,
  pc_vars_primary
)

## SA3: Excluding DO ####

# This analysis evaluates how the BDL-substituted DO values
# influence the multivariate results. It is not the primary
# analysis because DO censoring is strongly region-dependent.

season_permanova_five <- map(
  pc_primary_transformed,
  ~ run_season_permanova(.x, pc_vars_five)
)

season_permdisp_five <- map(
  pc_primary_transformed,
  ~ run_season_permdisp(.x, pc_vars_five)
)

region_tests_five <- map(
  pc_primary_transformed,
  ~ run_region_tests(.x, pc_vars_five)
)

season_permanova_five_results <- bind_rows(
  lapply(seq_along(season_permanova_five), function(i) {
    extract_permanova(
      season_permanova_five[[i]],
      "Season",
      i
    )
  })
)

season_permdisp_five_results <- bind_rows(
  lapply(seq_along(season_permdisp_five), function(i) {
    extract_permdisp(season_permdisp_five[[i]], i)
  })
)

region_permanova_five_results <- bind_rows(
  lapply(seq_along(region_tests_five), function(i) {
    extract_permanova(
      region_tests_five[[i]]$permanova,
      "Model",
      i
    )
  })
)

region_permdisp_five_results <- bind_rows(
  lapply(seq_along(region_tests_five), function(i) {
    extract_permdisp(region_tests_five[[i]]$permdisp, i)
  })
)

print(summarise_permutation_results(
  season_permanova_five_results
))
print(summarise_permutation_results(
  season_permdisp_five_results
))
print(summarise_permutation_results(
  region_permanova_five_results
))
print(summarise_permutation_results(
  region_permdisp_five_results
))


## Compare Primary And Sensitivity Results ####

# Compact summaries for reporting
permanova_comparison <- bind_rows(
  summarise_permutation_results(
    season_permanova_primary_results
  ) %>% mutate(Analysis = "Primary: cluster-aware MICE, all parameters",
               Test = "Season"),
  summarise_permutation_results(
    region_permanova_primary_results
  ) %>% mutate(Analysis = "Primary: cluster-aware MICE, all parameters",
               Test = "Region"),
  summarise_permutation_results(
    season_permanova_simple_results
  ) %>% mutate(Analysis = "Sensitivity: simple MICE, all parameters",
               Test = "Season"),
  summarise_permutation_results(
    region_permanova_simple_results
  ) %>% mutate(Analysis = "Sensitivity: simple MICE, all parameters",
               Test = "Region"),
  summarise_permutation_results(
    season_permanova_five_results
  ) %>% mutate(Analysis = "Sensitivity: cluster-aware MICE, w/o DO",
               Test = "Season"),
  summarise_permutation_results(
    region_permanova_five_results
  ) %>% mutate(Analysis = "Sensitivity: cluster-aware MICE, w/o DO",
               Test = "Region")
) %>%
  relocate(Analysis, Test)

print(season_permanova_complete)
print(season_permdisp_complete)
print(region_tests_complete$permanova)
print(region_tests_complete$permdisp)

print(permanova_comparison)

write_tsv(
  permanova_comparison,
  "supplementary_files/Supplementary_permanova_comparison.tsv",
)


