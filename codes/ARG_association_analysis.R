library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(readr)
library(ggplot2)
library(ComplexUpset)
library(forcats)
library(lubridate)

library(purrr)
library(lme4)
library(DHARMa)
library(broom.mixed)

library(flextable)
library(officer)

# ==== ARG Association ====
## Load Data ####
arg_cols <- c(
  "blaCTX-M", "blaSHV", "blaTEM",
  "blaNDM", "blaOXA-1", "tetA",
  "sul1", "sul2"
)

# Keep isolates with uidA confirmation and complete ARG results
arg_analysis_data <- read.table("ARPHI-isolate-database.tsv", sep = "\t", header = TRUE) %>%
  rename(
    `blaCTX-M` = blaCTX.M,
    `blaOXA-1` = blaOXA.1
  ) %>%
  filter(
    !is.na(pcr_confirmation),
    if_all(all_of(arg_cols), ~ !is.na(.x))
  ) %>%
  select(campaign_id, water_sample_id, isolate_id, Season, Site, Region, all_of(arg_cols))

# Remove rare genes: blaNDM and blaOXA-1
arg_cols <- c(
  "blaCTX-M", "blaSHV", "blaTEM",
  "tetA", "sul1", "sul2"
)

# Convert data to binary
arg_binary <- arg_analysis_data %>%
  select(-blaNDM, -"blaOXA-1") %>%
  mutate(
    across(
      all_of(arg_cols),
      ~ case_when(
        .x == "positive" ~ 1L,
        .x == "negative" ~ 0L,
        TRUE ~ NA_integer_
      )
    ),
    water_sample_id = factor(water_sample_id)
  )

# Check how many isolates come from each sample
sample_sizes <- arg_analysis_data %>%
  count(water_sample_id, name = "n_isolates")

# Check whether each gene pair has enough variation
arg_variation <- arg_binary %>%
  summarise(
    across(
      all_of(arg_cols),
      list(
        positive = ~ sum(.x == 1),
        negative = ~ sum(.x == 0)
      )
    )
  )

# Generate all unique gene pairs
gene_pairs <- combn(arg_cols, 2, simplify = FALSE)

## UpSet plot ####
gene_order <- c(
  "tetA",
  "blaTEM",
  "sul2",
  "sul1",
  "blaSHV",
  "blaCTX-M",
  "blaNDM",
  "blaOXA-1"
)

df_upset <- arg_analysis_data %>%
  select(7:last_col()) %>%
  mutate(
    across(
      everything(),
      ~ case_when(
        .x == "positive" ~ TRUE,
        .x == "negative" ~ FALSE,
        TRUE ~ NA_integer_
      )
    )
  ) %>%
  filter(
    if_any(everything(), ~ .x)
  )

gene_order <- gene_order[gene_order %in% names(df_upset)]
max_set <- max(colSums(df_upset, na.rm = TRUE))

p_upset <- ComplexUpset::upset(
  df_upset %>%
    select(all_of(gene_order)),
  
  intersect = gene_order,
  name = "Antibiotic Resistance Genes",
  sort_intersections_by = "cardinality",
  sort_intersections = "descending",
  sort_sets = FALSE,
  min_size = 1,
  
  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = TRUE,
        text = list(
          size = 3,
          angle = 90,
          vjust = 0.5,
          hjust = -0.5
        )
      ) +
      labs(
        y = "Number of isolates"
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  ),
  
  set_sizes =
    upset_set_size(
      geom = geom_bar(
        stat = "count"
      )
    ) +
    geom_text(
      stat = "count",
      aes(label = after_stat(count)),
      hjust = 1.1,
      size = 3
    ) +
    expand_limits(y = max_set * 1.15) +
    labs(
      x = "Set size"
    ),
  
  matrix = intersection_matrix(
    geom = geom_point(
      size = 2.5
    ),
    segment = geom_segment(
      linewidth = 0.5
    )
  )
)

ggsave("supplementary_files/UpSet_plot.png",
       plot = p_upset,
       width = 10, height = 5, units = "in", dpi = 600)

# Summarize
# Extract exact ARG combinations and their frequencies
intersection_summary <- df_upset %>%
  mutate(
    intersection = apply(
      select(., all_of(gene_order)),
      1,
      function(x) {
        genes <- gene_order[as.logical(x)]
        paste(genes, collapse = " + ")
      }
    )
  ) %>%
  count(intersection, sort = TRUE)

head(intersection_summary, 10)

#### Mixed effects analysis ####
# Candidate random-effects structures
model_formulas <- list(
  M1_sample = gene_a ~ gene_b +
    (1 | water_sample_id),
  
  M2_sample_site = gene_a ~ gene_b +
    (1 | water_sample_id) +
    (1 | Site),
  
  M3_sample_campaign = gene_a ~ gene_b +
    (1 | water_sample_id) +
    (1 | campaign_id),
  
  M4_sample_site_campaign = gene_a ~ gene_b +
    (1 | water_sample_id) +
    (1 | Site) +
    (1 | campaign_id)
)

# Convert variance components into a readable string
get_variance_components <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  
  paste(
    paste0(
      vc$grp, "=", signif(vc$vcov, 4)
    ),
    collapse = "; "
  )
}

# Extract convergence and singularity diagnostics
get_model_diagnostics <- function(model) {
  opt_code <- model@optinfo$conv$opt
  conv_messages <- model@optinfo$conv$lme4$messages
  
  if (is.null(opt_code)) {
    opt_code <- NA_integer_
  }
  
  if (is.null(conv_messages) || length(conv_messages) == 0) {
    conv_text <- "none"
  } else {
    conv_text <- paste(conv_messages, collapse = " | ")
  }
  
  paste0(
    "singular=", isSingular(model, tol = 1e-4),
    "; optimizer_code=", paste(opt_code, collapse = ","),
    "; convergence_messages=", conv_text
  )
}

# Fit one model for one gene pair
fit_one_pair <- function(pair, formula, model_name,
                         data = arg_binary, run_dharma = TRUE) {
  
  gene_a_name <- pair[1]
  gene_b_name <- pair[2]
  
  # Retain complete observations for the current pair and
  # grouping variables required by the candidate model.
  required_groups <- switch(
    model_name,
    M1_sample = "water_sample_id",
    M2_sample_site = c("water_sample_id", "Site"),
    M3_sample_campaign = c("water_sample_id", "campaign_id"),
    M4_sample_site_campaign =
      c("water_sample_id", "Site", "campaign_id")
  )
  
  df <- data %>%
    transmute(
      across(all_of(required_groups), ~ factor(.x)),
      gene_a = .data[[gene_a_name]],
      gene_b = .data[[gene_b_name]]
    ) %>%
    filter(complete.cases(.))
  
  # A standard result row for failed or non-estimable models
  empty_result <- function(note) {
    tibble(
      gene_a = gene_a_name,
      gene_b = gene_b_name,
      n = nrow(df),
      n_samples = if ("water_sample_id" %in% names(df)) {
        n_distinct(df$water_sample_id)
      } else {
        NA_integer_
      },
      AIC = NA_real_,
      OR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      variance_components = NA_character_,
      model_diagnostics = NA_character_,
      p_value = NA_real_,
      note = note
    )
  }
  
  # Check that both binary variables vary
  if (nrow(df) == 0) {
    return(empty_result("No complete observations"))
  }
  
  if (n_distinct(df$gene_a) < 2) {
    return(empty_result("Outcome has no variation"))
  }
  
  if (n_distinct(df$gene_b) < 2) {
    return(empty_result("Predictor has no variation"))
  }
  
  # Fit model while retaining warnings for the diagnostics
  fit_warnings <- character()
  
  model <- tryCatch(
    withCallingHandlers(
      glmer(
        formula,
        data = df,
        family = binomial(link = "logit"),
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 200000)
        )
      ),
      warning = function(w) {
        fit_warnings <<- c(
          fit_warnings,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(model, "error")) {
    return(empty_result(
      paste("Model error:", conditionMessage(model))
    ))
  }
  
  # ---- DHARMa diagnostics ----
  dharma_uniformity_statistic <- NA_real_
  dharma_uniformity_p <- NA_real_
  
  dharma_dispersion_statistic <- NA_real_
  dharma_dispersion_p <- NA_real_
  
  dharma_outlier_statistic <- NA_real_
  dharma_outlier_p <- NA_real_
  
  dharma_note <- "DHARMa diagnostics not run"
  
  if (run_dharma) {
    
    set.seed(123)
    
    sim_res <- DHARMa::simulateResiduals(
      fittedModel = model,
      n = 1000,
      plot = FALSE
    )
    
    # Save the DHARMa object for later exploration
    pair_key <- paste(gene_a_name, gene_b_name, sep = "__")
    dharma_objects[[pair_key]] <<- sim_res
    
    uniformity_test <- DHARMa::testUniformity(sim_res)
    dispersion_test <- DHARMa::testDispersion(sim_res)
    outlier_test <- DHARMa::testOutliers(sim_res)
    
    # Extract test statistics and p-values
    get_stat <- function(test) {
      if (is.null(test$statistic) || length(test$statistic) == 0) {
        return(NA_real_)
      }
      as.numeric(test$statistic[1])
    }
    
    get_p <- function(test) {
      if (is.null(test$p.value) || length(test$p.value) == 0) {
        return(NA_real_)
      }
      as.numeric(test$p.value[1])
    }
    
    dharma_uniformity_statistic <- get_stat(uniformity_test)
    dharma_uniformity_p <- get_p(uniformity_test)
    
    dharma_dispersion_statistic <- get_stat(dispersion_test)
    dharma_dispersion_p <- get_p(dispersion_test)
    
    dharma_outlier_statistic <- get_stat(outlier_test)
    dharma_outlier_p <- get_p(outlier_test)
    
    # Generate a concise diagnostic summary
    significant_tests <- c(
      if (!is.na(dharma_uniformity_p) &&
          dharma_uniformity_p < 0.05) "Uniformity test",
      
      if (!is.na(dharma_dispersion_p) &&
          dharma_dispersion_p < 0.05) "Dispersion test",
      
      if (!is.na(dharma_outlier_p) &&
          dharma_outlier_p < 0.05) "Outlier test"
    )
    
    unavailable_tests <- c(
      if (is.na(dharma_uniformity_p)) "Uniformity test",
      if (is.na(dharma_dispersion_p)) "Dispersion test",
      if (is.na(dharma_outlier_p)) "Outlier test"
    )
    
    dharma_note <- case_when(
      length(significant_tests) > 0 ~
        paste("Significant:", paste(significant_tests, collapse = ", ")),
      
      length(unavailable_tests) > 0 ~
        paste("Unavailable:", paste(unavailable_tests, collapse = ", ")),
      
      TRUE ~ "No significant deviations"
    )
  }
  
  # Extract fixed-effect estimate and Wald 95% CI
  coefs <- summary(model)$coefficients
  
  if (!"gene_b" %in% rownames(coefs)) {
    return(empty_result("gene_b coefficient unavailable"))
  }
  
  beta <- coefs["gene_b", "Estimate"]
  se <- coefs["gene_b", "Std. Error"]
  p <- coefs["gene_b", "Pr(>|z|)"]
  
  # Combine diagnostics and warnings
  diagnostics <- get_model_diagnostics(model)
  
  if (length(fit_warnings) > 0) {
    diagnostics <- paste0(
      diagnostics,
      "; warnings=",
      paste(unique(fit_warnings), collapse = " | ")
    )
  }
  
  # Flag noteworthy model issues
  notes <- character()
  
  if (isSingular(model, tol = 1e-4)) {
    notes <- c(notes, "Singular fit")
  }
  
  if (!is.null(model@optinfo$conv$opt) &&
      model@optinfo$conv$opt != 0) {
    notes <- c(notes, "Optimizer convergence issue")
  }
  
  if (length(model@optinfo$conv$lme4$messages) > 0) {
    notes <- c(notes, "Convergence message;")
  }
  
  if (length(fit_warnings) > 0) {
    notes <- c(notes, "Model warning;")
  }
  
  tibble(
    gene_a = gene_a_name,
    gene_b = gene_b_name,
    n = nobs(model),
    n_samples = n_distinct(df$water_sample_id),
    AIC = AIC(model),
    OR = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    variance_components = get_variance_components(model),
    model_diagnostics = diagnostics,
    model_status = if (length(notes) == 0) {
      "No convergence issues"
    } else {
      paste(unique(notes), collapse = "; ")
    },
    p_value = p,
    dharma_uniformity_statistic = dharma_uniformity_statistic,
    dharma_uniformity_p = dharma_uniformity_p,
    
    dharma_dispersion_statistic = dharma_dispersion_statistic,
    dharma_dispersion_p = dharma_dispersion_p,
    
    dharma_outlier_statistic = dharma_outlier_statistic,
    dharma_outlier_p = dharma_outlier_p,
    
    dharma_note = dharma_note
  )
}

# Fit all 15 gene pairs for one candidate model

fit_model_table <- function(model_name, formula,
                            data = arg_binary,
                            run_dharma = TRUE) {
  
  message("fit_model_table: run_dharma = ", run_dharma)
  
  map_dfr(
    gene_pairs,
    function(pair) {
      fit_one_pair(
        pair = pair,
        formula = formula,
        model_name = model_name,
        data = data,
        run_dharma = run_dharma
      )
    }
  ) %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = "BH"),
      .after = p_value
    ) %>%
    arrange(p_adjusted)
}

#### Run for all models ####
## Generate results for M1
# Initiate DHARMa list
dharma_objects <- list()
results_M1_sample <- fit_model_table(
  "M1_sample",
  model_formulas$M1_sample,
  run_dharma = TRUE
)

## Explore a specific flagged pair
sim_res <- dharma_objects[["blaSHV__tetA"]]
DHARMa::testOutliers(sim_res, margin = "both")

results_M2_sample_site <- fit_model_table(
  "M2_sample_site",
  model_formulas$M2_sample_site,
  run_dharma = FALSE
)

results_M3_sample_campaign <- fit_model_table(
  "M3_sample_campaign",
  model_formulas$M3_sample_campaign,
  run_dharma = FALSE
)

results_M4_sample_site_campaign <- fit_model_table(
  "M4_sample_site_campaign",
  model_formulas$M4_sample_site_campaign,
  run_dharma = FALSE
)


#### Compare AIC ####
aic_comparison <- bind_rows(
  results_M1_sample %>%
    mutate(model = "M1_sample"),
  results_M2_sample_site %>%
    mutate(model = "M2_sample_site"),
  results_M3_sample_campaign %>%
    mutate(model = "M3_sample_campaign"),
  results_M4_sample_site_campaign %>%
    mutate(model = "M4_sample_site_campaign")
) %>%
  select(gene_a, gene_b, model, n, AIC) %>%
  group_by(gene_a, gene_b) %>%
  mutate(
    delta_AIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  ungroup()

aic_plot_data <- aic_comparison %>%
  mutate(
    gene_pair = paste(gene_a, gene_b, sep = " – "),
    
    evidence = case_when(
      delta_AIC < 2 ~ "< 2: Indistinguishable",
      delta_AIC <= 4 ~ "2–4: Moderate evidence",
      delta_AIC > 4 ~ "> 4: Strong evidence",
      TRUE ~ NA_character_
    ),
    
    model = recode(
      model,
      M1_sample = "M1",
      M2_sample_site = "M2",
      M3_sample_campaign = "M3",
      M4_sample_site_campaign = "M4"
    )
  )

# Order gene pairs according to M1 delta AIC
aic_plot_data <- aic_plot_data %>%
  mutate(
    gene_pair = factor(
      gene_pair,
      levels = levels(plot_data_M1$gene_pair)
    ),
    model = factor(
      model,
      levels = c(
        "M1",
        "M2",
        "M3",
        "M4"
      )
    ),
    evidence = factor(
      evidence,
      levels = c(
        "< 2: Indistinguishable",
        "2–4: Moderate evidence",
        "> 4: Strong evidence"
      )
    )
  )

#### Create Supplementary Data ####
format_export_results <- function(data, filename,
                                  exclude_cols = c("n", "n_samples",
                                                   "model_diagnostics")) {
  
  formatted_results <- data %>%
    select(-any_of(exclude_cols)) %>%
    mutate(
      dharma_outlier_statistic = as.integer(dharma_outlier_statistic),
      across(
        where(is.double),
        ~ case_when(
          is.na(.x) ~ NA_character_,
          .x == 0 ~ "0.000",
          abs(.x) < 0.001 ~ sprintf("%.2e", .x),
          TRUE ~ sprintf("%.3f", .x)
        )
      )
    )
  
  write_tsv(
    formatted_results,
    filename
  )
  
  return(formatted_results)
}

format_export_results(
  results_M1_sample,
  "supplementary_files/Supplementary_ARG_assoc_M1.tsv"
)

format_export_results(
  results_M2_sample_site,
  "supplementary_files/Supplementary_ARG_assoc_M2.tsv"
)

format_export_results(
  results_M3_sample_campaign,
  "supplementary_files/Supplementary_ARG_assoc_M3.tsv"
)

format_export_results(
  results_M4_sample_site_campaign,
  "supplementary_files/Supplementary_ARG_assoc_M4.tsv"
)

#### PLOT ####
plot_data_M1 <- results_M1_sample %>%
  filter(!is.na(OR), !is.na(CI_low), !is.na(CI_high)) %>%
  mutate(
    gene_pair = paste(gene_a, gene_b, sep = " – "),
    significance = p_adjusted < 0.05,
    gene_pair = fct_reorder(gene_pair, OR)
  )

plot_data_M2 <- results_M2_sample_site %>%
  filter(!is.na(OR), !is.na(CI_low), !is.na(CI_high)) %>%
  mutate(
    gene_pair = paste(gene_a, gene_b, sep = " – "),
    significance = p_adjusted < 0.05,
    gene_pair = factor(
      gene_pair,
      levels = levels(plot_data_M1$gene_pair)
    )
  )

plot_data_M3 <- results_M3_sample_campaign %>%
  filter(!is.na(OR), !is.na(CI_low), !is.na(CI_high)) %>%
  mutate(
    gene_pair = paste(gene_a, gene_b, sep = " – "),
    significance = p_adjusted < 0.05,
    gene_pair = factor(
      gene_pair,
      levels = levels(plot_data_M1$gene_pair)
    )
  )

plot_data_M4 <- results_M4_sample_site_campaign %>%
  filter(!is.na(OR), !is.na(CI_low), !is.na(CI_high)) %>%
  mutate(
    gene_pair = paste(gene_a, gene_b, sep = " – "),
    significance = p_adjusted < 0.05,
    gene_pair = factor(
      gene_pair,
      levels = levels(plot_data_M1$gene_pair)
    )
  )

# AIC Comparison Plot
AIC_compare_plot <- ggplot(
  aic_plot_data,
  aes(
    x = delta_AIC,
    y = gene_pair,
    fill = evidence,
    shape = model,
    group = model
  )
) +
  geom_hline(
    yintercept = seq(
      0.5,
      length(levels(aic_plot_data$gene_pair)) + 0.5,
      by = 1
    ),
    color = "grey88",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(2, 4),
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.4
  ) +
  geom_point(
    size = 3,
    stroke = 0.9,
    color = "black",
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "< 2: Indistinguishable" = "green",
      "2–4: Moderate evidence" = "orange",
      "> 4: Strong evidence" = "red"
    )
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24)
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      ceiling(max(aic_plot_data$delta_AIC, na.rm = TRUE))
    ),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  labs(
    x = expression(Delta * "AIC"),
    y = NULL,
    fill = "Evidence against model",
    shape = "Model"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        fill = c("green", "orange", "red"),
        color = "black",
        alpha = 1
      )
    ),
    shape = guide_legend(
      override.aes = list(
        fill = "grey70",
        color = "black",
        alpha = 0.8
      )
    )
  )+
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

# Model Plots
main_plot_m1 <- ggplot(plot_data_M1, aes(x = OR, y = gene_pair)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(aes(shape = significance), size = 2.8) +
  scale_x_log10() +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "adjusted p < 0.05", `FALSE` = "adjusted p ≥ 0.05")
  ) +
  labs(
    x = "Odds ratio (log scale; 95% CI)",
    y = NULL,
    shape = "Adjusted significance"
  ) +
  theme_classic(base_size = 11)

plot_m1 <- ggplot(plot_data_M1, aes(x = OR, y = gene_pair)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(aes(shape = significance), size = 2.8) +
  scale_x_log10() +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "adjusted p < 0.05", `FALSE` = "adjusted p ≥ 0.05")
  ) +
  labs(
    x = "Odds ratio (log scale; 95% CI)",
    y = NULL,
    shape = "Adjusted significance",
    title = "Model M1"
  ) +
  theme_classic(base_size = 11)

plot_m2 <- ggplot(plot_data_M2, aes(x = OR, y = gene_pair)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(aes(shape = significance), size = 2.8) +
  scale_x_log10() +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "adjusted p < 0.05", `FALSE` = "adjusted p ≥ 0.05")
  ) +
  labs(
    x = "Odds ratio (log scale; 95% CI)",
    y = NULL,
    shape = "Adjusted significance",
    title = "Model M2"
  ) +
  theme_classic(base_size = 11)

plot_m3 <- ggplot(plot_data_M3, aes(x = OR, y = gene_pair)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(aes(shape = significance), size = 2.8) +
  scale_x_log10() +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "adjusted p < 0.05", `FALSE` = "adjusted p ≥ 0.05")
  ) +
  labs(
    x = "Odds ratio (log scale; 95% CI)",
    y = NULL,
    shape = "Adjusted significance",
    title = "Model M3"
  ) +
  theme_classic(base_size = 11)

plot_m4 <- ggplot(plot_data_M4, aes(x = OR, y = gene_pair)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(aes(shape = significance), size = 2.8) +
  scale_x_log10() +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "adjusted p < 0.05", `FALSE` = "adjusted p ≥ 0.05")
  ) +
  labs(
    x = "Odds ratio (log scale; 95% CI)",
    y = NULL,
    shape = "Adjusted significance",
    title = "Model M4"
  ) +
  theme_classic(base_size = 11)

# Save plots
ggsave("supplementary_files/aic_compare_plot.png", plot = AIC_compare_plot, width = 8, height = 4)

ggsave("Figure_plot_m1.png", plot = main_plot_m1, width = 8, height = 4)

ggsave("supplementary_files/plot_m1.png", plot = plot_m1, width = 6, height = 4)
ggsave("supplementary_files/plot_m2.png", plot = plot_m2, width = 6, height = 4)
ggsave("supplementary_files/plot_m3.png", plot = plot_m3, width = 6, height = 4)
ggsave("supplementary_files/plot_m4.png", plot = plot_m4, width = 6, height = 4)

