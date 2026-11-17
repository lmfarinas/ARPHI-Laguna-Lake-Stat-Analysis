library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(lubridate)
library(mice)
library(pan)
library(glmmTMB)
library(lme4)
library(broom.mixed)
library(DHARMa)

set.seed(2026)

output_dir <- "supplementary_files"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_imputations <- 20
n_iterations <- 20

# Physicochemical variables used in the primary PCA.
pc_vars_primary <- c(
  "log10_Ammonia", "log10_BOD", "DO",
  "log10_FC", "pH", "log10_Phosphate"
)

# Six ARGs retained for the association analysis.
arg_genes <- c(
  "blaCTX-M", "blaSHV", "blaTEM",
  "tetA", "sul1", "sul2"
)

# Environmental predictors tested against ARG prevalence.
env_predictors <- c(
  "PC1",
  "log10_Ammonia",
  "log10_BOD",
  "DO",
  "log10_FC",
  "pH",
  "log10_Phosphate",
  "log10_Nitrate"
)

# Random-effects structures compared during model selection.
random_effect_structures <- c(
  "none", "site", "campaign", "site_campaign"
)


# ----- Generate the primary MICE dataset ------------------------------------

## Load and prepare physicochemical data ####

pc_vars_all <- c(
  "Nitrate", "Ammonia", "BOD", "DO",
  "FC", "pH", "Phosphate"
)

sample_metadata_complete <- read.delim(
  "sample_metadata.tsv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    sample_collection_date = as.Date(sample_collection_date),
    Year = lubridate::year(sample_collection_date)
  )

pc_all_long <- sample_metadata_complete %>%
  select(
    campaign_id, water_sample_id, Site, Region, Season,
    sample_collection_date, Year,
    all_of(pc_vars_all)
  ) %>%
  mutate(across(all_of(pc_vars_all), as.character)) %>%
  pivot_longer(
    cols = all_of(pc_vars_all),
    names_to = "Parameter",
    values_to = "Raw_Value"
  ) %>%
  mutate(
    Raw_Value = str_trim(Raw_Value),
    Status = case_when(
      is.na(Raw_Value) | Raw_Value == "***" ~ "Missing",
      str_starts(Raw_Value, "<") ~ "BRL",
      TRUE ~ "Quantified"
    ),
    Reporting_Limit = if_else(
      Status == "BRL",
      as.numeric(str_remove(Raw_Value, "^<")),
      NA_real_
    ),
    Value = if_else(
      Status == "Quantified",
      as.numeric(Raw_Value),
      NA_real_
    )
  )

# Substitute values below the reporting limit (BRL) with reporting limit / sqrt(2).
pc_limit_sqrt2_wide <- pc_all_long %>%
  mutate(
    Value_processed = case_when(
      Status == "BRL" ~ Reporting_Limit / sqrt(2),
      Status == "Quantified" ~ Value,
      TRUE ~ NA_real_
    )
  ) %>%
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


## Prepare data for two-level MICE ####

# Nitrate is retained for the environmental association models but excluded
# from the MICE procedure and primary PCA.
pc_vars_log <- c("Ammonia", "BOD", "FC", "Phosphate")

mice_data <- pc_limit_sqrt2_wide %>%
  select(-Nitrate) %>%
  mutate(
    Site_cluster = as.integer(factor(Site)),
    Season_WET = as.integer(Season == "WET"),
    Year_c = Year - mean(Year, na.rm = TRUE),
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

# Standardize continuous variables for numerical stability during imputation.
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


## Run primary site-clustered MICE ####

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

# Site random-intercept grouping variable.
pred_cluster[
  c("log10_Ammonia", "log10_BOD"),
  "Site_cluster"
] <- -2

# Fixed-effect predictors.
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

imp_primary <- mice(
  mice_data_scaled,
  m = n_imputations,
  maxit = n_iterations,
  method = meth_cluster,
  predictorMatrix = pred_cluster,
  printFlag = TRUE
)


## Back-transform imputations and restore sample metadata ####

completed_primary <- complete(
  imp_primary,
  action = "all"
)

pc_analysis_primary <- lapply(
  completed_primary,
  function(dat) {

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
      ) %>%
      bind_cols(
        pc_limit_sqrt2_wide %>%
          select(
            campaign_id, water_sample_id, Site, Region, Season,
            sample_collection_date, Year, Nitrate
          )
      )
  }
)


# ----- Prepare PCA scores and ARG sample-level data -------------------------

## Transform physicochemical variables and extract PC1 ####

pc_primary_transformed <- lapply(
  pc_analysis_primary,
  function(dat) {
    dat %>%
      mutate(
        log10_Ammonia = log10(Ammonia),
        log10_BOD = log10(BOD),
        log10_FC = log10(FC),
        log10_Phosphate = log10(Phosphate),
        log10_Nitrate = log10(Nitrate)
      )
  }
)

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

region_levels <- c(
  "East Lake", "West Lake",
  "Quadrat 1", "Quadrat 2",
  "Quadrat 3", "Quadrat 4"
)

pc1_scores <- map2(
  pca_primary,
  pc_primary_transformed,
  function(pca_fit, transformed_data) {
    as.data.frame(pca_fit$x) %>%
      bind_cols(
        transformed_data %>%
          select(
            water_sample_id, campaign_id,
            Site, Region, Season
          )
      ) %>%
      transmute(
        water_sample_id,
        campaign_id,
        Site = factor(Site),
        Region = factor(Region, levels = region_levels),
        Season = factor(Season),
        PC1
      )
  }
)


## Load and prepare ARG data ####

arg_analysis_data <- read.table(
  "ARPHI-isolate-database.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  filter(
    !is.na(pcr_confirmation),
    if_all(
      all_of(c(
        "blaCTX-M", "blaSHV", "blaTEM",
        "blaNDM", "blaOXA-1", "tetA", "sul1", "sul2"
      )),
      ~ !is.na(.x)
    )
  ) %>%
  select(
    campaign_id, water_sample_id, isolate_id,
    Season, Site, Region,
    all_of(c(
      "blaCTX-M", "blaSHV", "blaTEM",
      "blaNDM", "blaOXA-1", "tetA", "sul1", "sul2"
    ))
  )

# Exclude rare genes blaNDM and blaOXA-1 from the association analysis.
arg_binary <- arg_analysis_data %>%
  select(-blaNDM, -`blaOXA-1`) %>%
  mutate(
    across(
      all_of(arg_genes),
      ~ case_when(
        .x == "positive" ~ 1L,
        .x == "negative" ~ 0L,
        TRUE ~ NA_integer_
      )
    ),
    water_sample_id = factor(water_sample_id)
  )

# Aggregate isolate-level ARG results to the water-sample level.
arg_sample_summary <- arg_binary %>%
  group_by(water_sample_id) %>%
  summarise(
    Site = first(Site),
    Region = first(Region),
    Season = first(Season),
    campaign_id = first(campaign_id),
    n_isolates = n(),
    across(
      all_of(arg_genes),
      list(
        positive = ~ sum(.x == 1, na.rm = TRUE),
        tested = ~ sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )


## Assemble ARG–PC analysis datasets ####

prepare_arg_pc_data <- function(iteration) {

  pc1_dat <- pc1_scores[[iteration]] %>%
    select(water_sample_id, PC1)

  pc_dat <- pc_primary_transformed[[iteration]] %>%
    select(
      water_sample_id,
      log10_Ammonia, log10_BOD, DO,
      log10_FC, pH, log10_Phosphate, log10_Nitrate
    )

  dat <- arg_sample_summary %>%
    inner_join(pc1_dat, by = "water_sample_id") %>%
    inner_join(pc_dat, by = "water_sample_id") %>%
    mutate(
      Site = factor(Site),
      Region = factor(Region, levels = region_levels),
      Season = factor(Season, levels = c("DRY", "WET")),
      campaign_id = factor(campaign_id)
    ) %>%
    mutate(
      across(
        all_of(env_predictors),
        ~ as.numeric(scale(.x))
      )
    )

  dat
}

arg_pc_data <- map(
  seq_along(pc1_scores),
  prepare_arg_pc_data
)



# ----PRIMARY ANALYSIS: BETA-BINOMIAL GLMMs----

## Fit beta-binomial models across imputations ####

fit_arg_betabinomial <- function(
    gene,
    predictor,
    random_effects = random_effect_structures,
    data_list = arg_pc_data
) {

  random_effects <- match.arg(random_effects)

  response <- paste0(
    "cbind(`", gene, "_positive`, ",
    "`", gene, "_tested` - `", gene, "_positive`)"
  )

  random_formula <- switch(
    random_effects,
    none = "",
    site = "+ (1 | Site)",
    campaign = "+ (1 | campaign_id)",
    site_campaign = "+ (1 | Site) + (1 | campaign_id)"
  )

  formula <- as.formula(
    paste(response, "~", predictor, random_formula)
  )

  fits <- map(data_list, function(dat) {

    dat <- dat %>%
      filter(
        .data[[paste0(gene, "_tested")]] > 0,
        is.finite(.data[[predictor]]),
        !is.na(Site),
        !is.na(campaign_id)
      )

    glmmTMB(
      formula,
      data = dat,
      family = betabinomial(link = "logit"),
      control = glmmTMBControl(
        optCtrl = list(
          iter.max = 10000,
          eval.max = 10000
        )
      )
    )
  })

  pooled_fit <- mice::pool(
    mice::as.mira(fits),
    dfcom = Inf
  )

  results <- summary(
    pooled_fit,
    conf.int = TRUE
  ) %>%
    mutate(
      OR = exp(estimate),
      OR_low = exp(`2.5 %`),
      OR_high = exp(`97.5 %`)
    )

  list(
    formula = formula,
    random_effects = random_effects,
    fits = fits,
    pooled = pooled_fit,
    results = results
  )
}


## Run all beta-binomial models ####

run_all_arg_betabinomial <- function(
    predictors = env_predictors,
    genes = arg_genes,
    random_effects = random_effect_structures
) {

  crossing(
    Gene = genes,
    Predictor = predictors,
    Random_effects = random_effects
  ) %>%
    mutate(
      model = pmap(
        list(Gene, Predictor, Random_effects),
        function(gene, predictor, random_effects) {
          fit_arg_betabinomial(
            gene = gene,
            predictor = predictor,
            random_effects = random_effects
          )
        }
      )
    )
}

arg_betabinomial_results <- run_all_arg_betabinomial()


## Extract pooled estimates and adjust for multiple testing ####

arg_betabinomial_summary <- arg_betabinomial_results %>%
  mutate(
    results = map2(
      model,
      Predictor,
      function(mod, predictor) {
        mod$results %>%
          filter(term == predictor) %>%
          select(
            estimate, std.error, statistic, df,
            p.value, OR, OR_low, OR_high
          )
      }
    )
  ) %>%
  select(Gene, Predictor, Random_effects, results) %>%
  unnest(results) %>%
  group_by(Random_effects) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()


## Compare random-effects structures using AIC ####

arg_betabinomial_aic <- arg_betabinomial_results %>%
  mutate(
    aic_by_imp = map(
      model,
      ~ map_dbl(
        .x$fits,
        ~ tryCatch(AIC(.x), error = function(e) NA_real_)
      )
    )
  ) %>%
  select(Gene, Predictor, Random_effects, aic_by_imp) %>%
  unnest_longer(
    aic_by_imp,
    indices_to = "imputation"
  ) %>%
  rename(AIC = aic_by_imp) %>%
  group_by(Gene, Predictor, imputation) %>%
  mutate(
    min_AIC = if (all(!is.finite(AIC))) NA_real_
    else min(AIC[is.finite(AIC)]),
    delta_AIC = AIC - min_AIC
  ) %>%
  ungroup()

arg_betabinomial_aic_summary <- arg_betabinomial_aic %>%
  group_by(Gene, Predictor, Random_effects) %>%
  summarise(
    n_successful = sum(is.finite(AIC)),
    median_AIC = if (all(!is.finite(AIC))) NA_real_
    else median(AIC, na.rm = TRUE),
    median_delta_AIC = if (all(!is.finite(delta_AIC))) NA_real_
    else median(delta_AIC, na.rm = TRUE),
    best_count = sum(
      delta_AIC == min(delta_AIC, na.rm = TRUE),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

arg_betabinomial_aic_summary_wide <- arg_betabinomial_aic_summary %>%
  select(-n_successful, -best_count) %>%
  pivot_wider(
    id_cols = c(Gene, Predictor),
    names_from = Random_effects,
    values_from = c(median_AIC, median_delta_AIC)
  ) %>%
  mutate(across(where(is.numeric), round, digits = 3))

write_tsv(
  arg_betabinomial_aic_summary_wide,
  "supplementary_files/Supplementary_arg_betabinomial_aic_summary.tsv")

arg_betabinomial_aic_winners <- arg_betabinomial_aic_summary %>%
  group_by(Gene, Predictor) %>%
  arrange(median_delta_AIC, .by_group = TRUE) %>%
  mutate(
    second_best_model = nth(Random_effects, 2),
    second_best_median_delta_AIC = nth(median_delta_AIC, 2)
  ) %>%
  filter(
    is.finite(median_delta_AIC),
    median_delta_AIC == min(median_delta_AIC, na.rm = TRUE)
  ) %>%
  ungroup()

# Retain the site-only random effect
selected_betabinomial_models <- arg_betabinomial_results %>%
  filter(Random_effects == "site")


## DHARMa diagnostics for selected beta-binomial models ####
run_dharma_one <- function(fit, n = 500) {

  tryCatch({

    sim <- DHARMa::simulateResiduals(
      fittedModel = fit,
      n = n,
      plot = FALSE
    )

    uniformity <- tryCatch(
      DHARMa::testUniformity(sim, plot = FALSE),
      error = function(e) NULL
    )

    dispersion <- tryCatch(
      DHARMa::testDispersion(sim, plot = FALSE),
      error = function(e) NULL
    )

    outlier <- tryCatch(
      DHARMa::testOutliers(sim, plot = FALSE),
      error = function(e) NULL
    )

    get_stat <- function(x) {
      if (is.null(x)) NA_real_
      else as.numeric(x$statistic[1])
    }

    get_p <- function(x) {
      if (is.null(x)) NA_real_
      else x$p.value
    }

    list(
      simulation = sim,
      statistics = tibble(
        dharma_uniformity_statistic = get_stat(uniformity),
        dharma_uniformity_p = get_p(uniformity),
        dharma_dispersion_statistic = get_stat(dispersion),
        dharma_dispersion_p = get_p(dispersion),
        dharma_outlier_statistic = get_stat(outlier),
        dharma_outlier_p = get_p(outlier),
        dharma_error = NA_character_
      )
    )

  }, error = function(e) {

    list(
      simulation = NULL,
      statistics = tibble(
        dharma_uniformity_statistic = NA_real_,
        dharma_uniformity_p = NA_real_,
        dharma_dispersion_statistic = NA_real_,
        dharma_dispersion_p = NA_real_,
        dharma_outlier_statistic = NA_real_,
        dharma_outlier_p = NA_real_,
        dharma_error = conditionMessage(e)
      )
    )
  })
}

dharma_betabinomial_results <- pmap_dfr(
  selected_betabinomial_models %>%
    select(Gene, Predictor, Random_effects, model),

  function(Gene, Predictor, Random_effects, model) {

    map_dfr(
      seq_along(model$fits),
      function(i) {

        result <- run_dharma_one(model$fits[[i]])

        bind_cols(
          tibble(
            Gene = Gene,
            Predictor = Predictor,
            Random_effects = Random_effects,
            imputation = i,
            dharma_simulation = list(result$simulation)
          ),
          result$statistics
        )
      }
    )
  }
)
dharma_betabinomial_summary <- dharma_betabinomial_results %>%
  group_by(Gene, Predictor, Random_effects) %>%
  summarise(
    median_uniformity_p = if (all(is.na(dharma_uniformity_p))) NA_real_
    else median(dharma_uniformity_p, na.rm = TRUE),
    
    median_dispersion_p = if (all(is.na(dharma_dispersion_p))) NA_real_
    else median(dharma_dispersion_p, na.rm = TRUE),
    
    median_outlier_p = if (all(is.na(dharma_outlier_p))) NA_real_
    else median(dharma_outlier_p, na.rm = TRUE),
    uniformity_problem =
      mean(dharma_uniformity_p < 0.05, na.rm = TRUE) > 0.5,
    dispersion_problem =
      mean(dharma_dispersion_p < 0.05, na.rm = TRUE) > 0.5,
    outlier_problem =
      mean(dharma_outlier_p < 0.05, na.rm = TRUE) > 0.5,
    n_passed = sum(
      !is.na(dharma_uniformity_p) &
        !is.na(dharma_dispersion_p) &
        !is.na(dharma_outlier_p) &
        dharma_uniformity_p >= 0.05 &
        dharma_dispersion_p >= 0.05 &
        dharma_outlier_p >= 0.05
    ),
    .groups = "drop"
  ) %>%
  mutate(
    dharma_note = case_when(
      n_passed == 0 ~
        "DHARMa diagnostics failed",
      uniformity_problem & dispersion_problem & outlier_problem ~
        "Evidence of non-uniformity, dispersion, and outlier problems",
      uniformity_problem & dispersion_problem ~
        "Evidence of non-uniformity and dispersion problems",
      uniformity_problem & outlier_problem ~
        "Evidence of non-uniformity and outlier problems",
      dispersion_problem & outlier_problem ~
        "Evidence of dispersion and outlier problems",
      uniformity_problem ~
        "Evidence of non-uniformity",
      dispersion_problem ~
        "Evidence of dispersion",
      outlier_problem ~
        "Evidence of outlier problems",
      TRUE ~
        "No consistent evidence of DHARMa diagnostic problems"
    )
  )

write_tsv(
  dharma_betabinomial_summary,
  "supplementary_files/Supplementary_arg_dharma_betabinomial_summary.tsv")


# Example: inspect a selected model's simulation if needed.
# Replace the gene and predictor with the association of interest.
# sim <- dharma_betabinomial_results %>%
#   filter(Gene == "blaSHV", Predictor == "DO", imputation == 1) %>%
#   pull(dharma_simulation) %>%
#   pluck(1)
# plot(sim)
# DHARMa::plotQQunif(sim)
# DHARMa::testUniformity(sim)
# DHARMa::testDispersion(sim)
# DHARMa::testOutliers(sim)


## Beta-binomial plots ####
plot_data_betabinomial_final <- arg_betabinomial_summary %>%
  filter(
    Random_effects=="site") %>%
  mutate(
    # Convert significance to a factor with all expected levels
    significance = factor(
      if_else(
        p_adj < 0.05,
        "Significant after BH",
        "Not significant"
      ),
      levels = c("Significant after BH", "Not significant")
    ),
    Gene = factor(Gene, levels = arg_genes),
    Predictor = factor(Predictor, levels = rev(env_predictors))
  )

plot_data_betabinomial_all <- arg_betabinomial_summary %>%
  mutate(
    # Convert significance to a factor with all expected levels
    significance = factor(
      if_else(
        p_adj < 0.05,
        "Significant after BH",
        "Not significant"
      ),
      levels = c("Significant after BH", "Not significant")
    ),
    Gene = factor(Gene, levels = arg_genes),
    Predictor = factor(Predictor, levels = env_predictors),
    Random_effects = factor(
      Random_effects,
      levels = random_effect_structures,
      labels = c(
        "No random effect", "Site",
        "Campaign", "Site + Campaign"
      )
    )
  )

write_tsv(
  plot_data_betabinomial_all %>%
    select(-df),
  "supplementary_files/Supplementary_arg_betabinomial_all_summary.tsv")

# Main beta-binomial plot: Site only random-effects structure.
plot_arg_pc_betabinom_final <- ggplot(
  plot_data_betabinomial_final,
  aes(x = OR, y = Predictor, color = significance)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = OR_low, xmax = OR_high),
    orientation = "y",
    width = 0.3,
    show.legend = TRUE
  ) +
  geom_point(size = 2.5, show.legend = TRUE) +
  facet_wrap(~ Gene, ncol = 3) +
  scale_color_manual(
    values = c(
      "Significant after BH" = "firebrick",
      "Not significant" = "black"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Odds ratio (95% CI)",
    y = "Environmental predictor",
    color = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  "supplementary_files/ARG_PC_betabinomial_final.png",
  plot = plot_arg_pc_betabinom_final,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600
)

# Supplementary model-comparison plot: all random-effects structures.
plot_arg_pc_betabinom_all <- ggplot(
  plot_data_betabinomial_all,
  aes(x = OR, y = Random_effects, color = significance)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = OR_low, xmax = OR_high),
    orientation = "y",
    width = 0.3,
    show.legend = TRUE
  ) +
  geom_point(size = 2, show.legend = TRUE) +
  facet_grid(Gene ~ Predictor) +
  scale_color_manual(
    values = c(
      "Significant after BH" = "firebrick",
      "Not significant" = "black"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Odds ratio (95% CI)",
    y = "Random-effects structure",
    color = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  "supplementary_files/ARG_PC_betabinomial_all.png",
  plot = plot_arg_pc_betabinom_all,
  width = 11,
  height = 6,
  units = "in",
  dpi = 1000
)


# ---- SUPPLEMENTARY ANALYSIS: BINOMIAL GLMMs ----

## Fit binomial models across imputations ####

fit_arg_glmm <- function(
    gene,
    predictor,
    random_effects = random_effect_structures,
    data_list = arg_pc_data
) {

  random_effects <- match.arg(random_effects)

  response <- paste0(
    "cbind(`", gene, "_positive`, ",
    "`", gene, "_tested` - `", gene, "_positive`)"
  )

  random_formula <- switch(
    random_effects,
    none = "",
    site = "+ (1 | Site)",
    campaign = "+ (1 | campaign_id)",
    site_campaign = "+ (1 | Site) + (1 | campaign_id)"
  )

  formula <- as.formula(
    paste(response, "~", predictor, random_formula)
  )

  fits <- map(data_list, function(dat) {

    dat <- dat %>%
      filter(
        .data[[paste0(gene, "_tested")]] > 0,
        is.finite(.data[[predictor]]),
        !is.na(Site),
        !is.na(campaign_id)
      )

    if (random_effects == "none") {

      glm(
        formula,
        data = dat,
        family = binomial(link = "logit")
      )

    } else {

      glmer(
        formula,
        data = dat,
        family = binomial(link = "logit"),
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 200000)
        )
      )
    }
  })

  pooled_fit <- mice::pool(
    mice::as.mira(fits),
    dfcom = Inf
  )

  results <- summary(
    pooled_fit,
    conf.int = TRUE
  ) %>%
    mutate(
      OR = exp(estimate),
      OR_low = exp(`2.5 %`),
      OR_high = exp(`97.5 %`)
    )

  list(
    formula = formula,
    random_effects = random_effects,
    fits = fits,
    pooled = pooled_fit,
    results = results
  )
}


## Run all binomial models ####

run_all_arg_glmms <- function(
    predictors = env_predictors,
    genes = arg_genes,
    random_effects = random_effect_structures
) {

  crossing(
    Gene = genes,
    Predictor = predictors,
    Random_effects = random_effects
  ) %>%
    mutate(
      model = pmap(
        list(Gene, Predictor, Random_effects),
        function(gene, predictor, random_effects) {
          fit_arg_glmm(
            gene = gene,
            predictor = predictor,
            random_effects = random_effects
          )
        }
      )
    )
}

arg_glmm_results <- run_all_arg_glmms()


## Extract pooled estimates and adjust for multiple testing ####

arg_glmm_summary <- arg_glmm_results %>%
  mutate(
    results = map2(
      model,
      Predictor,
      function(mod, predictor) {
        mod$results %>%
          filter(term == predictor) %>%
          select(
            estimate, std.error, statistic, df,
            p.value, OR, OR_low, OR_high
          )
      }
    )
  ) %>%
  select(Gene, Predictor, Random_effects, results) %>%
  unnest(results) %>%
  group_by(Random_effects) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()


## AIC-based selection for supplementary binomial plot ####

arg_glmm_aic <- arg_glmm_results %>%
  mutate(
    aic_by_imp = map(
      model,
      ~ map_dbl(
        .x$fits,
        ~ tryCatch(AIC(.x), error = function(e) NA_real_)
      )
    )
  ) %>%
  select(Gene, Predictor, Random_effects, aic_by_imp) %>%
  unnest_longer(
    aic_by_imp,
    indices_to = "imputation"
  ) %>%
  rename(AIC = aic_by_imp) %>%
  group_by(Gene, Predictor, imputation) %>%
  mutate(
    min_AIC = if (all(!is.finite(AIC))) NA_real_
    else min(AIC[is.finite(AIC)]),
    delta_AIC = AIC - min_AIC
  ) %>%
  ungroup()

arg_glmm_aic_summary <- arg_glmm_aic %>%
  group_by(Gene, Predictor, Random_effects) %>%
  summarise(
    n_successful = sum(is.finite(AIC)),
    median_AIC = if (all(!is.finite(AIC))) NA_real_
    else median(AIC, na.rm = TRUE),
    median_delta_AIC = if (all(!is.finite(delta_AIC))) NA_real_
    else median(delta_AIC, na.rm = TRUE),
    .groups = "drop"
  )

arg_glmm_aic_summary_wide <- arg_glmm_aic_summary %>%
  select(-n_successful) %>%
  pivot_wider(
    id_cols = c(Gene, Predictor),
    names_from = Random_effects,
    values_from = c(median_AIC, median_delta_AIC)
  ) %>%
  mutate(across(where(is.numeric), round, digits = 3))

write_tsv(
  arg_glmm_aic_summary_wide,
  "supplementary_files/Supplementary_arg_binomial_aic_summary.tsv")

arg_glmm_aic_winners <- arg_glmm_aic_summary %>%
  group_by(Gene, Predictor) %>%
  filter(
    is.finite(median_delta_AIC),
    median_delta_AIC == min(median_delta_AIC, na.rm = TRUE)
  ) %>%
  ungroup()


## Binomial supplementary plots ####

plot_data_binomial_final <- arg_glmm_summary %>%
  inner_join(
    arg_glmm_aic_winners %>%
      select(Gene, Predictor, Random_effects),
    by = c("Gene", "Predictor", "Random_effects")
  ) %>%
  mutate(
    # Convert significance to a factor with all expected levels
    significance = factor(
      if_else(
        p_adj < 0.05,
        "Significant after BH",
        "Not significant"
      ),
      levels = c("Significant after BH", "Not significant")
    ),
    Gene = factor(Gene, levels = arg_genes),
    Predictor = factor(Predictor, levels = rev(env_predictors))
  )

plot_data_binomial_all <- arg_glmm_summary %>%
  mutate(
    # Convert significance to a factor with all expected levels
    significance = factor(
      if_else(
        p_adj < 0.05,
        "Significant after BH",
        "Not significant"
      ),
      levels = c("Significant after BH", "Not significant")
    ),
    Gene = factor(Gene, levels = arg_genes),
    Predictor = factor(Predictor, levels = env_predictors),
    Random_effects = factor(
      Random_effects,
      levels = random_effect_structures,
      labels = c(
        "No random effect", "Site",
        "Campaign", "Site + Campaign"
      )
    )
  )

# Supplementary binomial plot using the AIC-selected model.
plot_arg_pc_binom_final <- ggplot(
  plot_data_binomial_final,
  aes(x = OR, y = Predictor, color = significance)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = OR_low, xmax = OR_high),
    orientation = "y",
    width = 0.3
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ Gene, ncol = 3) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Significant after BH" = "firebrick",
      "Not significant" = "black"
    )
  ) +
  labs(
    x = "Odds ratio (95% CI; log scale)",
    y = "Environmental predictor",
    color = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(output_dir, "ARG_PC_binomial_final.png"),
  plot = plot_arg_pc_binom_final,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600
)

# Supplementary plot showing all binomial random-effects structures.
plot_arg_pc_binom_all <- ggplot(
  plot_data_binomial_all,
  aes(x = OR, y = Random_effects, color = significance)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = OR_low, xmax = OR_high),
    orientation = "y",
    width = 0.3
  ) +
  geom_point(size = 2) +
  facet_grid(Gene ~ Predictor) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Significant after BH" = "firebrick",
      "Not significant" = "black"
    )
  ) +
  labs(
    x = "Odds ratio (95% CI; log scale)",
    y = "Random-effects structure",
    color = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  "supplementary_files/ARG_PC_binomial_all.png",
  plot = plot_arg_pc_binom_all,
  width = 11,
  height = 6,
  units = "in",
  dpi = 1000
)
