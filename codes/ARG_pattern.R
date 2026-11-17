library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)
library(lubridate)

library(purrr)
library(lme4)
library(emmeans)
library(DHARMa)
library(broom.mixed)

# ==== Spatial-Seasonal Effect ====
#### Load Data ####
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

# Convert data to binary and prepare for analysis
arg_spatial_data <- arg_analysis_data %>%
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
    Site = factor(Site),
    water_sample_id = factor(water_sample_id),
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
    Season = factor(
      Season,
      levels = c("DRY", "WET")
    )
  )

#### GLMM ####
# function to fit and compare models
fit_spatial_seasonal <- function(gene, data = arg_spatial_data) {
  
  # Retain complete observations
  df <- data %>%
    select(
      all_of(gene),
      Region,
      Season,
      Site,
      water_sample_id
    ) %>%
    filter(complete.cases(.)) %>%
    droplevels()
  
  # Check outcome variation
  if (n_distinct(df[[gene]]) < 2) {
    return(list(
      result = tibble(
        ARG = gene,
        n = nrow(df),
        n_samples = n_distinct(df$water_sample_id),
        AIC_M1 = NA_real_,
        AIC_M2 = NA_real_,
        LRT_statistic = NA_real_,
        df = NA_real_,
        p_value = NA_real_,
        selected_model = NA_character_,
        note = "Outcome has no variation"
      ),
      models = NULL,
      data = df
    ))
  }
  
  # Standardize outcome name
  df <- df %>%
    rename(outcome = all_of(gene))
  
  # M1: Additive model
  model_M1 <- glmer(
    outcome ~ Region + Season +
      (1 | water_sample_id),
    data = df,
    family = binomial(link = "logit"),
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
  
  # M2: Interaction model
  model_M2 <- glmer(
    outcome ~ Region * Season +
      (1 | water_sample_id),
    data = df,
    family = binomial(link = "logit"),
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
  
  # Likelihood-ratio test
  lrt <- anova(model_M1, model_M2, test = "Chisq")
  
  result <- tibble(
    ARG = gene,
    n = nobs(model_M1),
    n_samples = n_distinct(df$water_sample_id),
    AIC_M1 = AIC(model_M1),
    AIC_M2 = AIC(model_M2),
    LRT_statistic = lrt$Chisq[2],
    df = lrt$`Chi Df`[2],
    p_value = lrt$`Pr(>Chisq)`[2],
    selected_model = ifelse(
      lrt$`Pr(>Chisq)`[2] < 0.05,
      "M2: Interaction",
      "M1: Additive"
    )
  )
  
  # Return results, fitted models, and analysis data
  list(
    result = result,
    models = list(
      M1_additive = model_M1,
      M2_interaction = model_M2
    ),
    data = df
  )
}

# Run the analysis for all six ARGs
spatial_seasonal_fits <- set_names(
  map(arg_cols, fit_spatial_seasonal),
  arg_cols
)
summary(spatial_seasonal_fits[["sul1"]]$models$M1_additive)

#Extract the summary results:
spatial_seasonal_results <- map_dfr(
  spatial_seasonal_fits,
  "result"
) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adjusted)

# Check LRT for additive model
test_additive_effects <- function(fit_object) {
  
  # Extract additive model and its analysis data
  model_full <- fit_object$models$M1_additive
  df <- fit_object$data
  
  # Reduced model without Region
  model_no_region <- update(
    model_full,
    . ~ . - Region
  )
  
  # Reduced model without Season
  model_no_season <- update(
    model_full,
    . ~ . - Season
  )
  
  # LRTs
  lrt_region <- anova(
    model_no_region,
    model_full,
    test = "Chisq"
  )
  
  lrt_season <- anova(
    model_no_season,
    model_full,
    test = "Chisq"
  )
  
  tibble(
    ARG = fit_object$result$ARG,
    
    n = nobs(model_full),
    n_samples = n_distinct(df$water_sample_id),
    
    AIC = AIC(model_full),
    
    region_LRT = lrt_region$Chisq[2],
    region_df = lrt_region$`Chi Df`[2],
    region_p = lrt_region$`Pr(>Chisq)`[2],
    
    season_LRT = lrt_season$Chisq[2],
    season_df = lrt_season$`Chi Df`[2],
    season_p = lrt_season$`Pr(>Chisq)`[2]
  )
}

# Run across all six ARGs
additive_effect_results <- map_dfr(
  spatial_seasonal_fits,
  test_additive_effects
) %>%
  mutate(
    region_p_adjusted = p.adjust(region_p, method = "BH"),
    season_p_adjusted = p.adjust(season_p, method = "BH")
  ) %>%
  mutate(
    effect_classification = case_when(
      region_p_adjusted < 0.05 &
        season_p_adjusted < 0.05 ~ "Both Region and Season",
      
      region_p_adjusted < 0.05 ~ "Region only",
      
      season_p_adjusted < 0.05 ~ "Season only",
      
      TRUE ~ "Neither"
    )
  )

#### Plot ####
# Extract emmeans for all ARGs
emm_results <- map(
  arg_cols,
  function(gene) {
    
    model <- spatial_seasonal_fits[[gene]]$models$M1_additive
    
    emm <- emmeans(
      model,
      ~ Region * Season,
      type = "response"
    )
    
    as.data.frame(emm) %>%
      mutate(Gene = gene)
  }
)

# Combine results
emm_all <- bind_rows(emm_results)

emm_all$Region <- factor(
  emm_all$Region,
  levels = c(
    "East Lake", "West Lake",
    "Quadrat 1", "Quadrat 2",
    "Quadrat 3", "Quadrat 4"
  )
)

emm_all$Gene <- factor(
  emm_all$Gene,
  levels = arg_cols
)

ggplot(
  emm_all,
  aes(
    x = Region,
    y = prob,
    color = Season,
    group = Season
  )
) +
  geom_point(
    position = position_dodge(0.3),
    size = 2.5
  ) +
  geom_line(
    position = position_dodge(0.3)
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.2,
    position = position_dodge(0.3)
  ) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    x = "Region",
    y = "Predicted probability of ARG-positive isolates",
    color = "Season"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )



#-----New---------
#### Load Data ####
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

# Convert data to binary and prepare for analysis
arg_spatial_data <- arg_analysis_data %>%
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
    Site = factor(Site),
    campaign_id = factor(campaign_id),
    water_sample_id = factor(water_sample_id),
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
    Season = factor(
      Season,
      levels = c("DRY", "WET")
    )
  )

#### Function 1: Fit GLMM ####

fit_one_glmm <- function(
    gene,
    model_name,
    formula,
    data = arg_spatial_data
) {
  
  # Retain complete observations for consistent model comparisons
  df <- data %>%
    select(
      all_of(gene),
      Region,
      Season,
      Site,
      campaign_id,
      water_sample_id
    ) %>%
    filter(complete.cases(.)) %>%
    droplevels() %>%
    rename(outcome = all_of(gene))
  
  # Standard empty result for failed models
  empty_result <- function(note) {
    list(
      coefficients = tibble(),
      model = NULL,
      data = df,
      status = note
    )
  }
  
  # Check outcome variation
  if (nrow(df) == 0) {
    return(empty_result("No complete observations"))
  }
  
  if (n_distinct(df$outcome) < 2) {
    return(empty_result("Outcome has no variation"))
  }
  
  # Capture model warnings
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
  
  # Handle fitting errors
  if (inherits(model, "error")) {
    return(empty_result(
      paste("Model error:", conditionMessage(model))
    ))
  }
  
  # Model diagnostics
  singular_fit <- isSingular(model, tol = 1e-4)
  
  opt_code <- model@optinfo$conv$opt
  conv_messages <- model@optinfo$conv$lme4$messages
  
  notes <- character()
  
  if (singular_fit) {
    notes <- c(notes, "Singular fit")
  }
  
  if (
    length(opt_code) > 0 &&
    any(opt_code != 0, na.rm = TRUE)
  ) {
    notes <- c(notes, "Optimizer convergence issue")
  }
  
  if (!is.null(conv_messages) && length(conv_messages) > 0) {
    notes <- c(notes, paste(
      "Convergence message:",
      paste(conv_messages, collapse = " | ")
    ))
  }
  
  if (length(fit_warnings) > 0) {
    notes <- c(
      notes,
      paste(
        "Model warning:",
        paste(unique(fit_warnings), collapse = " | ")
      )
    )
  }
  
  model_status <- if (length(notes) == 0) {
    "No convergence issues"
  } else {
    paste(unique(notes), collapse = "; ")
  }
  
  # Extract fixed-effect estimates
  coefs <- as.data.frame(
    summary(model)$coefficients
  ) %>%
    rownames_to_column("term") %>%
    as_tibble()
  
  # Extract variance components
  vc <- as.data.frame(VarCorr(model))
  
  variance_components <- paste(
    paste0(
      vc$grp,
      "=",
      signif(vc$vcov, 4)
    ),
    collapse = "; "
  )
  
  # Extract model-level information first
  model_n <- stats::nobs(model)
  model_n_samples <- dplyr::n_distinct(df$water_sample_id)
  
  # Create coefficient-level results
  coefficient_results <- coefs %>%
    transmute(
      ARG = gene,
      model = model_name,
      term,
      n = model_n,
      n_samples = model_n_samples,
      OR = exp(Estimate),
      CI_low = exp(Estimate - 1.96 * `Std. Error`),
      CI_high = exp(Estimate + 1.96 * `Std. Error`),
      p_value = `Pr(>|z|)`,
      variance_components = variance_components,
      model_status = model_status
    )
  
  list(
    coefficients = coefficient_results,
    model = model,
    data = df,
    status = model_status
  )
}

#### Model specifications ####

primary_model_specs <- list(
  M1_additive = outcome ~ Region + Season + (1 | water_sample_id),
  M2_interaction = outcome ~ Region * Season + (1 | water_sample_id),
  M3_without_region = outcome ~ Season + (1 | water_sample_id),
  M4_without_season = outcome ~ Region + (1 | water_sample_id)
)

sensitivity_model_specs <- list(
  M1_additive = outcome ~ Region + Season + (1 | water_sample_id),
  M5_sample_site = outcome ~ Region + Season +
    (1 | water_sample_id) + (1 | Site),
  M6_sample_campaign = outcome ~ Region + Season +
    (1 | water_sample_id) + (1 | campaign_id),
  M7_sample_site_campaign = outcome ~ Region + Season +
    (1 | water_sample_id) + (1 | Site) + (1 | campaign_id)
)

#### Fit model sets ####

fit_model_set <- function(
    model_specs,
    genes = arg_cols,
    data = arg_spatial_data
) {
  
  fits <- map(
    genes,
    function(gene) {
      
      map(
        names(model_specs),
        function(model_name) {
          fit_one_glmm(
            gene = gene,
            model_name = model_name,
            formula = model_specs[[model_name]],
            data = data
          )
        }
      ) %>%
        set_names(names(model_specs))
    }
  ) %>%
    set_names(genes)
  
  fits
}

# Run the primary models:
primary_fits <- fit_model_set(
  primary_model_specs,
  arg_cols,
  arg_spatial_data
)

# Run the sensitivity models:
sensitivity_fits <- fit_model_set(
  sensitivity_model_specs
)


# access any fitted model using:
primary_fits[["blaCTX-M"]]$M1_additive$model
primary_fits[["blaCTX-M"]]$M3_without_region$model
sensitivity_fits[["blaCTX-M"]]$M5_sample_site$model

#### LRT comparisons ####

run_spatial_lrt <- function(fit_object) {
  
  model_M1 <- fit_object$M1_additive$model
  model_M2 <- fit_object$M2_interaction$model
  
  df <- fit_object$M1_additive$data
  
  if (is.null(model_M1) || is.null(model_M2)) {
    return(tibble(
      interaction_LRT = NA_real_,
      interaction_df = NA_real_,
      interaction_p = NA_real_,
      region_LRT = NA_real_,
      region_df = NA_real_,
      region_p = NA_real_,
      season_LRT = NA_real_,
      season_df = NA_real_,
      season_p = NA_real_
    ))
  }
  
  # Reduced model without Region
  model_no_region <- glmer(
    outcome ~ Season +
      (1 | water_sample_id),
    data = df,
    family = binomial(link = "logit"),
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
  
  # Reduced model without Season
  model_no_season <- glmer(
    outcome ~ Region +
      (1 | water_sample_id),
    data = df,
    family = binomial(link = "logit"),
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
  
  # Interaction test: M1 vs M2
  lrt_interaction <- anova(
    model_M1,
    model_M2,
    test = "Chisq"
  )
  
  # Overall Region effect
  lrt_region <- anova(
    model_no_region,
    model_M1,
    test = "Chisq"
  )
  
  # Overall Season effect
  lrt_season <- anova(
    model_no_season,
    model_M1,
    test = "Chisq"
  )
  
  tibble(
    interaction_LRT = lrt_interaction$Chisq[2],
    interaction_df = lrt_interaction$`Chi Df`[2],
    interaction_p = lrt_interaction$`Pr(>Chisq)`[2],
    
    region_LRT = lrt_region$Chisq[2],
    region_df = lrt_region$`Chi Df`[2],
    region_p = lrt_region$`Pr(>Chisq)`[2],
    
    season_LRT = lrt_season$Chisq[2],
    season_df = lrt_season$`Chi Df`[2],
    season_p = lrt_season$`Pr(>Chisq)`[2]
  )
}

# Run across all ARGs
spatial_lrt_results <- map_dfr(
  arg_cols,
  function(gene) {
    
    result <- run_spatial_lrt(
      primary_fits[[gene]]
    )
    
    result %>%
      mutate(ARG = gene, .before = 1)
  }
) %>%
  mutate(
    interaction_p_adjusted = p.adjust(
      interaction_p,
      method = "BH"
    ),
    region_p_adjusted = p.adjust(
      region_p,
      method = "BH"
    ),
    season_p_adjusted = p.adjust(
      season_p,
      method = "BH"
    )
  )

supplementary_lrt <- spatial_lrt_results %>%
  relocate(interaction_p_adjusted, .after = interaction_p) %>%
  relocate(region_p_adjusted, .after = region_p)


#### Create supplementary model tables ####

create_supplementary_table <- function(fits) {
  
  results <- map_dfr(
    names(fits),
    function(gene) {
      
      map_dfr(
        names(fits[[gene]]),
        function(model_name) {
          
          fits[[gene]][[model_name]]$coefficients
        }
      )
    }
  )
  
  results %>%
    group_by(model, term) %>%
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    ungroup() %>%
    select(
      ARG,
      model,
      term,
      n,
      n_samples,
      OR,
      CI_low,
      CI_high,
      p_value,
      p_adjusted,
      variance_components,
      model_status
    ) %>%
    arrange(model, term, ARG)
}

## Generate tables
# Primary models: M1, M2, M3, M4
supplementary_primary <- create_supplementary_table(
  primary_fits
)

supplementary_primary <- supplementary_primary %>%
  arrange(model, ARG, term) %>%
  mutate(
    across(
      where(is.double),
      ~ case_when(
        is.na(.x) ~ NA_character_,
        .x == 0 ~ "0.000",
        abs(.x) < 0.001 ~ sprintf("%.2e", .x),
        TRUE ~ sprintf("%.3f", .x)
      )
    ),
    term = recode(
      term,
      "(Intercept)" = "Intercept"
    )
  ) %>%
  relocate(model, ARG, variance_components, model_status) %>%
  select(-n, -n_samples)

# Sensitivity models: M1, M5, M6, M7
supplementary_sensitivity <- create_supplementary_table(
  sensitivity_fits
)

supplementary_sensitivity <- supplementary_sensitivity %>%
  arrange(model, ARG, term) %>%
  mutate(
    across(
      where(is.double),
      ~ case_when(
        is.na(.x) ~ NA_character_,
        .x == 0 ~ "0.000",
        abs(.x) < 0.001 ~ sprintf("%.2e", .x),
        TRUE ~ sprintf("%.3f", .x)
      )
    ),
    term = recode(
      term,
      "(Intercept)" = "Intercept"
    )
  ) %>%
  relocate(model, ARG, variance_components, model_status) %>%
  select(-n, n_samples)

# Export files
write_tsv(
  supplementary_lrt,
  "supplementary_files/Supplementary_ARG_spatial_seasonal_LRT_results.tsv")

write_tsv(
  supplementary_primary,
  "supplementary_files/Supplementary_ARG_spatial_seasonal_primary.tsv"
)

write_tsv(
  supplementary_sensitivity,
  "supplementary_files/Supplementary_ARG_spatial_seasonal_sensitivity.tsv"
)

#### Main-text figure: Estimated marginal means ####

emm_results <- map(
  arg_cols,
  function(gene) {
    
    model <- primary_fits[[gene]]$M1_additive$model
    
    emm <- emmeans(
      model,
      ~ Region * Season,
      type = "response"
    )
    
    as.data.frame(emm) %>%
      mutate(Gene = gene)
  }
)

emm_all <- bind_rows(emm_results) %>%
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
    Gene = factor(
      Gene,
      levels = arg_cols
    )
  )

emm_annotations <- tibble(
  Gene = c("blaCTX-M", "blaSHV"),
  label = c(
    "Season effect: adjusted p = 0.020",
    "Season effect: adjusted p = 0.020"
  )
)

# Plot
p_emmeans <- ggplot(
  emm_all,
  aes(
    x = Region,
    y = prob,
    color = Season,
    group = Season
  )
) +
  geom_point(
    position = position_dodge(0.3),
    size = 2.5
  ) +
  geom_line(
    position = position_dodge(0.3)
  ) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.2,
    position = position_dodge(0.3)
  ) +
  facet_wrap(
    ~ Gene,
    scales = "free_y"
  ) +
  labs(
    x = "Region",
    y = "Predicted probability of ARG-positive isolates",
    color = "Season"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  geom_text(
    data = emm_annotations,
    aes(
      x = Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.5,
    size = 3.0
  )

ggsave(
  "supplementary_files/ARG_prevalence_emmeans.png",
  plot = p_emmeans,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600
)

#### DHARMa diagnostics: M1 additive model only #### 

dharma_simulations <- list()

dharma_results <- purrr::map_dfr(seq_along(arg_cols), function(i) {
  
  gene <- arg_cols[i]
  
  model <- primary_fits[[gene]][["M1_additive"]][["model"]]
  
  if (is.null(model) || !inherits(model, "merMod")) {
    return(tibble(
      ARG = gene,
      dharma_uniformity_statistic = NA_real_,
      dharma_uniformity_p = NA_real_,
      dharma_dispersion_statistic = NA_real_,
      dharma_dispersion_p = NA_real_,
      dharma_outlier_statistic = NA_real_,
      dharma_outlier_p = NA_real_,
      dharma_note = "DHARMa diagnostics unavailable: invalid model"
    ))
  }
  
  set.seed(123 + i)
  
  sim_res <- DHARMa::simulateResiduals(
    fittedModel = model,
    n = 1000,
    plot = FALSE
  )
  
  # Save simulation object for later exploration
  dharma_simulations[[gene]] <<- sim_res
  
  # Run diagnostic tests
  uniformity_test <- DHARMa::testUniformity(sim_res)
  dispersion_test <- DHARMa::testDispersion(sim_res)
  outlier_test <- DHARMa::testOutliers(sim_res)
  
  # Extract statistics and p-values
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
  
  uniformity_stat <- get_stat(uniformity_test)
  uniformity_p <- get_p(uniformity_test)
  
  dispersion_stat <- get_stat(dispersion_test)
  dispersion_p <- get_p(dispersion_test)
  
  outlier_stat <- get_stat(outlier_test)
  outlier_p <- get_p(outlier_test)
  
  # Summarize diagnostic results
  significant_tests <- c(
    if (!is.na(uniformity_p) && uniformity_p < 0.05)
      "Uniformity test",
    if (!is.na(dispersion_p) && dispersion_p < 0.05)
      "Dispersion test",
    if (!is.na(outlier_p) && outlier_p < 0.05)
      "Outlier test"
  )
  
  unavailable_tests <- c(
    if (is.na(uniformity_p)) "Uniformity test",
    if (is.na(dispersion_p)) "Dispersion test",
    if (is.na(outlier_p)) "Outlier test"
  )
  
  dharma_note <- dplyr::case_when(
    length(significant_tests) > 0 ~
      paste("Significant:", paste(significant_tests, collapse = ", ")),
    
    length(unavailable_tests) > 0 ~
      paste("Unavailable:", paste(unavailable_tests, collapse = ", ")),
    
    TRUE ~ "No significant deviations"
  )
  
  tibble(
    ARG = gene,
    dharma_uniformity_statistic = uniformity_stat,
    dharma_uniformity_p = uniformity_p,
    dharma_dispersion_statistic = dispersion_stat,
    dharma_dispersion_p = dispersion_p,
    dharma_outlier_statistic = outlier_stat,
    dharma_outlier_p = outlier_p,
    dharma_note = dharma_note
  )
})


#### Sensitivity model comparison ####

calculate_sensitivity_aic <- function(sensitivity_fits) {
  
  purrr::map_dfr(names(sensitivity_fits), function(gene) {
    
    purrr::map_dfr(names(sensitivity_fits[[gene]]), function(model_name) {
      
      fit_object <- sensitivity_fits[[gene]][[model_name]]
      fitted_model <- fit_object[["model"]]
      
      if (is.null(fitted_model) ||
          !inherits(fitted_model, "merMod")) {
        
        return(tibble(
          ARG = gene,
          model = model_name,
          AIC = NA_real_,
          n = NA_integer_,
          model_status = "Model unavailable or invalid"
        ))
      }
      
      tibble(
        ARG = gene,
        model = model_name,
        AIC = stats::AIC(fitted_model),
        n = stats::nobs(fitted_model),
        model_status = fit_object$status
      )
    })
  }) %>%
    group_by(ARG) %>%
    mutate(
      delta_AIC = AIC - min(AIC, na.rm = TRUE)
    ) %>%
    ungroup()
}

sensitivity_aic <- calculate_sensitivity_aic(sensitivity_fits)

supplementary_sensitivity_aic <- sensitivity_aic %>%
  select(-n) %>%
  left_join(
    dharma_results,
    by = "ARG"
  ) %>%
  mutate(
    across(
      starts_with("dharma_"),
      ~ if_else(
        model %in% c(
          "M5_sample_site",
          "M6_sample_campaign",
          "M7_sample_site_campaign"
        ),
        NA,
        .
      )
    ),
    dharma_note = if_else(
      model %in% c(
        "M5_sample_site",
        "M6_sample_campaign",
        "M7_sample_site_campaign"
      ),
      "DHARMa diagnostics not run",
      dharma_note
    ),
    dharma_outlier_statistic = as.integer(dharma_outlier_statistic)
  ) %>%
  arrange(ARG) %>%
  mutate(
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

## Explore a specific flagged model
sim_res <- dharma_simulations[["tetA"]]
DHARMa::testOutliers(sim_res, margin = "both")


write_tsv(
  supplementary_sensitivity_aic,
  "supplementary_files/Supplementary_ARG_spatial_seasonal_AIC_DHARMa.tsv"
)

  