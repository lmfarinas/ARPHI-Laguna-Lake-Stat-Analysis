library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(lubridate)

# ==== Standardize ARG Data ====
## LOAD ORIGINAL DATA ####
sites <-read.table("raw_data/sites.tsv", header=TRUE, sep="\t")
Lake_wet <- read.table("raw_data/ARPHI-ARG-Database - Lake-Wet.tsv", header=TRUE, sep="\t",
                       na.strings = c("", "N/A"), stringsAsFactors = FALSE)
Lake_dry <- read.table("raw_data/ARPHI-ARG-Database - Lake-Dry.tsv", header=TRUE, sep="\t",
                       na.strings = c("", "N/A"), stringsAsFactors = FALSE) 
Trib_wet <-read.table("raw_data/ARPHI-ARG-Database - Tributary-Wet.tsv", header=TRUE, sep="\t",
                      na.strings = c("", "N/A"), stringsAsFactors = FALSE)
Trib_dry <-read.table("raw_data/ARPHI-ARG-Database - Tributary-Dry.tsv", header=TRUE, sep="\t",
                      na.strings = c("", "N/A"), stringsAsFactors = FALSE)

## COMBINE AND STANDARDIZE #### 

Combined_arg_data <- rbind(Lake_wet,Lake_dry,Trib_wet,Trib_dry) %>%
  left_join(sites,by="Site") %>%
  relocate(Region, .after = Site) %>%
  rename(
    isolate_id = Code,
    `blaCTX-M` = blaCTX.M,
    `blaOXA-1` = blaOXA.1,
    colony_characteristics = Colony.characteristics..EMBA.,
    sample_collection_date = Sample.Collection.Date,
    pcr_confirmation = Status
  ) %>%
  mutate(
    isolate_id = case_when(
      Season == "WET" ~ paste0("W-", isolate_id),
      Season == "DRY" ~ paste0("D-", isolate_id),
      TRUE ~ isolate_id
    ),
    sample_collection_date = as.Date(
      sample_collection_date,
      format = "%m/%d/%Y"
    ),
    water_sample_id = paste(
      substr(Season, 1, 1),
      Site,
      format(sample_collection_date, "%Y%m%d"),
      sep = "-"
    )
  ) %>%
  relocate(water_sample_id, isolate_id)


# Build a chronological table of sampling dates and seasons
campaign_map <- Combined_arg_data %>%
  distinct(sample_collection_date, Season) %>%
  mutate(
    Year = format(sample_collection_date, "%Y"),
    Season = toupper(Season)
  ) %>%
  arrange(sample_collection_date) %>%
  group_by(Year) %>%
  mutate(
    season_block = cumsum(
      Season != lag(Season, default = first(Season))
    ) + 1L
  ) %>%
  ungroup() %>%
  group_by(Year, Season, season_block) %>%
  summarise(
    campaign_start = min(sample_collection_date),
    campaign_end = max(sample_collection_date),
    .groups = "drop"
  ) %>%
  arrange(Year, campaign_start) %>%
  group_by(Year, Season) %>%
  mutate(campaign_number = row_number()) %>%
  ungroup() %>%
  mutate(
    campaign_id = paste(Year, Season, campaign_number, sep = "-")
  )

date_campaign_map <- Combined_arg_data %>%
  distinct(sample_collection_date, Season) %>%
  mutate(
    Year = format(sample_collection_date, "%Y"),
    Season = toupper(Season)
  ) %>%
  arrange(sample_collection_date) %>%
  group_by(Year) %>%
  mutate(
    season_block = cumsum(
      Season != lag(Season, default = first(Season))
    ) + 1L
  ) %>%
  ungroup() %>%
  left_join(
    campaign_map %>%
      select(Year, Season, season_block, campaign_id),
    by = c("Year", "Season", "season_block")
  )

Combined_arg_data <- Combined_arg_data %>%
  select(-any_of(c("campaign_id", "Year", "season_block"))) %>%
  mutate(
    Season = toupper(Season),
    Year = format(sample_collection_date, "%Y")
  ) %>%
  left_join(
    date_campaign_map %>%
      select(sample_collection_date, Season, campaign_id),
    by = c("sample_collection_date", "Season")
  ) %>%
  relocate(Year, campaign_id)

# ARG columns
arg_cols <- c(
  "blaCTX-M", "blaSHV", "blaTEM",
  "blaNDM", "blaOXA-1", "tetA",
  "sul1", "sul2"
)

## CREATE SUPPLEMENTARY FILES ####

### Sample Metadata Summary ####
sample_metadata <- Combined_arg_data %>%
  group_by(
    campaign_id,
    water_sample_id,
    sample_collection_date,
    Season,
    Site,
    Region
  ) %>%
  summarise(
    n_isolates = n_distinct(isolate_id, na.rm = TRUE),
    
    n_uidA_positive = n_distinct(
      isolate_id[pcr_confirmation %in% "uidA positive"],
      na.rm = TRUE
    ),
    
    n_uidA_negative = n_distinct(
      isolate_id[pcr_confirmation %in% "uidA negative"],
      na.rm = TRUE
    ),
    
    n_uidA_untested = n_distinct(
      isolate_id[is.na(pcr_confirmation)],
      na.rm = TRUE
    ),
    
    n_isolates_arg_analysis = n_distinct(
      isolate_id[
        !is.na(pcr_confirmation) &
          if_all(all_of(arg_cols), ~ !is.na(.x))
      ],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(sample_collection_date, Site)


### Site-Season Distribution ####
# Define explicit region order
region_levels <- c(
  "East Lake", 
  "West Lake", 
  "Quadrat 1", 
  "Quadrat 2", 
  "Quadrat 3", 
  "Quadrat 4"
)

site_season_summary <- sample_metadata %>%
  group_by(Region, Site, Season) %>%
  summarise(
    n_water_samples = n_distinct(water_sample_id),
    n_samples_arg_analysis = n_distinct(
      water_sample_id[n_isolates_arg_analysis > 0]
    ),
    n_isolates = sum(n_isolates),
    n_isolates_arg_analysis = sum(n_isolates_arg_analysis),
    .groups = "drop"
  ) %>%
  # Set explicit Region factor order
  mutate(
    Region = factor(Region, levels = region_levels)
  ) %>%
  # Sort by Region level, then site prefix/number
  arrange(Region, substr(Site, 1, 1), parse_number(Site)) %>%
  # Set Site factor levels according to sorted order
  mutate(
    Site = factor(Site, levels = unique(Site))
  ) %>%
  mutate(
    Region = forcats::fct_recode(
      Region,
      "East\nLake" = "East Lake",
      "West\nLake" = "West Lake",
      "Quadrat\n1" = "Quadrat 1",
      "Quadrat\n2" = "Quadrat 2",
      "Quadrat\n3" = "Quadrat 3",
      "Quadrat\n4" = "Quadrat 4"
    )
  )


plot_samples <- ggplot(site_season_summary, aes(x = Site)) +
  
  # Background: total collected
  geom_col(
    aes(y = n_water_samples, fill = "Total collected"),
    color = "grey40",
    width = 0.7,
    linewidth = 0.3
  ) +
  
  # Foreground: included in ARG analysis
  geom_col(
    aes(y = n_samples_arg_analysis, fill = "Included in ARG analysis"),
    color = "black",
    width = 0.7,
    linewidth = 0.3
  ) +
  
  # Total count labels
  geom_text(
    aes(y = n_water_samples, label = n_water_samples),
    vjust = 0.5,
    hjust = -0.5,
    angle = 90,
    size = 3
  ) +
  
  # ARG analysis count labels
  geom_text(
    aes(y = n_samples_arg_analysis, label = n_samples_arg_analysis),
    vjust = 0.5,
    hjust = 1.5,
    angle = 90,
    color = "white",
    fontface = "bold",
    size = 3
  ) +
  
  # Nested faceting
  ggh4x::facet_nested(
    ~ Season + Region,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    name = NULL,
    values = c(
      "Total collected" = "grey85",
      "Included in ARG analysis" = "black"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, by = 2),
    labels = as.integer,
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    x = "Sampling site",
    y = "Number of water samples"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom",
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    
    strip.text.x = element_text(
      face = "bold",
      size = 9
    ),
    
    # Rotate only the Region labels
    ggh4x.facet.nestline.x = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    # Give the x-axis and facet strips more room
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 20,
      l = 15
    )
  )

plot_isolates <- ggplot(site_season_summary, aes(x = Site)) +
  
  # Background: total collected
  geom_col(
    aes(y = n_isolates, fill = "Total collected"),
    color = "grey40",
    width = 0.7,
    linewidth = 0.3
  ) +
  
  # Foreground: included in ARG analysis
  geom_col(
    aes(y = n_isolates_arg_analysis, fill = "Included in ARG analysis"),
    color = "black",
    width = 0.7,
    linewidth = 0.3
  ) +
  
  # Total count labels
  geom_text(
    aes(y = n_isolates, label = n_isolates),
    vjust = 0.5,
    hjust = -0.2,
    angle = 90,
    size = 3
  ) +
  
  # ARG analysis count labels
  geom_text(
    aes(y = n_isolates_arg_analysis, label = n_isolates_arg_analysis),
    vjust = 0.5,
    hjust = 1.2,
    angle = 90,
    color = "white",
    fontface = "bold",
    size = 3
  ) +
  
  # Nested faceting: Season and Region
  ggh4x::facet_nested(
    ~ Season + Region,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    name = NULL,
    values = c(
      "Total collected" = "grey85",
      "Included in ARG analysis" = "black"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  
  labs(
    x = "Sampling site",
    y = "Number of isolates"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom",
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    
    strip.text.x = element_text(
      face = "bold",
      size = 9
    ),
    
    ggh4x.facet.nestline.x = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 20,
      l = 15
    )
  )

ggsave("supplementary_files/sample_distribution.png",
       plot = plot_samples,
       width = 8.5, height = 4, dpi = 600)
ggsave("supplementary_files/isolate_distribution.png",
       plot = plot_isolates,
       width = 8.5, height = 4.5, dpi = 600)


site_season_wide <- site_season_summary %>%
  pivot_wider(
    names_from = Season,
    values_from = c(
      n_water_samples,
      n_isolates,
      n_isolates_arg_analysis
    ),
    names_glue = "{.value}_{Season}",
    values_fill = 0
  )


### Sampling Coverage ####
site_campaign_summary <- sample_metadata %>%
  group_by(campaign_id, Site) %>%
  summarise(
    n_water_samples = n_distinct(water_sample_id),
    .groups = "drop"
  )

# Include unsampled site-campaign combinations as zero
sampling_coverage <- site_campaign_summary %>%
  complete(
    campaign_id,
    Site,
    fill = list(n_water_samples = 0)
  ) %>%
  mutate(
    campaign_id = factor(
      campaign_id,
      levels = unique(sample_metadata$campaign_id[
        order(sample_metadata$sample_collection_date)
      ])
    ),
    Site = factor(
      Site,
      levels = unique(rev(site_season_wide$Site))
    )
  )

plot_sampling_cov <- ggplot(
  sampling_coverage,
  aes(x = campaign_id, y = Site, fill = n_water_samples)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "C",
    name = "Collected water samples",
    breaks = scales::pretty_breaks()
  ) +
  labs(
    x = "Sampling campaign",
    y = "Sampling site"
    # title = "Sampling coverage across sites and campaigns"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave("supplementary_files/sampling_coverage.png",
       plot = plot_sampling_cov,
       width = 7, height = 4, dpi = 600)

## Output tables

write_tsv(Combined_arg_data, file = "ARPHI-isolate-database.tsv")


# ==== Standardize Physicochemical Data ====

#### LOAD ORIGINAL DATA ####
trib_dry <- read.table("raw_data/Physico-chem Data - Trib DRY v.1.tsv", header=TRUE, sep="\t") %>%
  rename(sample_collection_date=2, BOD=4, DO=5, FC=6, pH=7, Ammonia=8, Nitrate=9, Phosphate=10) %>%
  mutate(Season="DRY",
         FC = as.numeric(gsub(",", "", FC)),
         sample_collection_date = as.Date(
           sample_collection_date,
           format = "%m/%d/%Y"
         )) %>%
  mutate(across(c("BOD", "DO", "FC", "pH", "Ammonia", "Nitrate", "Phosphate"), ~ {
    .x %>%
      as.character()               # Convert to string
  }))
trib_wet <- read.table("raw_data/Physico-chem Data - Trib WET v.1.tsv", header=TRUE, sep="\t") %>%
  rename(sample_collection_date=2, BOD=4, DO=5, FC=6, pH=7, Ammonia=8, Nitrate=9, Phosphate=10) %>%
  mutate(Site = str_replace(Site, "R", ""),
         Season="WET",
         FC = as.numeric(gsub(",", "", FC)),
         sample_collection_date = as.Date(
           sample_collection_date,
           format = "%m/%d/%Y"
         )) %>%
  mutate(across(c("BOD", "DO", "FC", "pH", "Ammonia", "Nitrate", "Phosphate"), ~ {
    .x %>%
      as.character()               # Convert to string
  }))

lake_dry <- read.table("raw_data/Physico-chem Data - Lake DRY v.1.tsv", header=TRUE, sep="\t") %>%
  rename(sample_collection_date=2, BOD=4, DO=5, FC=6, pH=7, Ammonia=8, Nitrate=9, Phosphate=10) %>%
  mutate(Season="DRY",
         FC = as.numeric(gsub(",", "", FC)),
         sample_collection_date = as.Date(
           sample_collection_date,
           format = "%m/%d/%Y"
         )) %>%
  select(-TSS) %>%
  mutate(across(c("BOD", "DO", "FC", "pH", "Ammonia", "Nitrate", "Phosphate"), ~ {
    .x %>%
      as.character()               # Convert to string
  }))

lake_wet <- read.table("raw_data/Physico-chem Data - Lake WET v.1.tsv", header=TRUE, sep="\t")%>%
  rename(sample_collection_date=2, BOD=4, DO=5, FC=6, pH=7, Ammonia=8, Nitrate=9, Phosphate=10) %>%
  mutate(Season="WET",
         FC = as.numeric(gsub(",", "", FC)),
         sample_collection_date = as.Date(
           sample_collection_date,
           format = "%m/%d/%Y"
         )) %>%
  select(-TSS) %>%
  mutate(across(c("BOD", "DO", "FC", "pH", "Ammonia", "Nitrate", "Phosphate"), ~ {
    .x %>%
      as.character()               # Convert to string
  }))

# Combine Data
combined_pc_data_raw <- bind_rows(lake_dry, lake_wet, trib_dry, trib_wet) %>%
  left_join(sites, by="Site") %>%
  arrange(sample_collection_date)

## Fix colelction dates for physicochemical data in tributary sites
# Since only one sample was collected per site per month,
# match the month with the sample collection date from sample_metadata
exact_date_lookup <- sample_metadata %>%
  filter(grepl("^T\\d+$", Site)) %>%
  mutate(
    Year = year(sample_collection_date),
    Month = month(sample_collection_date)
  ) %>%
  distinct(Site, Year, Month, sample_collection_date)

combined_pc_data_date_corrected <- combined_pc_data_raw %>%
  mutate(
    original_date = sample_collection_date,
    Year = lubridate::year(sample_collection_date),
    Month = lubridate::month(sample_collection_date)
  ) %>%
  left_join(
    exact_date_lookup %>%
      rename(exact_date = sample_collection_date),
    by = c("Site", "Year", "Month"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    sample_collection_date = if_else(
      grepl("^T\\d+$", Site) & !is.na(exact_date),
      exact_date,
      sample_collection_date
    )
  ) %>%
  select(-Year, -Month, -exact_date) %>%
  relocate(original_date)

# Physicochemical variables
pc_vars <- c(
  "BOD", "DO", "FC", "pH",
  "Ammonia", "Nitrate", "Phosphate"
)

#### BDL THRESHOLD HISTORY ####
# check data when a below detection (BDL) limit was first recorded
bdl_threshold_history <- combined_pc_data_raw %>%
  select(Year, sample_collection_date, Site, all_of(pc_vars)) %>%
  pivot_longer(
    cols = all_of(pc_vars),
    names_to = "Parameter",
    values_to = "Value"
  ) %>%
  mutate(Value = trimws(as.character(Value))) %>%
  filter(str_starts(Value, "<")) %>%
  group_by(Parameter, Value) %>%
  summarise(
    n = n(),
    first_date = min(sample_collection_date, na.rm = TRUE),
    first_year = min(Year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Parameter, first_date)

#### RECODE DATA ####
# Some recorded data used half-limit values.
# Will change those to the appropriate BDL thresholds

# function
recode_pc_value <- function(parameter, value, date) {
  
  x <- trimws(as.character(value))
  num <- suppressWarnings(as.numeric(x))
  
  case_when(
    
    # Ammonia: half of 0.01 before first reported BDL
    parameter == "Ammonia" &
      !is.na(num) & abs(num - 0.005) < 1e-8 &
      date < as.Date("2023-07-04") ~ "<0.01",
    
    # Ammonia: half of 0.002 before first reported BDL
    parameter == "Ammonia" &
      !is.na(num) & abs(num - 0.001) < 1e-8 &
      date < as.Date("2024-09-10") ~ "<0.002",
    
    # BOD
    parameter == "BOD" &
      !is.na(num) & abs(num - 0.5) < 1e-8 &
      date < as.Date("2023-12-01") ~ "<1",
    
    # Dissolved oxygen
    parameter == "DO" &
      !is.na(num) & abs(num - 0.05) < 1e-8 &
      date < as.Date("2023-07-01") ~ "<0.1",
    
    # Nitrate: half of 0.1 before first reported BDL
    parameter == "Nitrate" &
      !is.na(num) & abs(num - 0.05) < 1e-8 &
      date < as.Date("2023-07-01") ~ "<0.1",
    
    # Nitrate: half of 0.01 before first reported BDL
    parameter == "Nitrate" &
      !is.na(num) & abs(num - 0.005) < 1e-8 &
      date < as.Date("2024-04-01") ~ "<0.01",
    
    # Phosphate
    parameter == "Phosphate" &
      !is.na(num) & abs(num - 0.005) < 1e-8 &
      date < as.Date("2024-10-08") ~ "<0.01",
    
    # Retain all other values
    TRUE ~ x
  )
}

pc_recoding_audit <- combined_pc_data_raw %>%
  select(Year, sample_collection_date, Site, all_of(pc_vars)) %>%
  pivot_longer(
    cols = all_of(pc_vars),
    names_to = "Parameter",
    values_to = "Original_Value"
  ) %>%
  mutate(
    Original_Value = trimws(as.character(Original_Value)),
    Recoded_Value = recode_pc_value(
      Parameter,
      Original_Value,
      sample_collection_date
    ),
    Change = Original_Value != Recoded_Value,
    Reason = case_when(
      Change ~ paste0(
        "Half-detection-limit value recoded to ",
        Recoded_Value,
        " based on date-specific BDL reporting history"
      ),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Change) %>%
  select(
    Year,
    sample_collection_date,
    Site,
    Parameter,
    Original_Value,
    Recoded_Value,
    Reason
  ) %>%
  arrange(Parameter, sample_collection_date, Site)

# Summarize the number of changes by parameter and value
pc_recoding_summary <- pc_recoding_audit %>%
  count(
    Parameter,
    Original_Value,
    Recoded_Value,
    name = "n_changed"
  ) %>%
  arrange(Parameter, Original_Value)

# Recode the dataset
combined_pc_data_recoded <- combined_pc_data_date_corrected %>%
  mutate(
    across(
      all_of(pc_vars),
      ~ recode_pc_value(
        cur_column(),
        .x,
        sample_collection_date
      )
    )
  )

# Count BDL and missing values per parameter
pc_data_quality <- combined_pc_data_recoded %>%
  summarise(
    across(
      all_of(pc_vars),
      list(
        BDL = ~ sum(startsWith(trimws(.x), "<"), na.rm = TRUE),
        Missing = ~ sum(trimws(.x) == "***" | is.na(.x), na.rm = TRUE)
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Parameter", ".value"),
    names_pattern = "(.*)_(BDL|Missing)"
  ) %>%
  mutate(
    n_total = nrow(combined_pc_data_raw),
    n_quantified = n_total - BDL - Missing,
    quantified_pct = round(n_quantified / n_total * 100, 2),
    n_BDL = BDL,
    BDL_pct = round(BDL / n_total * 100, 2),
    n_missing = Missing,
    missing_pct = round(Missing / n_total * 100, 2)
  ) %>%
  select(-BDL, -Missing) %>%
  arrange(n_quantified)

#### CREATE SUPPLEMENTARY TABLES ####

# Finalize sample metadata
sample_metadata_complete <- sample_metadata %>%
  left_join(
    combined_pc_data_recoded %>%
      select(Site, sample_collection_date, all_of(pc_vars)),
    by = c("Site", "sample_collection_date")
  )

write_tsv(
  sample_metadata_complete,
  "sample_metadata.tsv"
)
write_tsv(
  bdl_threshold_history,
  "supplementary_files/Supplementary_BDL_threshold_history.tsv"
)
write_tsv(
  pc_recoding_audit,
  "supplementary_files/Supplementary_BDL_recoding_audit.tsv"
)
write_tsv(
  pc_recoding_summary,
  "supplementary_files/Supplementary_BDL_recoding_summary.tsv"
)
write_tsv(
  pc_data_quality,
  "supplementary_files/Supplementary_physicochemical_data_quality.tsv"
)
write_tsv(
  combined_pc_data_recoded,
  "physicochemical_data_recoded.tsv"
)


