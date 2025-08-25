# =========================================================
# NUCLEOTIDE DIVERSITY (π) — CLEANED SCRIPT (LOGIC UNCHANGED)
# =========================================================

# Libraries
# ---------------------------
library(dplyr)
library(tidyr)
library(readr)
library(effsize)
library(rstatix)
library(ggplot2)
library(gghalves)
library(scales)

# -------------------------------------
# Config: base folder and species map
# -------------------------------------
file_dir <- "~/SCHOOL/PhD/CLASSES/DISSERTATION/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3/Nucleotide"

species_list <- c(
  "Asian tiger mosquito"    = "Mosq",
  "Longhorned tick"         = "Ticks",
  "Goldfish"                = "Fish",
  "Late blight fungus"      = "Fungus",
  "Two-spotted spider mite" = "Mite",
  "Chinese mitten crab"     = "Crab",
  "Wild boar"               = "Pig",
  "Colorado potato beetle"  = "Beetle",
  "Rainbow trout"           = "RainbowFish",
  "Old World bollworm"      = "Old_World_Bollworm",
  "Silver carp"             = "SilverCarp"
)

# =========================================================
# 1) Read & combine windowed π data (all species)
#    - Adds BIN_SIZE (for weighted means and sanity checks)
# =========================================================
all_pi_data <- data.frame()
for (sp in names(species_list)) {
  code <- species_list[[sp]]
  native_file   <- file.path(file_dir, paste0("N_",  code, ".filtered.variants.10000.pi.windowed.pi"))
  invasive_file <- file.path(file_dir, paste0("IN_", code, ".filtered.variants.10000.pi.windowed.pi"))
  
  native_pi <- if (file.exists(native_file)) {
    df <- read.table(native_file, header = TRUE)
    df$Status <- "Native"; df$Species <- sp; df
  } else NULL
  
  invasive_pi <- if (file.exists(invasive_file)) {
    df <- read.table(invasive_file, header = TRUE)
    df$Status <- "Invasive"; df$Species <- sp; df
  } else NULL
  
  all_pi_data <- bind_rows(all_pi_data, native_pi, invasive_pi)
}

all_pi_data <- all_pi_data %>%
  mutate(BIN_SIZE = BIN_END - BIN_START + 1)

# =========================================================
# 2) Pooled (ALL WINDOWS) summaries & tests
#    - Wilcoxon on ALL windows
#    - Cliff’s Δ on down-sampled 5,000 per status (as you did)
# =========================================================
pi_summary <- all_pi_data %>%
  group_by(Status) %>%
  summarise(
    Mean_Pi          = mean(PI, na.rm = TRUE),
    Weighted_Mean_Pi = weighted.mean(PI, w = BIN_SIZE, na.rm = TRUE),
    Median_Pi        = median(PI, na.rm = TRUE),
    SD_Pi            = sd(PI, na.rm = TRUE),
    N                = n(),
    .groups = "drop"
  )
print(pi_summary)

# Wilcoxon (all windows)
pooled_wilcox <- wilcox.test(PI ~ Status, data = all_pi_data)
print(pooled_wilcox)

# Cliff’s Δ on down-sampled 5,000 per status
set.seed(42)
native_sub   <- all_pi_data %>% filter(Status == "Native")   %>% sample_n(5000)
invasive_sub <- all_pi_data %>% filter(Status == "Invasive") %>% sample_n(5000)
pooled_cliff <- cliff.delta(PI ~ Status, data = bind_rows(native_sub, invasive_sub))
print(pooled_cliff)

# =========================================================
# 3) Per-species summaries + Δπ + species-level Cliff’s Δ
#    - Weighted means per species & status
#    - Down-sampled Cliff’s Δ per species (5k/status)
# =========================================================
species_pi_summary <- all_pi_data %>%
  group_by(Species, Status) %>%
  summarise(
    Mean_Pi          = mean(PI, na.rm = TRUE),
    Weighted_Mean_Pi = weighted.mean(PI, w = BIN_SIZE, na.rm = TRUE),
    Median_Pi        = median(PI, na.rm = TRUE),
    SD_Pi            = sd(PI, na.rm = TRUE),
    N                = n(),
    .groups = "drop"
  )
print(species_pi_summary)

effect_sizes <- data.frame()
for (sp in unique(all_pi_data$Species)) {
  sp_data <- all_pi_data %>% filter(Species == sp)
  if (length(unique(sp_data$Status)) == 2) {
    native   <- sp_data %>% filter(Status == "Native")   %>% sample_n(min(5000, n()), replace = FALSE)
    invasive <- sp_data %>% filter(Status == "Invasive") %>% sample_n(min(5000, n()), replace = FALSE)
    combined <- bind_rows(native, invasive)
    
    cd <- cliff.delta(PI ~ Status, data = combined)
    
    effect_sizes <- bind_rows(effect_sizes, data.frame(
      Species         = sp,
      Cliff_Delta     = cd$estimate,
      Delta_Magnitude = cd$magnitude
    ))
  }
}

species_pi_summary <- left_join(species_pi_summary, effect_sizes, by = "Species")

pi_wide <- species_pi_summary %>%
  select(Species, Status, Weighted_Mean_Pi) %>%
  pivot_wider(names_from = Status, values_from = Weighted_Mean_Pi) %>%
  mutate(Delta_Pi = Invasive - Native)

# (Removed duplicate heatmap stats block: it recomputed per-species p-values and merged_stats but
#  wasn’t used anywhere else. Keeping Δπ & Cliff’s Δ above — the stats used.)

# =========================================================
# 4) Figure 2 — Half-violin with points, box, and mean dot
#    - 5k/status
# =========================================================
set.seed(123)
plot_data <- all_pi_data %>%
  group_by(Status) %>%
  slice_sample(n = 5000) %>%
  ungroup()

ggplot(plot_data, aes(x = Status, y = PI, fill = Status, color = Status)) +
  gghalves::geom_half_violin(side = "l", alpha = 0.8, trim = TRUE, width = 0.9) +
  gghalves::geom_half_point(
    side = "r", shape = 21, size = 1.2, stroke = 0.3, alpha = 0.3, width = 0.2
  ) +
  geom_boxplot(width = 0.1, linewidth = 1, outlier.shape = NA, fill = "white", color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = c("Native" = "#458B00", "Invasive" = "goldenrod3")) +
  scale_color_manual(values = c("Native" = "#458B00", "Invasive" = "goldenrod3")) +
  labs(x = "", y = "Nucleotide diversity (π)", title = NULL) +
  coord_cartesian(ylim = c(0, 0.006)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold")
  )

ggsave(
  filename = "Nucleotide_Diversity_HalfViolinBigger.png",
  path = "C:/Users/Owner/OneDrive - Auburn University/Documents/SCHOOL/PhD/CLASSES/DISSERTATION/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3/Graphs",
  dpi = 600, width = 10, height = 7, units = "in"
)

# =========================================================
# 5) Quick sanity checks (counts and bin sizes)
#    - Kept as diagnostics; remove if not needed
# =========================================================
all_pi_data %>% group_by(Species, Status) %>% summarise(N_windows = n(), .groups = "drop")
all_pi_data %>% group_by(Species) %>% summarise(
  Mean_Bin_Size = mean(BIN_SIZE), Median_Bin_Size = median(BIN_SIZE),
  Min_Bin_Size = min(BIN_SIZE),   Max_Bin_Size = max(BIN_SIZE),
  .groups = "drop"
)
all_pi_data %>% group_by(Species, Status) %>% summarise(
  Mean_Bin_Size = mean(BIN_SIZE), Median_Bin_Size = median(BIN_SIZE),
  Min_Bin_Size = min(BIN_SIZE),   Max_Bin_Size = max(BIN_SIZE),
  .groups = "drop"
)

# =========================================================
# 6) Regional analyses (6 species): telomere vs centromere
#    - Same logic; I only renamed the mapping to avoid
#      overwriting main species_list.
# =========================================================
setwd("C:/Users/Owner/OneDrive - Auburn University/Documents/SCHOOL/PhD/CLASSES/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3/Nucleotide")

region_species <- list(
  Crab = "Chinese mitten crab",
  Mosq = "Asian tiger mosquito",
  Fish = "Goldfish",
  RainbowFish = "Rainbow trout",
  Old_World_Bollworm = "Old World Bollworm",
  SilverCarp = "Silver carp"
)

# 6a) Build per-species combined data & per-panel means 
all_combined_data <- data.frame()
region_means <- data.frame()

for (species_code in names(region_species)) {
  
  native <- read_tsv(paste0("N_",  species_code, ".filtered.variants.10000.pi.windowed.pi"), show_col_types = FALSE) %>%
    mutate(Status = "Native")
  invasive <- read_tsv(paste0("IN_", species_code, ".filtered.variants.10000.pi.windowed.pi"), show_col_types = FALSE) %>%
    mutate(Status = "Invasive")
  
  combined <- bind_rows(native, invasive)
  chrom_lengths <- combined %>%
    group_by(CHROM) %>%
    summarize(chr_len = max(BIN_END), .groups = "drop")
  
  combined <- combined %>%
    left_join(chrom_lengths, by = "CHROM") %>%
    mutate(
      REL_POS = (BIN_START + BIN_END) / 2 / chr_len,
      Species = region_species[[species_code]],
      Region  = case_when(
        REL_POS <= 0.05 ~ "Telomere",
        REL_POS >= 0.95 ~ "Telomere",
        REL_POS >= 0.45 & REL_POS <= 0.55 ~ "Centromere",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Region))
  
  all_combined_data <- bind_rows(all_combined_data, combined)
  
  summary_df <- combined %>%
    group_by(Species, Status, Region) %>%
    summarize(
      Mean_pi = mean(PI, na.rm = TRUE),
      SD_pi   = sd(PI, na.rm = TRUE),
      N       = n(),
      SE_pi   = SD_pi / sqrt(N),
      .groups = "drop"
    )
  
  region_means <- bind_rows(region_means, summary_df)
}

# 6b) Figure: per-species regional bars with SE
region_pi_plot <- ggplot(region_means, aes(x = Region, y = Mean_pi, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Mean_pi - SE_pi, ymax = Mean_pi + SE_pi),
                position = position_dodge(width = 0.7), width = 0.2) +
  facet_wrap(~Species, scales = "free_y") +
  labs(y = expression("Mean " * pi), x = "Region", fill = "Status") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y  = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "black", linewidth = 0.6)
  ) +
  scale_fill_manual(values = c("Native" = "#458B00", "Invasive" = "goldenrod3"))

ggsave(
  filename = "C:/Users/Owner/OneDrive - Auburn University/Documents/SCHOOL/PhD/CLASSES/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3/Graphs/PerSpecies_Region_Status_Comparison.png",
  plot = region_pi_plot, width = 12, height = 7, dpi = 600
)

# 6c) Per-species/status regional Wilcoxon tests + r 
stats_results  <- data.frame()
effect_sizes_r <- data.frame()
group_means    <- data.frame()

for (sp in unique(all_combined_data$Species)) {
  for (status in c("Native", "Invasive")) {
    subset_df <- all_combined_data %>%
      filter(Species == sp, Status == status, Region %in% c("Telomere", "Centromere"))
    
    if (length(unique(subset_df$Region)) == 2) {
      test <- wilcox_test(subset_df, PI ~ Region, paired = FALSE)
      test$Species <- sp; test$Status <- status
      stats_results <- bind_rows(stats_results, test)
      
      eff <- wilcox_effsize(subset_df, PI ~ Region, ci = TRUE)  # r effect size
      eff$Species <- sp; eff$Status <- status
      effect_sizes_r <- bind_rows(effect_sizes_r, eff)
      
      means <- subset_df %>%
        group_by(Region) %>%
        summarise(
          Mean = mean(PI, na.rm = TRUE), SD = sd(PI, na.rm = TRUE), N = n(),
          SE = SD / sqrt(N),
          CI_low  = Mean - qt(0.975, df = N - 1) * SE,
          CI_high = Mean + qt(0.975, df = N - 1) * SE,
          .groups = "drop"
        ) %>% mutate(Species = sp, Status = status)
      
      group_means <- bind_rows(group_means, means)
    }
  }
}

print(stats_results)
print(effect_sizes_r)
print(group_means)

write_csv(stats_results,  "Wilcoxon_Telomere_vs_Centromere_by_Species.csv")
write_csv(effect_sizes_r, "Effect_Sizes_Telomere_vs_Centromere.csv")
write_csv(group_means,    "GroupMeans_Telomere_vs_Centromere.csv")

# 6d) Within-region (telomere or centromere): Native vs Invasive per species 
comparison_results <- data.frame()
region_means2 <- data.frame()

for (species_code in names(region_species)) {
  native <- read_tsv(paste0("N_",  species_code, ".filtered.variants.10000.pi.windowed.pi"), show_col_types = FALSE) %>%
    mutate(Status = "Native")
  invasive <- read_tsv(paste0("IN_", species_code, ".filtered.variants.10000.pi.windowed.pi"), show_col_types = FALSE) %>%
    mutate(Status = "Invasive")
  
  combined <- bind_rows(native, invasive)
  
  chrom_lengths <- combined %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(BIN_END), .groups = "drop")
  
  combined <- combined %>%
    left_join(chrom_lengths, by = "CHROM") %>%
    mutate(
      REL_POS = (BIN_START + BIN_END) / 2 / chr_len,
      Species = region_species[[species_code]],
      Region  = case_when(
        REL_POS <= 0.05 ~ "Telomere",
        REL_POS >= 0.95 ~ "Telomere",
        REL_POS >= 0.45 & REL_POS <= 0.55 ~ "Centromere",
        TRUE ~ NA_character_
      )
    ) %>% filter(!is.na(Region))
  
  for (reg in c("Telomere", "Centromere")) {
    reg_data <- combined %>% filter(Region == reg)
    if (length(unique(reg_data$Status)) == 2) {
      test <- wilcox_test(reg_data, PI ~ Status)
      test$Species <- region_species[[species_code]]
      test$Region  <- reg
      comparison_results <- bind_rows(comparison_results, test)
      
      reg_means <- reg_data %>%
        group_by(Status) %>%
        summarise(
          Mean = mean(PI, na.rm = TRUE), SD = sd(PI, na.rm = TRUE), N = n(),
          SE = SD / sqrt(N),
          CI_low  = Mean - qt(0.975, df = N - 1) * SE,
          CI_high = Mean + qt(0.975, df = N - 1) * SE,
          .groups = "drop"
        ) %>%
        mutate(Species = region_species[[species_code]], Region = reg)
      
      region_means2 <- bind_rows(region_means2, reg_means)
    }
  }
}

print(comparison_results)
print(region_means2)

write_csv(comparison_results, "status_comparison_within_region_by_species.csv")
write_csv(region_means2,      "GroupMeans_Status_by_Region.csv")

#  6e) Pooled regional contrasts (across species) 
pooled_data <- all_combined_data %>% filter(Region %in% c("Telomere", "Centromere"))

get_delta <- function(df, group_col) {
  means <- df %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean_pi = mean(PI, na.rm = TRUE)) %>%
    pivot_wider(names_from = !!sym(group_col), values_from = mean_pi)
  delta <- means[[2]] - means[[1]]
  names(delta) <- NULL
  delta
}

get_means <- function(df, group_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean_pi = mean(PI, na.rm = TRUE)) %>%
    arrange(.data[[group_col]]) %>%
    pull(mean_pi)
}

telomere_data   <- pooled_data %>% filter(Region == "Telomere")
telomere_means  <- get_means(telomere_data, "Status")
telomere_test   <- wilcox_test(telomere_data, PI ~ Status) %>%
  mutate(
    Test    = "Telomere: Native vs Invasive",
    Mean_1  = telomere_means[1],
    Mean_2  = telomere_means[2],
    Delta   = get_delta(telomere_data, "Status"),
    EffectSize = cliff.delta(PI ~ Status, data = telomere_data)$estimate,
    CI_low     = cliff.delta(PI ~ Status, data = telomere_data)$conf.int[1],
    CI_high    = cliff.delta(PI ~ Status, data = telomere_data)$conf.int[2],
    Magnitude  = cliff.delta(PI ~ Status, data = telomere_data)$magnitude
  )

centromere_data  <- pooled_data %>% filter(Region == "Centromere")
centromere_means <- get_means(centromere_data, "Status")
centromere_test  <- wilcox_test(centromere_data, PI ~ Status) %>%
  mutate(
    Test    = "Centromere: Native vs Invasive",
    Mean_1  = centromere_means[1],
    Mean_2  = centromere_means[2],
    Delta   = get_delta(centromere_data, "Status"),
    EffectSize = cliff.delta(PI ~ Status, data = centromere_data)$estimate,
    CI_low     = cliff.delta(PI ~ Status, data = centromere_data)$conf.int[1],
    CI_high    = cliff.delta(PI ~ Status, data = centromere_data)$conf.int[2],
    Magnitude  = cliff.delta(PI ~ Status, data = centromere_data)$magnitude
  )

native_data  <- pooled_data %>% filter(Status == "Native")
native_means <- get_means(native_data, "Region")
native_region_test <- wilcox_test(native_data, PI ~ Region) %>%
  mutate(
    Test    = "Native: Telomere vs Centromere",
    Mean_1  = native_means[1],
    Mean_2  = native_means[2],
    Delta   = get_delta(native_data, "Region"),
    EffectSize = cliff.delta(PI ~ Region, data = native_data)$estimate,
    CI_low     = cliff.delta(PI ~ Region, data = native_data)$conf.int[1],
    CI_high    = cliff.delta(PI ~ Region, data = native_data)$conf.int[2],
    Magnitude  = cliff.delta(PI ~ Region, data = native_data)$magnitude
  )

invasive_data  <- pooled_data %>% filter(Status == "Invasive")
invasive_means <- get_means(invasive_data, "Region")
invasive_region_test <- wilcox_test(invasive_data, PI ~ Region) %>%
  mutate(
    Test    = "Invasive: Telomere vs Centromere",
    Mean_1  = invasive_means[1],
    Mean_2  = invasive_means[2],
    Delta   = get_delta(invasive_data, "Region"),
    EffectSize = cliff.delta(PI ~ Region, data = invasive_data)$estimate,
    CI_low     = cliff.delta(PI ~ Region, data = invasive_data)$conf.int[1],
    CI_high    = cliff.delta(PI ~ Region, data = invasive_data)$conf.int[2],
    Magnitude  = cliff.delta(PI ~ Region, data = invasive_data)$magnitude
  )

overall_results <- bind_rows(
  native_region_test,
  invasive_region_test,
  telomere_test,
  centromere_test
) %>%
  select(Test, Mean_1, Mean_2, Delta, EffectSize, CI_low, CI_high, Magnitude, statistic, p)

print(overall_results)
write_csv(overall_results, "Region_Status_Comparisons_with_Means_and_CI.csv")

# 6f) Optional pooled bar plot
summary_data <- pooled_data %>%
  filter(Region %in% c("Telomere", "Centromere")) %>%
  mutate(
    Comparison = case_when(
      Region == "Telomere"   & Status %in% c("Native", "Invasive") ~ "Telomere: Native vs Invasive",
      Region == "Centromere" & Status %in% c("Native", "Invasive") ~ "Centromere: Native vs Invasive",
      Status == "Native"   & Region %in% c("Telomere", "Centromere") ~ "Native: Telomere vs Centromere",
      Status == "Invasive" & Region %in% c("Telomere", "Centromere") ~ "Invasive: Telomere vs Centromere",
      TRUE ~ NA_character_
    ),
    FacetLabel = case_when(
      Comparison %in% c("Telomere: Native vs Invasive", "Centromere: Native vs Invasive") ~ Region,
      Comparison %in% c("Native: Telomere vs Centromere", "Invasive: Telomere vs Centromere") ~ Status,
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(FacetLabel, Region, Status) %>%
  summarise(
    Mean = mean(PI, na.rm = TRUE),
    SD   = sd(PI, na.rm = TRUE),
    N    = n(),
    SE   = SD / sqrt(N),
    .groups = "drop"
  )

colors <- c("Native" = "#458B00", "Invasive" = "goldenrod3")

ggplot(summary_data, aes(x = Region, y = Mean, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(x = "Region", y = expression(bold("Mean " * pi)), fill = "Status") +
  theme_minimal(base_size = 14) +
  theme(
    axis.line   = element_line(color = "black", linewidth = 0.6),
    axis.ticks  = element_line(color = "black"),
    axis.text   = element_text(face = "bold", color = "black"),
    axis.title  = element_text(face = "bold", color = "black"),
    legend.title= element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.grid  = element_blank(),
    panel.border= element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = "C:/Users/Owner/OneDrive - Auburn University/Documents/SCHOOL/PhD/CLASSES/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3/Graphs/Region_Status_Comparison_Plot.png",
  plot = last_plot(), width = 10, height = 6, dpi = 600
)
