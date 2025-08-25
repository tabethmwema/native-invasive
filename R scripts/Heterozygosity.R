# =========================================================
# OBSERVED HETEROZYGOSITY (ΔHO)
# =========================================================

# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(car)       # Levene's test
library(effsize)   # Cliff's delta


setwd("~/SCHOOL/PhD/CLASSES/DISSERTATION/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3")

# Load data
In_N_data <- read.csv("Merged_Heterozygosity_Standardized.csv", header = TRUE)

# Keep factor order so Cliff's delta sign matches text:
# (positive δ means Invasive > Native)
In_N_data$Status <- factor(In_N_data$Status, levels = c("Invasive", "Native"))

# Treat INDV as the common name label for plotting/paired comparison
In_N_data <- In_N_data %>%
  mutate(Species = str_trim(INDV))

# =========================================================
# (A) Assumption checks (for context; stats below use nonparametrics)
# =========================================================
# 1) Shapiro–Wilk on the raw observed heterozygosity values
shapiro_result <- shapiro.test(In_N_data$Observed_Het.Percent)

# 2) Levene’s test for equal variances across Status
levene_result <- car::leveneTest(Observed_Het.Percent ~ Status, data = In_N_data)

print(shapiro_result)
print(levene_result)

# =========================================================
# (B) Paired comparison by species (used for Figure 4 bars)
#     ΔHO = HO_invasive – HO_native per species
# =========================================================
het_pairs <- In_N_data %>%
  select(Species, Status, Observed_Het.Percent) %>%
  pivot_wider(names_from = Status, values_from = Observed_Het.Percent) %>%
  drop_na() %>%
  mutate(Diff = Invasive - Native)

# Paired Wilcoxon (signed-rank) across species-level pairs
paired_result <- wilcox.test(het_pairs$Invasive, het_pairs$Native, paired = TRUE)
print(paired_result)

# =========================================================
# (C) Figure 4 — per-species ΔHO bars (percentage points)
# =========================================================
het_plot <- ggplot(het_pairs, aes(x = reorder(Species, Diff), y = Diff)) +
  geom_bar(
    stat  = "identity",
    fill  = ifelse(het_pairs$Diff < 0, "plum3", "orchid4"),
    color = "black"
  ) +
  geom_text(
    aes(label = round(Diff, 2)),
    hjust = ifelse(het_pairs$Diff > 0, -0.2, 1.1),
    color = "black", size = 4, fontface = "bold"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Species",
    y = expression(Delta ~ "Observed Heterozygosity (%)")
  ) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(
    axis.text  = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.ticks = element_line(color = "black"),
    axis.line  = element_line(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
print(het_plot)

# Save figure (same path & settings)
ggsave(
  filename = "C:/Users/Owner/OneDrive - Auburn University/Documents/SCHOOL/PhD/CLASSES/DISSERTATION/DISSERTATION CHAPTERS/In use/NEW CHAPTERS/VCFs/New_use this/25_march_2025/All 3/Graphs/Delta_Heterozygosity_High.png",
  plot = het_plot, width = 13, height = 8, dpi = 600, units = "in"
)

# =========================================================
# (D) Unpaired pooled comparison across all samples
#     Wilcoxon rank-sum + group summaries + Cliff’s δ
# =========================================================
# Wilcoxon rank-sum
wilcox_result <- wilcox.test(Observed_Het.Percent ~ Status, data = In_N_data)
print(wilcox_result)

# Group-wise mean/median (the numbers you quoted)
summary_stats <- In_N_data %>%
  group_by(Status) %>%
  summarise(
    Mean_Het   = mean(Observed_Het.Percent, na.rm = TRUE),
    Median_Het = median(Observed_Het.Percent, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# Effect size: Cliff’s delta (positive = Invasive > Native)
cd <- effsize::cliff.delta(Observed_Het.Percent ~ Status, data = In_N_data)
print(cd)

# printout of headline values
cat("\n--- HEADLINE SUMMARY ---\n")
print(summary_stats %>% mutate(across(where(is.numeric), ~ round(.x, 1))))
cat("Wilcoxon rank-sum: W =", wilcox_result$statistic, " p =", signif(wilcox_result$p.value, 3), "\n")
cat("Cliff's delta (Invasive vs Native): δ =", round(unname(cd$estimate), 2),
    " (", cd$magnitude, ")\n", sep = "")
cat("Shapiro-Wilk: W =", round(unname(shapiro_result$statistic), 2),
    ", p <", ifelse(shapiro_result$p.value < 0.001, "0.001", round(shapiro_result$p.value, 3)), "\n", sep = "")
cat("Levene's: F =", round(levene_result$`F value`[1], 2),
    ", p =", signif(levene_result$`Pr(>F)`[1], 3), "\n", sep = "")
