### Which metrics contribute to introduction

# --- Packages ----------------------------------------------------------------
library(tidyverse)   # dplyr, ggplot2, readr, forcats, etc.
library(scales)      # rescale()

# --- Paths -------------------------------------------------------------------
in_path  <- "C:/Users/Owner/OneDrive/ILP.csv"
out_dir <- "C:/Users/Owner/OneDrive/Graphs"
out_csv  <- file.path(out_dir, "ILP_results.csv")
out_pngA <- file.path(out_dir, "ILP_panelA_ILP_bar.png")
out_pngB <- file.path(out_dir, "ILP_panelB_contributions.png")

# --- Read CSV (force Windows encoding to avoid invalid UTF-8) ----------------
df_raw <- readr::read_csv(
  in_path,
  show_col_types = FALSE,
  locale = readr::locale(encoding = "Windows-1252")
)

# --- Clean column names to safe snake_case -----------------------------------
# (replace spaces with underscores, collapse repeats, trim)
clean_headers <- function(v){
  v <- gsub("\\s+", "_", v)
  v <- gsub("[^A-Za-z0-9_]", "_", v)
  v <- gsub("_+", "_", v)
  v <- gsub("^_|_$", "", v)
  v
}
names(df_raw) <- clean_headers(names(df_raw))

# Expecting at least these columns after cleaning:
# Species, pi_native, pi_invasive, HO_native, HO_invasive,
# fROH_native, fROH_invasive, LongROH_native, LongROH_invasive

# --- Clean species strings (normalize to UTF-8, fix NBSP) --------------------
clean_species <- function(x){
  x <- iconv(x, from = "", to = "UTF-8", sub = " ")
  x <- gsub("\u00A0", " ", x, perl = TRUE)  # non-breaking space -> space
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}
df_raw$Species <- clean_species(df_raw$Species)

# --- Sanity check: make sure required columns exist --------------------------
required_cols <- c("Species","pi_native","pi_invasive","HO_native","HO_invasive",
                   "fROH_native","fROH_invasive","LongROH_invasive")
missing <- setdiff(required_cols, names(df_raw))
if (length(missing)){
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# --- Compute deltas & keep invasive long-ROH proportion ----------------------
df <- df_raw %>%
  mutate(
    delta_pi   = pi_invasive - pi_native,
    delta_HO   = HO_invasive - HO_native,       # (percentage points)
    delta_fROH = fROH_invasive - fROH_native,
    longROH_prop_invasive = LongROH_invasive
  )

# --- z-score helper that won’t blow up on constant vectors -------------------
z <- function(x){
  sx <- sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# --- ILP calculation, ranks, normalized 0-1 ----------------------------------
df_ilp <- df %>%
  mutate(
    z_dpi  = z(delta_pi),
    z_dho  = z(delta_HO),
    z_dfh  = z(delta_fROH),
    z_lroh = z(longROH_prop_invasive),
    ILP    = z_dpi + z_dho - z_dfh - z_lroh,
    ILP_rank   = rank(-ILP, ties.method = "min"),
    ILP_norm01 = rescale(ILP, to = c(0,1), from = range(ILP, na.rm = TRUE))
  ) %>%
  arrange(desc(ILP))

# --- Save table & view in RStudio -------------------------------------------
readr::write_csv(df_ilp, out_csv)
print(df_ilp, n = nrow(df_ilp))
if (interactive()) View(df_ilp)

# --- Plot A: ILP (normalized 0–1) by species --------------------------------
df_ilp <- df_ilp |> mutate(Species = fct_reorder(Species, ILP_norm01))

pA <- ggplot(df_ilp, aes(Species, ILP_norm01)) +
  # outline-only bars
  geom_col(fill = NA, color = "goldenrod3", linewidth = 1, width = 0.75) +
  # value labels
  geom_text(aes(label = number(ILP_norm01, accuracy = 0.01)),
            hjust = -0.12, size = 3) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "A. Introduction-Likelihood Proxy (ILP)",
       x = NULL, y = "ILP (normalized 0–1)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),   # remove all gridlines
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave(out_pngA, pA, width = 8, height = 6.5, dpi = 300)
pA  # show in viewer

# --- Plot B: Signed contributions to ILP for each species --------------------
contrib_long <- df_ilp %>%
  transmute(Species,
            `+ z(Δπ)`       = z_dpi,
            `+ z(ΔHO)`      = z_dho,
            `− z(ΔfROH)`    = -z_dfh,     # minus to match ILP formula
            `− z(Long ROH)` = -z_lroh) %>%
  pivot_longer(-Species, names_to = "Component", values_to = "Contribution") %>%
  mutate(
    Species   = fct_reorder(Species, df_ilp$ILP[match(Species, df_ilp$Species)]),
    Component = factor(Component, levels = c("+ z(Δπ)", "+ z(ΔHO)", "− z(ΔfROH)", "− z(Long ROH)"))
  )

pB <- ggplot(contrib_long, aes(Species, Contribution, fill = Component)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
  coord_flip() +
  labs(title = "B. Metric contributions to ILP",
       x = NULL, y = "Signed z-score contribution") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank()) +     # <-- removes all gridlines
  guides(fill = guide_legend(title = NULL)) # keep legend title off

ggsave(out_pngB, pB, width = 8, height = 6.5, dpi = 300)
pB  # show in viewer
       
       
       
       
       

