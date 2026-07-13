# =========================================================
# Sankey diagram: Predator -> Victim type -> Outcome
# Data source: victim.csv (event-level counts)
# =========================================================

# If needed, install packages once:
# install.packages(c("dplyr", "ggplot2", "tidyr", "stringr"))
# install.packages("ggsankey")
# If ggsankey is unavailable on your CRAN mirror:
# install.packages("devtools")
# devtools::install_github("davidsjoberg/ggsankey")

required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr", "ggsankey")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      "\nPlease install them and rerun this script."
    )
  )
}

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggsankey)

input_csv <- "C:/Users/yuanj/Desktop/人兽冲突/victim_analysis/victim.csv"
out_dir <- "C:/Users/yuanj/Desktop/人兽冲突/victim_analysis"

normalize_victim <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "Na", "Nan")] <- NA_character_
  x <- str_to_title(x)
  x <- ifelse(x == "Dog", "Dog", x)
  x
}

collapse_victim_type <- function(x) {
  case_when(
    x %in% c("Chicken", "Duck", "Goose", "Poultry") ~ "Poultry",
    x %in% c("Cow", "Cattle") ~ "Cattle",
    x %in% c("Human", "Sheep", "Dog", "Cattle", "Poultry") ~ x,
    TRUE ~ "Other stock"
  )
}

predator_name_map <- c(
  "ABB" = "Black bear",
  "AGC" = "Golden cat",
  "BB" = "Brown bear",
  "EL" = "Leopard",
  "GW" = "Grey wolf",
  "L" = "Leopard",
  "LE" = "Leopard",
  "SL" = "Snow leopard",
  "T" = "Tiger",
  "TI" = "Tiger",
  "CL" = "Clouded leopard",
  "D" = "Dhole"
)

# ---------------------------------------------------------
# 1) Read and clean event-level data
# ---------------------------------------------------------
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("No.", "Victem", "Number of victems", "Number of deaths")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

df <- df %>%
  mutate(
    victim_raw = normalize_victim(Victem),
    n_victims = suppressWarnings(as.integer(round(as.numeric(`Number of victems`)))),
    n_deaths = suppressWarnings(as.integer(round(as.numeric(`Number of deaths`)))),
    n_victims = ifelse(is.na(n_victims), 0L, pmax(n_victims, 0L)),
    n_deaths = ifelse(is.na(n_deaths), 0L, pmax(n_deaths, 0L)),
    n_deaths = pmin(n_deaths, n_victims),
    predator_code = toupper(str_extract(as.character(`No.`), "^[A-Za-z]+")),
    predator_code = ifelse(is.na(predator_code) | predator_code == "", "UNK", predator_code),
    predator = ifelse(predator_code %in% names(predator_name_map), predator_name_map[predator_code], predator_code),
    victim_type = collapse_victim_type(victim_raw)
  ) %>%
  filter(!is.na(victim_raw), n_victims > 0)

if (nrow(df) == 0) {
  stop("No valid harmed-individual records found after cleaning.")
}

# Keep top predator groups for readability, merge the rest as "Other species".
predator_top_n <- 6
predator_rank <- df %>%
  group_by(predator) %>%
  summarise(total_harmed = sum(n_victims), .groups = "drop") %>%
  arrange(desc(total_harmed))

keep_predators <- head(predator_rank$predator, predator_top_n)
df <- df %>%
  mutate(predator = ifelse(predator %in% keep_predators, predator, "Other species"))

# ---------------------------------------------------------
# 2) Expand to individual-level rows (one row = one victim)
# ---------------------------------------------------------
df_individual <- bind_rows(
  df %>%
    transmute(predator, victim_type, outcome = "Died", n = n_deaths) %>%
    filter(n > 0),
  df %>%
    transmute(predator, victim_type, outcome = "Survived", n = n_victims - n_deaths) %>%
    filter(n > 0)
) %>%
  tidyr::uncount(weights = n)

if (nrow(df_individual) == 0) {
  stop("Expanded individual-level table is empty.")
}

# ---------------------------------------------------------
# 3) Ordering and labels
# ---------------------------------------------------------
predator_levels <- df_individual %>%
  count(predator, sort = TRUE) %>%
  pull(predator)

victim_base_order <- c("Human", "Sheep", "Cattle", "Other stock", "Dog", "Poultry")
victim_present <- unique(df_individual$victim_type)
victim_levels <- c(
  victim_base_order[victim_base_order %in% victim_present],
  setdiff(sort(victim_present), victim_base_order)
)

outcome_levels <- c("Died", "Survived")
node_levels <- c(predator_levels, victim_levels, outcome_levels)

outcome_stats <- df_individual %>%
  count(outcome) %>%
  mutate(pct = n / sum(n))

outcome_label_map <- setNames(
  sprintf("%s\n%.0f%% of total", outcome_stats$outcome, 100 * outcome_stats$pct),
  outcome_stats$outcome
)

# ---------------------------------------------------------
# 4) Build sankey data and plot
# ---------------------------------------------------------
df_long <- df_individual %>%
  make_long(predator, victim_type, outcome) %>%
  mutate(
    node = as.character(node),
    next_node = as.character(next_node),
    node = factor(node, levels = node_levels),
    next_node = factor(next_node, levels = node_levels),
    fill_group = case_when(
      as.character(node) == "Died" ~ "Died",
      as.character(node) == "Survived" ~ "Survived",
      as.character(node) %in% victim_levels ~ "Victim",
      TRUE ~ "Predator"
    ),
    node_label = ifelse(
      as.character(node) %in% names(outcome_label_map),
      outcome_label_map[as.character(node)],
      as.character(node)
    )
  )

fill_colors <- c(
  "Predator" = "#9CB68A",
  "Victim" = "#E3D7BD",
  "Died" = "#D5A2A2",
  "Survived" = "#8FBFBC"
)

p <- ggplot(
  df_long,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = fill_group,
    label = node_label
  )
) +
  scale_x_discrete(labels = c("predator" = "Predator", "victim_type" = "Victim", "outcome" = "Outcome")) +
  geom_sankey(
    flow.alpha = 0.55,
    node.color = "#8B8B8B",
    width = 0.12,
    smooth = 8
  ) +
  geom_sankey_label(
    size = 4,
    color = "#2E2E2E",
    fill = alpha("white", 0.85),
    family = "sans"
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_sankey(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(size = 16, face = "bold", color = "#2E2E2E"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p)

png_path <- file.path(out_dir, "predator_victim_sankey_from_csv.png")
pdf_path <- file.path(out_dir, "predator_victim_sankey_from_csv.pdf")
ggsave(png_path, plot = p, width = 12, height = 8, dpi = 320, bg = "white")
ggsave(pdf_path, plot = p, width = 12, height = 8, bg = "white")

cat("Saved files:\n")
cat(" - ", png_path, "\n", sep = "")
cat(" - ", pdf_path, "\n", sep = "")

summary_tbl <- df_individual %>%
  count(outcome) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("\nOutcome summary (individual-level):\n")
print(summary_tbl)
