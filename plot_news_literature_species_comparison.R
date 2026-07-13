if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(ggplot2))

species_levels <- c(
  "Asiatic Black Bear",
  "Grey Wolf",
  "Brown Bear",
  "Tiger",
  "Snow Leopard",
  "Leopard",
  "Eurasian Lynx",
  "Dhole",
  "Asiatic Golden Cat",
  "Wolverine",
  "Clouded Leopard"
)

# Raw counts from the pasted tables.
# Literature shares are based on 87 papers; one paper may involve multiple species,
# so the proportions can sum to > 1.
df <- data.frame(
  Species = c(
    "Wolverine", "Dhole", "Eurasian Lynx", "Leopard", "Tiger",
    "Asiatic Black Bear", "Snow Leopard", "Grey Wolf", "Brown Bear",
    "Asiatic Black Bear", "Grey Wolf", "Tiger", "Brown Bear", "Snow Leopard",
    "Leopard", "Eurasian Lynx", "Asiatic Golden Cat", "Clouded Leopard", "Dhole"
  ),
  Method = c(
    rep("Literature", 9),
    rep("News", 10)
  ),
  Count = c(
    1, 3, 7, 10, 14, 21, 26, 29, 30,
    187, 79, 78, 64, 44, 44, 16, 4, 2, 2
  ),
  stringsAsFactors = FALSE
)

# Add missing species as zero counts so the two methods share a common 11-species universe.
all_species <- expand.grid(
  Species = species_levels,
  Method = c("Literature", "News"),
  stringsAsFactors = FALSE
)
df <- merge(all_species, df, by = c("Species", "Method"), all.x = TRUE, sort = FALSE)
df$Count[is.na(df$Count)] <- 0

source_totals <- c(Literature = 87, News = 520)
df$Method <- factor(df$Method, levels = c("Literature", "News"))
df$Species <- factor(df$Species, levels = species_levels)
df$Proportion <- df$Count / source_totals[as.character(df$Method)] * 100

species_order <- with(
  aggregate(Count ~ Species, data = df, FUN = sum),
  as.character(Species[order(-Count, Species)])
)
df$Species <- factor(df$Species, levels = species_order)

contingency_all <- xtabs(Count ~ Method + Species, data = df)
chi_result <- chisq.test(contingency_all, correct = FALSE)

# Standardized residuals are used to derive species-level significance markers.
resid_p <- 2 * pnorm(-abs(chi_result$stdres["Literature", ]))

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

sig_df <- data.frame(
  Species = species_order,
  TestMethod = "Chi-square residual test (from 2x11 table)",
  PValue = as.numeric(resid_p[species_order]),
  SignifLabel = vapply(as.numeric(resid_p[species_order]), p_to_stars, character(1)),
  stringsAsFactors = FALSE
)

species_max <- aggregate(Proportion ~ Species, data = df, FUN = max)
sig_df <- merge(sig_df, species_max, by = "Species", sort = FALSE)
sig_df$Species <- factor(sig_df$Species, levels = species_order)
sig_df <- sig_df[order(sig_df$Species), ]
sig_df$YPos <- sig_df$Proportion + 3.0

y_upper <- max(sig_df$YPos) + 2

p <- ggplot(df, aes(x = Species, y = Proportion, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.70,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.1f", Proportion)),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 4.2,
    family = "Arial"
  ) +
  geom_text(
    data = sig_df,
    aes(x = Species, y = YPos, label = SignifLabel),
    inherit.aes = FALSE,
    family = "Arial",
    size = 5.8,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Literature" = "#1F77B4", "News" = "#E15759"),
    name = "Data Source"
  ) +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = seq(0, ceiling(y_upper / 10) * 10, by = 10),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = "Species",
    y = "Share of source units (%)"
  ) +
  theme_bw(base_size = 15, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

args <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- sub(file_flag, "", args[grep(file_flag, args)])

if (length(script_path) > 0) {
  output_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  output_dir <- getwd()
}

gg_file_png <- file.path(output_dir, "news_literature_species_comparison.png")
gg_file_pdf <- file.path(output_dir, "news_literature_species_comparison.pdf")
result_file_species <- file.path(output_dir, "species_significance_tests.csv")
result_file_global <- file.path(output_dir, "global_chi_square_test_species.csv")

# Use a modestly taller figure so the 11 species and significance labels stay readable.
ggsave(filename = gg_file_png, plot = p, width = 12, height = 7.2, dpi = 600)
ggsave(filename = gg_file_pdf, plot = p, width = 12, height = 7.2, device = cairo_pdf)

write.csv(
  sig_df[, c("Species", "TestMethod", "PValue", "SignifLabel")],
  result_file_species,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    TestMethod = "Chi-square test of independence (2x11)",
    Statistic = unname(chi_result$statistic),
    DF = unname(chi_result$parameter),
    PValue = chi_result$p.value
  ),
  result_file_global,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Figure saved to:")
message(" - ", gg_file_png)
message(" - ", gg_file_pdf)
message("Species-level table saved to:")
message(" - ", result_file_species)
message("Global test table saved to:")
message(" - ", result_file_global)
