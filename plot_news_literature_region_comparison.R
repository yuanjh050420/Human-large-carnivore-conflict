if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(ggplot2))

regions <- c(
  "North China",
  "Northeast China",
  "East China",
  "Central China",
  "South China",
  "Southwest China",
  "Northwest China"
)

df <- data.frame(
  Region = rep(regions, times = 2),
  Method = rep(c("Literature", "News"), each = length(regions)),
  Frequency = c(
    5, 12, 2, 1, 0, 21, 37,
    19, 197, 38, 20, 2, 98, 146
  ),
  stringsAsFactors = FALSE
)

df$Method <- factor(df$Method, levels = c("Literature", "News"))
df$Proportion <- ave(df$Frequency, df$Method, FUN = function(x) x / sum(x) * 100)

literature_order <- df[df$Method == "Literature", c("Region", "Proportion")]
regions <- literature_order$Region[order(literature_order$Proportion, decreasing = TRUE)]
df$Region <- factor(df$Region, levels = regions)

contingency_all <- xtabs(Frequency ~ Method + Region, data = df)
chi_result <- chisq.test(contingency_all, correct = FALSE)

# Region-level markers are derived from chi-square standardized residuals.
region_p <- 2 * pnorm(-abs(chi_result$stdres["Literature", ]))

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

sig_df <- data.frame(
  Region = regions,
  TestMethod = "Chi-square residual test (from 2x7 table)",
  PValue = as.numeric(region_p),
  SignifLabel = NA_character_,
  stringsAsFactors = FALSE
)

sig_df$SignifLabel <- vapply(sig_df$PValue, p_to_stars, character(1))

region_max <- aggregate(Proportion ~ Region, data = df, FUN = max)
sig_df <- merge(sig_df, region_max, by = "Region", sort = FALSE)
sig_df$Region <- factor(sig_df$Region, levels = regions)
sig_df <- sig_df[order(sig_df$Region), ]
sig_df$YPos <- sig_df$Proportion + 2.2

y_upper <- max(sig_df$YPos) + 1

p <- ggplot(df, aes(x = Region, y = Proportion, fill = Method)) +
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
    size = 3.6,
    family = "Arial"
  ) +
  scale_fill_manual(
    values = c("Literature" = "#1F77B4", "News" = "#E15759"),
    name = "Data Source"
  ) +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = seq(0, 55, by = 5),
    expand = expansion(mult = c(0, 0.04))
  ) +
  geom_text(
    data = sig_df,
    aes(x = Region, y = YPos, label = SignifLabel),
    inherit.aes = FALSE,
    family = "Arial",
    size = 5,
    fontface = "bold"
  ) +
  labs(
    x = "Geographic Region",
    y = "Proportion (%)"
  ) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
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

png_file <- file.path(output_dir, "news_literature_region_comparison.png")
pdf_file <- file.path(output_dir, "news_literature_region_comparison.pdf")
pvalue_file <- file.path(output_dir, "regional_significance_tests.csv")
global_file <- file.path(output_dir, "global_chi_square_test.csv")

ggsave(filename = png_file, plot = p, width = 11, height = 6.5, dpi = 600)
ggsave(filename = pdf_file, plot = p, width = 11, height = 6.5, device = cairo_pdf)
write.csv(sig_df[, c("Region", "TestMethod", "PValue", "SignifLabel")], pvalue_file, row.names = FALSE)
write.csv(
  data.frame(
    TestMethod = "Chi-square test of independence (2x7)",
    Statistic = unname(chi_result$statistic),
    DF = unname(chi_result$parameter),
    PValue = chi_result$p.value
  ),
  global_file,
  row.names = FALSE
)

message("Figure saved to:")
message(" - ", png_file)
message(" - ", pdf_file)
message("Significance table saved to:")
message(" - ", pvalue_file)
message("Global test table saved to:")
message(" - ", global_file)
