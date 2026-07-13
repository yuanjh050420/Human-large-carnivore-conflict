# SKS-SKP scatter plots for predator and victim groups
# Input: victim.csv
# Output:
#   sks_skp_predator_k2.pdf
#   sks_skp_victim_k2.pdf
#   sks_skp_metrics_predator_k2.csv
#   sks_skp_metrics_victim_k2.csv

input_csv <- "C:/Users/yuanj/Desktop/人兽冲突/victim_analysis/victim.csv"
out_dir <- "C:/Users/yuanj/Desktop/人兽冲突/victim_analysis"
k <- 2

library(ggplot2)

raw <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("No.", "Victem", "Number of deaths")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

raw[["No."]] <- trimws(raw[["No."]])
raw <- raw[!is.na(raw[["No."]]) & raw[["No."]] != "", , drop = FALSE]

raw[["Number of deaths"]] <- suppressWarnings(as.numeric(raw[["Number of deaths"]]))
raw[["Number of deaths"]][is.na(raw[["Number of deaths"]])] <- 0

# Predator code from event ID prefix, e.g. ABB001 -> ABB
raw$PredatorCode <- toupper(sub("[0-9].*$", "", raw[["No."]]))
raw$PredatorCode <- trimws(raw$PredatorCode)
raw$PredatorCode[is.na(raw$PredatorCode) | raw$PredatorCode == ""] <- "UNK"

predator_name_map <- c(
  ABB = "Asiatic black bear",
  AGC = "Asiatic golden cat",
  BB  = "Brown bear",
  CL  = "Clouded leopard",
  D   = "Dhole",
  EL  = "Eurasian lynx",
  GW  = "Grey wolf",
  L   = "Leopard",
  SL  = "Snow leopard",
  T   = "Tiger"
)
raw$Predator <- unname(predator_name_map[raw$PredatorCode])
raw$Predator[is.na(raw$Predator)] <- raw$PredatorCode[is.na(raw$Predator)]

raw[["Victem"]] <- trimws(gsub("\\s+", " ", raw[["Victem"]]))
raw[["Victem"]] <- tools::toTitleCase(tolower(raw[["Victem"]]))
raw[["Victem"]][raw[["Victem"]] == "Cow"] <- "Cattle"
raw <- raw[!is.na(raw[["Victem"]]) & raw[["Victem"]] != "", , drop = FALSE]

calc_metrics <- function(deaths_vec, threshold = 2) {
  n <- length(deaths_vec)
  over_mask <- deaths_vec > threshold

  skp <- if (n > 0) mean(over_mask) else NA_real_                 # SKP = P(D > k)
  sks <- if (any(over_mask)) mean(deaths_vec[over_mask] - threshold) else 0  # SKS = E(D-k | D > k)
  skb <- if (n > 0) mean(pmax(deaths_vec - threshold, 0)) else NA_real_      # SKB = E[(D-k)+]

  c(
    EventCount = n,
    SKEvents = sum(over_mask),
    SKP_k = skp,
    SKS_k = sks,
    SKB_k = skb
  )
}

summarise_group <- function(df, group_col, min_events, threshold = 2) {
  event_level <- aggregate(
    df[["Number of deaths"]],
    by = list(EventID = df[["No."]], Group = df[[group_col]]),
    FUN = sum,
    na.rm = TRUE
  )
  names(event_level)[names(event_level) == "x"] <- "DeathsPerEvent"

  event_counts <- aggregate(event_level$EventID,
                            by = list(Group = event_level$Group),
                            FUN = length)
  names(event_counts)[2] <- "EventCount"

  selected <- event_counts$Group[event_counts$EventCount > min_events]
  if (length(selected) == 0) {
    stop(paste("No group has event count >", min_events))
  }

  sub_events <- event_level[event_level$Group %in% selected, , drop = FALSE]
  group_list <- split(sub_events$DeathsPerEvent, sub_events$Group)

  metric_mat <- t(sapply(group_list, calc_metrics, threshold = threshold))
  metrics <- data.frame(Group = rownames(metric_mat), metric_mat,
                        row.names = NULL, check.names = FALSE)

  for (col_name in setdiff(names(metrics), "Group")) {
    metrics[[col_name]] <- as.numeric(metrics[[col_name]])
  }

  metrics <- metrics[order(-metrics$SKB_k, -metrics$SKP_k, -metrics$SKS_k), ]
  metrics
}

plot_sks_skp <- function(metrics, x_label, out_pdf) {
  group_order <- as.character(metrics$Group)
  metrics$Group <- factor(group_order, levels = group_order)

  predator_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                        "#E69F00", "#56B4E9", "#F0E442", "#000000")
  victim_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                      "#66A61E", "#E6AB02", "#A6761D", "#666666")
  palette_base <- if (identical(x_label, "Victim")) victim_palette else predator_palette
  palette_vals <- setNames(palette_base[seq_len(min(length(palette_base), nrow(metrics)))], group_order)

  label_nudge_y <- max(diff(range(metrics$SKS_k, na.rm = TRUE)) * 0.04, 0.06)
  label_df <- metrics
  label_df$base_hjust <- ifelse(label_df$SKP_k > 0.85, 1, 0)

  # Species-specific label distance tuning for manuscript readability.
  offset_mult <- rep(1, nrow(label_df))
  group_name <- as.character(label_df$Group)
  offset_mult[group_name %in% c("Human", "Cattle", "Brown bear", "Brown Bear")] <- 0.60
  offset_mult[group_name %in% c("Grey wolf", "Gray wolf")] <- 1.75
  is_tiger <- group_name == "Tiger"
  is_leopard <- group_name == "Leopard"

  x_step <- 0.015 * offset_mult
  y_step <- label_nudge_y * offset_mult
  label_df$label_x <- ifelse(label_df$base_hjust == 1, label_df$SKP_k - x_step, label_df$SKP_k + x_step)
  label_df$label_x[is_tiger] <- label_df$SKP_k[is_tiger] - 0.02
  label_df$label_x[is_leopard] <- label_df$SKP_k[is_leopard] + 0.03
  label_df$label_x <- pmin(pmax(label_df$label_x, 0.02), 0.98)
  label_df$label_y <- label_df$SKS_k + y_step +
    ((rank(label_df$SKP_k, ties.method = "first") %% 2) * y_step * 0.35)
  label_df$label_y[is_leopard] <- label_df$SKS_k[is_leopard] - y_step[is_leopard] * 0.35
  label_df$hjust <- ifelse(label_df$base_hjust == 1, 1, 0)
  label_df$hjust[is_tiger] <- 1
  label_df$hjust[is_leopard] <- 0

  skb_upper <- max(20, max(metrics$SKB_k, na.rm = TRUE))

  p <- ggplot(metrics, aes(x = SKP_k, y = SKS_k, size = SKB_k, color = Group)) +
    geom_point(alpha = 0.85) +
    geom_text(
      data = label_df,
      aes(x = label_x, y = label_y, label = Group, hjust = hjust),
      size = 3.5,
      show.legend = FALSE,
      check_overlap = FALSE
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0.01, 0.03))) +
    scale_color_manual(values = palette_vals) +
    scale_size_continuous(
      range = c(3.2, 10),
      breaks = c(5, 10, 15, 20),
      limits = c(0, skb_upper)
    ) +
    guides(
      color = guide_legend(order = 1),
      size = guide_legend(order = 2, override.aes = list(color = "grey45", alpha = 1))
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      caption = NULL,
      x = "Surplus killing probability (SKP)",
      y = "Surplus killing severity (SKS)",
      size = "SKB",
      color = x_label
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 12, r = 26, b = 12, l = 12)
    )

  ggsave(out_pdf, p, device = cairo_pdf, width = 7.2, height = 5.6, units = "in")
}

pred_metrics <- summarise_group(raw, group_col = "Predator", min_events = 20, threshold = k)
vic_metrics  <- summarise_group(raw, group_col = "Victem",  min_events = 10, threshold = k)
vic_metrics$Group <- factor(vic_metrics$Group, levels = as.character(vic_metrics$Group))

if (nrow(pred_metrics) != 6) {
  warning(paste("Predator groups selected:", nrow(pred_metrics), "(expected 6)"))
}
if (nrow(vic_metrics) != 6) {
  warning(paste("Victim groups selected:", nrow(vic_metrics), "(expected 6)"))
}

out_pred_csv <- file.path(out_dir, paste0("sks_skp_metrics_predator_k", k, ".csv"))
out_vic_csv  <- file.path(out_dir, paste0("sks_skp_metrics_victim_k", k, ".csv"))
out_pred_pdf <- file.path(out_dir, paste0("sks_skp_predator_k", k, ".pdf"))
out_vic_pdf  <- file.path(out_dir, paste0("sks_skp_victim_k", k, ".pdf"))

write.csv(pred_metrics, out_pred_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(vic_metrics, out_vic_csv, row.names = FALSE, fileEncoding = "UTF-8")

plot_sks_skp(pred_metrics, x_label = "Predator", out_pdf = out_pred_pdf)
plot_sks_skp(vic_metrics,  x_label = "Victim",   out_pdf = out_vic_pdf)

cat("Done.\n")
cat("Threshold k =", k, "\n")
cat("Saved:\n")
cat(" -", out_pred_pdf, "\n")
cat(" -", out_vic_pdf, "\n")
cat(" -", out_pred_csv, "\n")
cat(" -", out_vic_csv, "\n")
