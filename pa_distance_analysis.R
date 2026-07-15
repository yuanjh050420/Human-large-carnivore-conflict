# 加载必要的库
library(sf)
library(dplyr)
library(ggplot2)
if (!require(patchwork)) {
  install.packages("patchwork")
  library(patchwork)
}  # 用于组合图形

# ==================== 1. Data Reading ====================
# Read China protected area polygon layer
reserves <- st_read("C:\\Users\\yuanj\\Desktop\\ArcGIS\\China Reserve\\China Reserve.shp")

# Read human-wildlife conflict point layer
layers <- st_layers("C:\\Users\\yuanj\\Desktop\\人兽冲突\\大型食肉动物冲突事件\\新闻事件汇总\\ALL.gdb")
print(layers$name)  # View available layer names
if (length(layers$name) == 0) {
  stop("No layers found in ALL.gdb. Please check the path and file.")
}
conflicts <- st_read("C:\\Users\\yuanj\\Desktop\\人兽冲突\\大型食肉动物冲突事件\\新闻事件汇总\\ALL.gdb", 
                     layer = layers$name[1])

# Unified coordinate system (using conflict point CRS)
if (st_crs(reserves) != st_crs(conflicts)) {
  reserves <- st_transform(reserves, crs = st_crs(conflicts))
}

# ==================== 2. Calculate shortest distance from conflict points to protected area ====================
# Vector method (efficient, no loops needed)
nearest_idx <- st_nearest_feature(conflicts, reserves)  # Index of nearest protected area for each conflict point
if (length(nearest_idx) == 0) {
  stop("No nearest protected area found. Please check the data.")
}

# Get boundaries of nearest reserves
nearest_boundaries <- st_boundary(reserves[nearest_idx, , drop = FALSE])

# Calculate distances to boundaries (pairwise, avoids building a full matrix)
distances <- st_distance(conflicts, nearest_boundaries, by_element = TRUE)
distances <- as.numeric(distances)  # Convert to numeric vector (unit: meters)

# Check if points are inside the nearest reserve
is_inside <- st_within(conflicts, reserves[nearest_idx, , drop = FALSE], sparse = FALSE)
is_inside <- apply(is_inside, 1, any)  # Convert to logical vector

# For points inside, make distance negative
distances[is_inside] <- -distances[is_inside]

# Add distance column
conflicts$distance <- distances

# View distance distribution summary
summary(conflicts$distance)

# ==================== 3. Generate background random points (for comparison) ====================
# Read China boundary from GDB as study area
china_layers <- st_layers("C:\\Users\\yuanj\\Desktop\\ArcGIS\\CN\\china.gdb")
china <- st_read("C:\\Users\\yuanj\\Desktop\\ArcGIS\\CN\\china.gdb", layer = china_layers$name[1])
# Unify CRS with conflict points
if (st_crs(china) != st_crs(conflicts)) {
  china <- st_transform(china, crs = st_crs(conflicts))
}
study_area <- st_union(china)  # Merge all features into a single polygon
if (nrow(conflicts) == 0) {
  stop("Conflict point data is empty. Please check the data reading.")
}
set.seed(123)  # Set random seed for reproducibility
bg_multiplier <- 10
background_points <- st_sample(study_area, size = bg_multiplier * nrow(conflicts), type = "random") %>%
  st_as_sf()

# Calculate distance from background points to protected area
bg_nearest <- st_nearest_feature(background_points, reserves)

# Get boundaries of nearest reserves for background points
bg_nearest_boundaries <- st_boundary(reserves[bg_nearest, , drop = FALSE])

# Calculate distances to boundaries (pairwise, avoids building a full matrix)
bg_distances <- st_distance(background_points, bg_nearest_boundaries, by_element = TRUE)
bg_distances <- as.numeric(bg_distances)

# Check if background points are inside the nearest reserve
bg_is_inside <- st_within(background_points, reserves[bg_nearest, , drop = FALSE], sparse = FALSE)
bg_is_inside <- apply(bg_is_inside, 1, any)  # Convert to logical vector

# For points inside, make distance negative
bg_distances[bg_is_inside] <- -bg_distances[bg_is_inside]

background_points$distance <- bg_distances

# Check for NA or invalid distances
if (any(is.na(distances)) || any(is.na(bg_distances))) {
  warning("NA values detected in distance calculation, possibly due to geometry issues.")
}

# ==================== 4. Visualization method 1: Density curve comparison (recommended) ====================
# This method directly shows the probability density distribution, making it easy to compare distribution differences
p1 <- ggplot() +
  geom_density(data = conflicts, aes(x = distance, color = "Conflict events", fill = "Conflict events"), 
               alpha = 0.4, bw = 3000) +  # bw = 3000 meters, bandwidth adjustable
  geom_density(data = background_points, aes(x = distance, color = "Random background points", fill = "Random background points"), 
               alpha = 0.4, bw = 3000) +
  scale_color_manual(values = c("Conflict events" = "red", "Random background points" = "blue")) +
  scale_fill_manual(values = c("Conflict events" = "red", "Random background points" = "blue")) +
  scale_x_continuous(breaks = seq(-50000, 250000, by = 50000), labels = function(x) paste0(x/1000, "km")) +
  coord_cartesian(xlim = c(-50000, 250000)) +
  scale_y_continuous(labels = scales::scientific) +
  labs(
    x = "Distance to protected area (negative for inside)",
    y = "Probability density",
    color = "Type",
    fill = "Type"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = c(0.85, 0.85)
  )

# ==================== 5. Visualization method 2: Grouped frequency + area correction ====================
# Group distances and calculate conflict density for each group (using background points as proxy for accessible area)

# Define distance bins (2 km interval) over all points
bin_width <- 2000

conflict_bins <- data.frame(distance = conflicts$distance) %>%
  mutate(bin_start = floor(distance / bin_width) * bin_width)

bg_bins <- data.frame(distance = background_points$distance) %>%
  mutate(bin_start = floor(distance / bin_width) * bin_width)

# Calculate frequencies by bin and convert to proportions using all points as denominator
min_bg_count <- 3
group_summary <- full_join(
  count(conflict_bins, bin_start, name = "conflict_count"),
  count(bg_bins, bin_start, name = "bg_count"),
  by = "bin_start"
) %>%
  mutate(
    conflict_count = ifelse(is.na(conflict_count), 0, conflict_count),
    bg_count = ifelse(is.na(bg_count), 0, bg_count),
    mid_dist = bin_start + bin_width / 2,
    conflict_prop = conflict_count / nrow(conflicts),
    bg_prop = bg_count / nrow(background_points),
    # Avoid unstable RR in bins with too few background points
    relative_risk = ifelse(bg_count >= min_bg_count, conflict_prop / bg_prop, NA_real_)
  ) %>%
  arrange(mid_dist)

# Use robust upper limit so a few extreme bins do not flatten the whole curve
rr_ymax <- as.numeric(quantile(group_summary$relative_risk, probs = 0.99, na.rm = TRUE))
if (!is.finite(rr_ymax) || rr_ymax <= 0) {
  rr_ymax <- 3
}

# Plot relative risk curve
p2 <- ggplot(group_summary, aes(x = mid_dist, y = relative_risk)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +  # Baseline risk=1
  geom_point(size = 3, color = "darkred", alpha = 0.7, na.rm = TRUE) +
  geom_smooth(method = "loess", se = TRUE, span = 0.2, color = "red", fill = "pink", alpha = 0.3, na.rm = TRUE) +
  labs(
    x = "Distance to protected area (negative for inside)",
    y = "Relative Risk Ratio"
  ) +
  scale_x_continuous(breaks = seq(-50000, 250000, by = 50000), labels = function(x) paste0(x/1000, "km")) +
  coord_cartesian(xlim = c(-50000, 250000), ylim = c(0, rr_ymax)) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  )

# ==================== 6. Combined output ====================
# Combine two graphs
combined_plot <- (p1 / p2) + patchwork::plot_layout(ncol = 1, heights = c(1, 1))

# 显示组合图
if (interactive()) {
  print(combined_plot)
}

# ==================== 7. Statistical testing ====================
# Test whether distance distributions of conflict points and background points are significantly different
ks_test <- ks.test(conflicts$distance, background_points$distance)
print(ks_test)

# Save plots to desktop
output_plot <- "C:\\Users\\yuanj\\Desktop\\conflict_distance_plot_最终.png"
ggsave(
  output_plot,
  plot = combined_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white",
  device = "png",
  limitsize = FALSE
)
message("Saved updated plot to: ", output_plot)
