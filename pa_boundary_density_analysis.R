# =============================================================================
# =============================================================================

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

reserves <- st_read("../ArcGIS/China Reserve/China Reserve.shp")

layers <- st_layers("大型食肉动物冲突事件/新闻事件汇总/ALL.gdb")
conflicts <- st_read("大型食肉动物冲突事件/新闻事件汇总/ALL.gdb",
                     layer = layers$name[1])

if (st_crs(reserves) != st_crs(conflicts)) {
  reserves <- st_transform(reserves, crs = st_crs(conflicts))
}

reserves <- st_make_valid(reserves)

print("正在按几何连通关系合并保护区...")
original_reserve_count <- nrow(reserves)
neighbors <- st_intersects(reserves)

component_id <- integer(length(neighbors))
current_component <- 0

for (i in seq_along(neighbors)) {
  if (component_id[i] != 0) next

  current_component <- current_component + 1
  stack <- i
  component_id[i] <- current_component

  while (length(stack) > 0) {
    node <- stack[length(stack)]
    stack <- stack[-length(stack)]
    next_nodes <- neighbors[[node]]

    unvisited <- next_nodes[component_id[next_nodes] == 0]
    if (length(unvisited) > 0) {
      component_id[unvisited] <- current_component
      stack <- c(stack, unvisited)
    }
  }
}

reserves <- reserves %>%
  mutate(merge_id = component_id) %>%
  group_by(merge_id) %>%
  summarise(do_union = TRUE, .groups = "drop")

print(paste("合并完成：由", original_reserve_count, "个保护区合并为", nrow(reserves), "个连通保护区块"))

reserves <- st_simplify(reserves, dTolerance = 100)

reserves$boundary_length_km <- as.numeric(st_length(st_boundary(reserves))) / 1000

buffer_10km <- st_buffer(reserves, dist = 10000)

buffer_inner <- st_buffer(reserves, dist = -10000)

buffer_outer_geom <- st_geometry(buffer_10km)
buffer_inner_geom <- st_geometry(buffer_inner)
inner_empty <- st_is_empty(buffer_inner_geom)

buffer_ring_geom <- st_sfc(lapply(seq_along(buffer_outer_geom), function(i) {
  if (inner_empty[i]) {
    buffer_outer_geom[[i]]
  } else {
    st_difference(buffer_outer_geom[i], buffer_inner_geom[i])[[1]]
  }
}), crs = st_crs(buffer_10km))

buffer_ring_10km <- st_set_geometry(buffer_10km, buffer_ring_geom)
inner_empty_count <- sum(inner_empty)
print(paste("向内收缩10km后核心区为空的保护区数量：", inner_empty_count))

conflicts_in_ring <- st_intersects(buffer_ring_10km, conflicts, sparse = FALSE)
conflict_count_10km <- rowSums(conflicts_in_ring)

buffer_area_10km <- as.numeric(st_area(buffer_ring_10km)) / 1e6

reserves$conflict_count_10km <- conflict_count_10km
reserves$buffer_area_10km <- buffer_area_10km

reserves$area_km2 <- as.numeric(st_area(reserves)) / 1e6

reserves$conflict_density_10km <- ifelse(
  buffer_area_10km > 0,
  (conflict_count_10km / buffer_area_10km) * 100,
  NA_real_
)

reserves$boundary_density <- reserves$boundary_length_km / reserves$area_km2

reserves_filtered <- reserves %>% filter(area_km2 > 1)
print(paste("筛选后保护区数量：", nrow(reserves_filtered), "/ 总数：", nrow(reserves)))

correlation_data <- reserves_filtered %>%
  st_drop_geometry() %>%
  filter(
    is.finite(boundary_density),
    is.finite(conflict_density_10km)
  )

plot_data <- correlation_data %>%
  filter(boundary_density > 0)

plot_data_positive <- plot_data %>%
  filter(conflict_density_10km > 0)

y_floor <- min(plot_data_positive$conflict_density_10km, na.rm = TRUE) / 5
plot_data_points <- plot_data %>%
  mutate(conflict_density_plot = ifelse(conflict_density_10km > 0, conflict_density_10km, y_floor))

zero_conflict_count <- sum(plot_data$conflict_density_10km == 0)

print(paste("用于相关性分析的保护区数量：", nrow(correlation_data)))
print(paste("用于作图的保护区数量：", nrow(plot_data), "/ 因x轴非正值被排除：", nrow(correlation_data) - nrow(plot_data)))
print(paste("其中冲突密度为0的保护区数量：", zero_conflict_count))
print(paste("用于对数作图与平滑拟合的正值保护区数量：", nrow(plot_data_positive)))
print(paste("用于显示0值点的纵轴底线：", signif(y_floor, 3)))

cor_test <- cor.test(
  correlation_data$boundary_density,
  correlation_data$conflict_density_10km,
  method = "spearman",
  exact = FALSE
)
print("边界密度与冲突密度的Spearman相关性（面积>1km²）：")
print(cor_test)

p <- ggplot(plot_data_points, aes(x = boundary_density, y = conflict_density_plot)) +
  geom_point(size = 3, alpha = 0.7, color = "purple") +
  geom_smooth(
    data = plot_data_positive,
    aes(x = boundary_density, y = conflict_density_10km),
    method = "loess",
    span = 0.9,
    se = TRUE,
    color = "red",
    fill = "pink",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(
    labels = scales::comma,
    breaks = c(0.01, 0.1, 1, 10, 100)
  ) +
  coord_cartesian(ylim = c(y_floor, NA)) +
  labs(
    x = expression("Boundary density (km km"^-2*", log scale)"),
    y = expression("Conflict density (events per 100 km"^2*", log scale)")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  )

output_plot <- "protected_area/boundary_vs_conflict_density.png"
ggsave(output_plot, plot = p, width = 8, height = 6, dpi = 300)

output_info <- file.info(output_plot)
print(paste("图已保存至:", normalizePath(output_plot, winslash = "/", mustWork = FALSE)))
print(paste("文件是否存在:", file.exists(output_plot)))
print(paste("文件大小（字节）:", output_info$size))
