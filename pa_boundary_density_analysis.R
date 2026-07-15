# =============================================================================
# 保护区边界密度与冲突密度关系分析
# 目标：建立10km缓冲区，计算每个保护区的边界密度和冲突密度，分析二者的关系
# 边界密度 = 保护区边界长度 / 保护区面积（km/km²），衡量保护区形状的复杂程度
# 冲突密度 = 缓冲区内冲突数 / 缓冲区面积 × 100（每100km²的冲突数）
# =============================================================================

# 加载必要的库
library(sf)       # 空间数据处理
library(dplyr)    # 数据操作
library(tidyr)    # 数据整理
library(ggplot2)  # 可视化
library(scales)   # 坐标轴格式化

# ==================== 1. 数据读取 ====================

# 读取中国保护区面图层
reserves <- st_read("C:\\Users\\yuanj\\Desktop\\ArcGIS\\China Reserve\\China Reserve.shp")

# 读取人兽冲突点图层（从GDB中读取第一个图层）
layers <- st_layers("C:\\Users\\yuanj\\Desktop\\人兽冲突\\大型食肉动物冲突事件\\新闻事件汇总\\ALL.gdb")
conflicts <- st_read("C:\\Users\\yuanj\\Desktop\\人兽冲突\\大型食肉动物冲突事件\\新闻事件汇总\\ALL.gdb",
                     layer = layers$name[1])

# 统一坐标系统（以冲突点CRS为准）
if (st_crs(reserves) != st_crs(conflicts)) {
  reserves <- st_transform(reserves, crs = st_crs(conflicts))
}

# 修复无效几何
reserves <- st_make_valid(reserves)

# 将几何上相连/相交的保护区合并为同一块
print("正在按几何连通关系合并保护区...")
original_reserve_count <- nrow(reserves)
neighbors <- st_intersects(reserves)

# 基于邻接关系计算连通分量（不依赖额外包）
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

# 简化几何（公差100米，加速后续缓冲区运算）
reserves <- st_simplify(reserves, dTolerance = 100)

# ==================== 2. 计算边界长度 ====================

# 提取保护区边界线并计算长度（单位：km）
reserves$boundary_length_km <- as.numeric(st_length(st_boundary(reserves))) / 1000

# ==================== 3. 创建10km缓冲区 ====================

# 向外扩展10km（包含保护区本体+外侧10km）
buffer_10km <- st_buffer(reserves, dist = 10000)

# 向内收缩10km（保护区内距边界>10km的核心区域）
buffer_inner <- st_buffer(reserves, dist = -10000)

# 逐保护区构造边界两侧10km环带，避免整层差集产生笛卡尔结果
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

# ==================== 4. 统计10km环带内冲突事件数 ====================

# 直接统计落入环带内的冲突点
conflicts_in_ring <- st_intersects(buffer_ring_10km, conflicts, sparse = FALSE)
conflict_count_10km <- rowSums(conflicts_in_ring)

# ==================== 5. 计算缓冲区面积 ====================

# 直接使用环带几何的实际面积（单位：km²）
buffer_area_10km <- as.numeric(st_area(buffer_ring_10km)) / 1e6

# ==================== 6. 合并数据并计算指标 ====================

reserves$conflict_count_10km <- conflict_count_10km
reserves$buffer_area_10km <- buffer_area_10km

# 计算保护区面积（km²）
reserves$area_km2 <- as.numeric(st_area(reserves)) / 1e6

# 冲突密度：每100km²缓冲区面积的冲突事件数
reserves$conflict_density_10km <- ifelse(
  buffer_area_10km > 0,
  (conflict_count_10km / buffer_area_10km) * 100,
  NA_real_
)

# 边界密度：保护区边界长度 / 保护区面积（km/km²），反映边界形状复杂度
reserves$boundary_density <- reserves$boundary_length_km / reserves$area_km2

# 筛选面积大于1km²的保护区
reserves_filtered <- reserves %>% filter(area_km2 > 1)
print(paste("筛选后保护区数量：", nrow(reserves_filtered), "/ 总数：", nrow(reserves)))

# 相关性分析使用所有有限值；作图保留0值，仅排除x轴非正值
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

# 为严格对数纵轴准备绘图值：0值仅用于散点显示时压到最底线
y_floor <- min(plot_data_positive$conflict_density_10km, na.rm = TRUE) / 5
plot_data_points <- plot_data %>%
  mutate(conflict_density_plot = ifelse(conflict_density_10km > 0, conflict_density_10km, y_floor))

zero_conflict_count <- sum(plot_data$conflict_density_10km == 0)

print(paste("用于相关性分析的保护区数量：", nrow(correlation_data)))
print(paste("用于作图的保护区数量：", nrow(plot_data), "/ 因x轴非正值被排除：", nrow(correlation_data) - nrow(plot_data)))
print(paste("其中冲突密度为0的保护区数量：", zero_conflict_count))
print(paste("用于对数作图与平滑拟合的正值保护区数量：", nrow(plot_data_positive)))
print(paste("用于显示0值点的纵轴底线：", signif(y_floor, 3)))

# ==================== 7. 相关性分析 ====================

# Spearman秩相关（适用于非正态分布数据），仅对筛选后数据
cor_test <- cor.test(
  correlation_data$boundary_density,
  correlation_data$conflict_density_10km,
  method = "spearman",
  exact = FALSE
)
print("边界密度与冲突密度的Spearman相关性（面积>1km²）：")
print(cor_test)

# ==================== 8. 可视化 ====================

# 绘图：横纵坐标均取对数；0值点贴在纵轴底线显示，平滑拟合仅使用冲突密度大于0的样本
p <- ggplot(plot_data_points, aes(x = boundary_density, y = conflict_density_plot)) +
  geom_point(size = 3, alpha = 0.7, color = "purple") +                                      # 散点
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
  ) +       # LOESS平滑趋势线
  scale_x_log10(labels = scales::comma) +                                                     # 横坐标取对数
  scale_y_log10(
    labels = scales::comma,
    breaks = c(0.01, 0.1, 1, 10, 100)
  ) +                                                     # 纵坐标严格对数
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

# 保存图片（300dpi高清输出）
output_plot <- "C:\\Users\\yuanj\\Desktop\\人兽冲突\\protected_area\\boundary_vs_conflict_density.png"
ggsave(output_plot, plot = p, width = 8, height = 6, dpi = 300)

output_info <- file.info(output_plot)
print(paste("图已保存至:", normalizePath(output_plot, winslash = "/", mustWork = FALSE)))
print(paste("文件是否存在:", file.exists(output_plot)))
print(paste("文件大小（字节）:", output_info$size))
