suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("包 %s 安装或加载失败，请手动安装后重试。", pkg))
  }
}

ensure_package("terra")
ensure_package("sf")

risk_path <- "C:/Users/yuanj/Desktop/人兽冲突/MaxEnt/maxent_main/03_prediction/maxent_logistic_risk.tif"
sdm_path <- "C:/Users/yuanj/Desktop/人兽冲突/MaxEnt/maxent_main/01_aligned_predictors/sdm_aligned.tif"
border_gdb <- "C:/Users/yuanj/Desktop/ArcGIS/CN/china.gdb"
out_dir <- "C:/Users/yuanj/Desktop/人兽冲突/整体作图/丰富度加MaxEnt风险"
target_res_m <- 1000
h_km <- 20
h_m <- h_km * 1000

# 若仅修改制图样式/阈值逻辑，可复用已生成核平滑栅格以避免整段重跑
force_recompute_smoothing <- FALSE

out_risk_s <- file.path(out_dir, "risk_smoothed_1km_h20km.tif")
out_sdm_s <- file.path(out_dir, "sdm_smoothed_1km_h20km.tif")
out_class <- file.path(out_dir, "risk_richness_4class_1km_h20km.tif")

for (p in c(risk_path, sdm_path, border_gdb)) {
  if (!file.exists(p)) stop(sprintf("输入不存在: %s", p))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tmp_dir <- file.path(out_dir, "_terra_tmp")
if (dir.exists(tmp_dir)) {
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = tmp_dir, memfrac = 0.7, progress = 1)

wopt_flt <- list(gdal = c("COMPRESS=LZW"), datatype = "FLT4S")
wopt_int <- list(gdal = c("COMPRESS=LZW"), datatype = "INT1U")

message("[1/7] 读取输入栅格")

load_border <- function(gdb_path, target_crs) {
  lyr <- sf::st_layers(gdb_path)
  if (is.null(lyr) || nrow(lyr) == 0) stop("边境库中未找到图层。")
  gt <- toupper(as.character(lyr$geomtype))
  idx <- which(gt %in% c("POLYGON", "MULTIPOLYGON"))
  layer_name <- if (length(idx) > 0) lyr$name[idx[1]] else lyr$name[1]
  b <- terra::vect(gdb_path, layer = layer_name)
  b <- terra::project(b, target_crs)
  terra::as.lines(b)
}

q75 <- function(r) {
  as.numeric(terra::global(r, fun = quantile, probs = 0.75, na.rm = TRUE)[1, 1])
}

risk <- terra::rast(risk_path)
sdm <- terra::rast(sdm_path)

reuse_smoothed <-
  !force_recompute_smoothing &&
  file.exists(out_risk_s) &&
  file.exists(out_sdm_s)

if (reuse_smoothed) {
  in_times <- file.info(c(risk_path, sdm_path))$mtime
  out_times <- file.info(c(out_risk_s, out_sdm_s))$mtime
  reuse_smoothed <- all(out_times >= in_times, na.rm = TRUE)
}

if (reuse_smoothed) {
  message("[2/7] 检测到已生成核平滑文件，直接复用")
  risk_s <- terra::rast(out_risk_s)
  sdm_s <- terra::rast(out_sdm_s)
  rm(risk, sdm)
  gc()
} else {
  if (force_recompute_smoothing) {
    message("[2/7] 已开启强制重算，忽略已有核平滑文件")
  } else {
    message("[2/7] 未找到可复用核平滑文件（或输入更新），执行重算")
  }

  # 以 MaxEnt 风险图为模板对齐 SDM
  if (!identical(trimws(terra::crs(risk)), trimws(terra::crs(sdm)))) {
    message("[3/7] SDM CRS 不同，执行投影对齐（写临时文件）")
    sdm <- terra::project(
      sdm,
      risk,
      method = "bilinear",
      filename = file.path(tmp_dir, "sdm_projected.tif"),
      overwrite = TRUE,
      wopt = wopt_flt
    )
  } else {
    same_grid <-
      nrow(risk) == nrow(sdm) &&
      ncol(risk) == ncol(sdm) &&
      all(abs(terra::res(risk) - terra::res(sdm)) < 1e-9) &&
      all(abs(as.vector(terra::ext(risk)) - as.vector(terra::ext(sdm))) < 1e-6)
    if (!same_grid) {
      message("[3/7] SDM 网格不同，执行重采样对齐（写临时文件）")
      sdm <- terra::resample(
        sdm,
        risk,
        method = "bilinear",
        filename = file.path(tmp_dir, "sdm_resampled.tif"),
        overwrite = TRUE,
        wopt = wopt_flt
      )
    }
  }

  # 仅保留两层均有值的像元
  message("[4/7] 构建共同有效区并掩膜（写临时文件）")
  common_mask <- terra::ifel(
    is.na(risk) | is.na(sdm),
    NA,
    1,
    filename = file.path(tmp_dir, "common_mask.tif"),
    overwrite = TRUE,
    wopt = wopt_int
  )
  risk <- terra::mask(
    risk,
    common_mask,
    filename = file.path(tmp_dir, "risk_masked.tif"),
    overwrite = TRUE,
    wopt = wopt_flt
  )
  sdm <- terra::mask(
    sdm,
    common_mask,
    filename = file.path(tmp_dir, "sdm_masked.tif"),
    overwrite = TRUE,
    wopt = wopt_flt
  )
  gc()

  # 高斯核平滑，d 使用地图单位（米）
  message(sprintf("[5/7] 高斯核平滑（h=%d km）", h_km))
  w <- terra::focalMat(risk, d = h_m, type = "Gauss")
  risk_s <- terra::focal(
    risk,
    w = w,
    fun = mean,
    na.rm = TRUE
  )
  sdm_s <- terra::focal(
    sdm,
    w = w,
    fun = mean,
    na.rm = TRUE
  )
  risk_s <- terra::mask(
    risk_s,
    common_mask,
    filename = file.path(tmp_dir, "risk_smoothed_masked.tif"),
    overwrite = TRUE,
    wopt = wopt_flt
  )
  sdm_s <- terra::mask(
    sdm_s,
    common_mask,
    filename = file.path(tmp_dir, "sdm_smoothed_masked.tif"),
    overwrite = TRUE,
    wopt = wopt_flt
  )
  rm(risk, sdm)
  gc()
}

message("[5/7] 计算Q75阈值并生成四分类")
risk_thr <- q75(risk_s)
sdm_thr <- q75(sdm_s)

risk_high <- risk_s >= risk_thr
sdm_high <- sdm_s >= sdm_thr

# 1=低风险低丰富度, 2=低风险高丰富度, 3=高风险低丰富度, 4=高风险高丰富度
class4 <- terra::ifel(
  risk_high & sdm_high, 4,
  terra::ifel(
    risk_high & !sdm_high, 3,
    terra::ifel(!risk_high & sdm_high, 2, 1)
  ),
  filename = file.path(tmp_dir, "class4_tmp.tif"),
  overwrite = TRUE,
  wopt = wopt_int
)
names(class4) <- "risk_richness_class"
rm(risk_high, sdm_high)
gc()

message("[6/7] 写出结果栅格")
if (!reuse_smoothed || force_recompute_smoothing) {
  terra::writeRaster(risk_s, out_risk_s, overwrite = TRUE, wopt = wopt_flt)
  terra::writeRaster(sdm_s, out_sdm_s, overwrite = TRUE, wopt = wopt_flt)
}
terra::writeRaster(class4, out_class, overwrite = TRUE, wopt = wopt_int)

labels_df <- data.frame(
  code = 1:4,
  label = c("Low risk - Low richness", "Low risk - High richness", "High risk - Low richness", "High risk - High richness"),
  color = c("#cfe8ff", "#69b3a2", "#f9a66c", "#d7301f"),
  stringsAsFactors = FALSE
)

border_line <- load_border(border_gdb, terra::crs(class4))

message("[7/7] 输出地图和统计")
map_png <- file.path(out_dir, "risk_richness_4class_1km_h20km.png")
png(map_png, width = 5200, height = 3600, res = 400)
layout(matrix(c(1, 2), nrow = 1), widths = c(8.0, 2.8))

par(mar = c(1.8, 3.5, 6.2, 1.0))
# 使用全像元绘图，避免terra为显示速度而抽样导致PNG看起来发糊
plot(class4, col = labels_df$color, type = "classes", axes = FALSE, box = FALSE, legend = FALSE,
  maxcell = ncell(class4), main = "")
lines(border_line, col = "#2b2b2b", lwd = 0.9)

par(mar = c(1.5, 0.5, 1.5, 0.5))
plot.new()
legend(
  "left",
  title = "Class",
  legend = labels_df$label,
  fill = labels_df$color,
  border = "#666666",
  bty = "n",
  cex = 1.1,
  ncol = 1,
  x.intersp = 0.8,
  y.intersp = 1.2
)
layout(1)
dev.off()

freq_df <- terra::freq(class4, digits = 0, value = TRUE)
pixel_area_km2 <- abs(prod(terra::res(class4))) / 1e6
area_df <- data.frame(
  code = as.integer(freq_df$value),
  area_km2 = as.numeric(freq_df$count) * pixel_area_km2,
  stringsAsFactors = FALSE
)
area_df <- merge(labels_df, area_df, by = "code", all.x = TRUE)
area_df$area_km2[is.na(area_df$area_km2)] <- 0

write.csv(area_df, file.path(out_dir, "risk_richness_4class_1km_h20km_area_km2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(
  data.frame(
    metric = c("output_resolution_m", "kernel_bandwidth_km", "risk_q75_threshold", "sdm_q75_threshold"),
    value = c(target_res_m, h_km, risk_thr, sdm_thr)
  ),
  file.path(out_dir, "risk_richness_thresholds_1km_h20km.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary_file <- file.path(out_dir, sprintf("run_summary_%dkm_h%dkm.txt", as.integer(target_res_m / 1000), h_km))

writeLines(
  c(
    "Risk-Richness overlay completed.",
    sprintf("Risk input: %s", risk_path),
    sprintf("SDM input: %s", sdm_path),
    sprintf("Output resolution: %.0f m", target_res_m),
    sprintf("Kernel bandwidth: %d km", h_km),
    sprintf("Risk Q75 threshold: %.6f", risk_thr),
    sprintf("SDM Q75 threshold: %.6f", sdm_thr),
    sprintf("Class tif: %s", out_class),
    sprintf("Map png: %s", map_png)
  ),
  con = summary_file,
  useBytes = TRUE
)

message("Done: Risk-Richness overlay outputs are saved.")
