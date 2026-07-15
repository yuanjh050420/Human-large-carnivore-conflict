suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(sf)
})

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("检测到缺少包 %s，尝试自动安装...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("包 %s 安装或加载失败，请手动安装后重试。", pkg))
  }
}

ensure_rjava <- function() {
  conda_java_home <- Sys.getenv("CONDA_PREFIX", unset = "")
  if (!nzchar(Sys.getenv("JAVA_HOME")) && nzchar(conda_java_home) && dir.exists(conda_java_home)) {
    Sys.setenv(JAVA_HOME = conda_java_home)
  }

  ensure_package("rJava")

  ok <- tryCatch(
    {
      suppressPackageStartupMessages(library(rJava))
      rJava::.jinit()
      TRUE
    },
    error = function(e) {
      message(sprintf("rJava JVM initialization failed: %s", conditionMessage(e)))
      FALSE
    }
  )

  if (!ok) {
    stop(
      paste(
        "rJava is installed but JVM initialization failed.",
        "Please ensure JDK is installed, or set JAVA_HOME and retry.",
        "For example: Sys.setenv(JAVA_HOME = Sys.getenv('CONDA_PREFIX'))",
        sep = "\n"
      )
    )
  }
}

ensure_package("dismo")
ensure_package("sf")
ensure_package("png")
ensure_rjava()

library(dismo)

u <- function(...) intToUtf8(strtoi(c(...), 16L))

root_dir <- "."
output_dir <- file.path(root_dir, "MaxEnt", "maxext_without_livestock")
aligned_dir <- file.path(output_dir, "01_aligned_predictors")
model_dir <- file.path(output_dir, "02_maxent_model")
prediction_dir <- file.path(output_dir, "03_prediction")
report_dir <- file.path(output_dir, "00_reports")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aligned_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

pop_dem_root <- file.path(root_dir, paste0("pop_dem_gdp", u("6570", "636e")))

paths <- list(
  conflict_gdb = file.path(root_dir, "大型食肉动物冲突事件", "新闻事件汇总", "ALL.gdb"),
  conflict_layer = "All_species_mainland",
  hm = file.path(root_dir, "HM_analysis", "HM2015.gdb"),
  sdm = file.path(root_dir, u("6574", "4f53", "4f5c", "56fe"), u("4e30", "5bcc", "5ea6", "52a0", "51b2", "7a81"), "carnivore_richness.tif"),
  pop_density = file.path(
    pop_dem_root,
    paste0(u("4e2d", "56fd"), "100m", u("4eba", "53e3", "5bc6", "5ea6", "6570", "636e", "96c6"), u("ff08"), "2015", u("5e74", "ff09")),
    "worldpop_2015_merged.tif"
  ),
  dis_to_pa = file.path(root_dir, "GLMM", u("56fe", "5c42"), "reserve_signed_distance_1km.tif"),
  dem = file.path(pop_dem_root, "dem_250m", "dem_250m"),
  tri = file.path(pop_dem_root, "dem_250m", "dem_250m_TRI.tif"),
  gdp = file.path(root_dir, "GDP", "GDP", u("6570", "636e", "0047", "0044", "0050", "0032", "0030", "0031", "0035", "5e74"), "gdp2015"),
  poaching = file.path(root_dir, "poaching", paste0(u("6bcf", "0031", "0030", "4e07", "4eba", "76d7", "730e", "6805", "683c"), ".gdb")),
  china_mask_gdb = "../ArcGIS/CN/Mainland.gdb",
  china_border_gdb = "../ArcGIS/CN/china.gdb"
)

required_paths <- c(
  paths$conflict_gdb,
  paths$hm,
  paths$sdm,
  paths$pop_density,
  paths$dis_to_pa,
  paths$dem,
  paths$tri,
  paths$gdp,
  paths$poaching,
  paths$china_mask_gdb,
  paths$china_border_gdb
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(sprintf("以下输入不存在，请检查路径：\n%s", paste(missing_paths, collapse = "\n")))
}

load_vector <- function(path, layer = NULL) {
  if (is.null(layer)) {
    terra::vect(path)
  } else {
    terra::vect(path, layer = layer)
  }
}

load_raster <- function(path) {
  terra::rast(path)
}

crs_txt <- function(x) {
  trimws(terra::crs(x))
}

pick_mask_layer <- function(gdb_path) {
  lyr <- sf::st_layers(gdb_path)
  if (is.null(lyr) || nrow(lyr) == 0) stop("Mainland.gdb 中没有图层。")
  gt <- toupper(as.character(lyr$geomtype))
  idx <- which(gt %in% c("POLYGON", "MULTIPOLYGON"))
  if (length(idx) > 0) return(lyr$name[idx[1]])
  lyr$name[1]
}

fill_na_by_neighborhood <- function(r, window = 3, max_iter = 3) {
  w <- matrix(1, window, window)
  out <- r
  for (i in seq_len(max_iter)) {
    na_now <- terra::global(is.na(out), "sum", na.rm = TRUE)[1, 1]
    if (is.na(na_now) || na_now == 0) break
    out <- terra::focal(out, w = w, fun = mean, na.rm = TRUE, na.policy = "only")
  }
  out
}

sanitize_non_finite <- function(r) {
  terra::app(r, fun = function(v) {
    v[!is.finite(v)] <- NA_real_
    v
  })
}

apply_mask_and_fill <- function(r, mask_vect, fill_window = 3, fill_iter = 3) {
  x <- terra::crop(r, mask_vect)
  x <- terra::mask(x, mask_vect)
  x <- fill_na_by_neighborhood(x, window = fill_window, max_iter = fill_iter)
  x <- sanitize_non_finite(x)
  x <- terra::mask(x, mask_vect)
  x
}

ensure_writable_output <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    unlink(path, force = TRUE)
    if (file.exists(path)) {
      stop(sprintf("输出文件被占用或无写权限，请关闭占用程序后重试：%s", path))
    }
  }
}

target_res_m <- 1000
species_raw <- load_vector(paths$conflict_gdb, layer = paths$conflict_layer)
target_crs_txt <- crs_txt(species_raw)
mask_layer <- pick_mask_layer(paths$china_mask_gdb)
china_mask_raw <- load_vector(paths$china_mask_gdb, layer = mask_layer)
china_mask <- terra::project(china_mask_raw, target_crs_txt)
border_layer <- pick_mask_layer(paths$china_border_gdb)
china_border_raw <- load_vector(paths$china_border_gdb, layer = border_layer)
china_border <- terra::project(china_border_raw, target_crs_txt)
china_border_line <- terra::as.lines(china_border)

report_source <- function(name, obj, kind) {
  if (kind == "vector") {
    data.frame(
      source = name,
      kind = kind,
      crs = crs_txt(obj),
      res_x = NA_real_,
      res_y = NA_real_,
      extent = paste(as.vector(ext(obj)), collapse = ", "),
      action = if (identical(crs_txt(obj), target_crs_txt)) "keep_crs" else "transform",
      stringsAsFactors = FALSE
    )
  } else {
    rr <- res(obj)
    data.frame(
      source = name,
      kind = kind,
      crs = crs_txt(obj),
      res_x = rr[1],
      res_y = rr[2],
      extent = paste(as.vector(ext(obj)), collapse = ", "),
      action = if (identical(crs_txt(obj), target_crs_txt) && all(abs(rr - target_res_m) < 1e-6)) "keep_grid" else "reproject_or_resample",
      stringsAsFactors = FALSE
    )
  }
}

source_info <- list(
  conflict_points = report_source("conflict_points", species_raw, "vector"),
  hm = report_source("hm", load_raster(paths$hm), "raster"),
  sdm = report_source("sdm", load_raster(paths$sdm), "raster"),
  pop_density = report_source("pop_density", load_raster(paths$pop_density), "raster"),
  dis_to_pa = report_source("dis_to_pa", load_raster(paths$dis_to_pa), "raster"),
  dem = report_source("dem", load_raster(paths$dem), "raster"),
  tri = report_source("tri", load_raster(paths$tri), "raster"),
  gdp = report_source("gdp", load_raster(paths$gdp), "raster"),
  poaching = report_source("poaching", load_raster(paths$poaching), "raster")
)

data_check <- do.call(rbind, source_info)
write.csv(data_check, file.path(report_dir, "00_data_check.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(
  c(
    "Data check summary",
    sprintf("Target CRS from conflict points: %s", target_crs_txt),
    sprintf("Target resolution: %sm", target_res_m),
    "",
    capture.output(print(data_check))
  ),
  con = file.path(report_dir, "00_data_check.txt"),
  useBytes = TRUE
)

template_ref <- load_raster(paths$poaching)
template <- terra::project(template_ref, target_crs_txt, res = target_res_m, method = "bilinear")
template <- apply_mask_and_fill(template, china_mask, fill_window = 3, fill_iter = 3)
template_file <- file.path(aligned_dir, "template_1km.tif")
processing_signature <- "maxent_main_v5_all_species_mainland_mainland_mask_hm_sdm_pop_density_dis_to_pa_dem_tri_gdp_poaching_mask_fill_v2_window3_iter3_finite_only"
signature_file <- file.path(aligned_dir, "aligned_signature.txt")
reuse_aligned_predictors <- TRUE
signature_ok <- file.exists(signature_file) && {
  sig_lines <- readLines(signature_file, warn = FALSE)
  length(sig_lines) >= 1 && identical(trimws(sig_lines[1]), processing_signature)
}
if (!signature_ok) {
  message("检测到对齐规则更新或缺少签名文件，将强制重建模板与变量栅格。")
  reuse_aligned_predictors <- FALSE
}

if (file.exists(template_file) && reuse_aligned_predictors) {
  message("检测到已有模板栅格，直接复用。")
  template <- load_raster(template_file)
} else {
  terra::writeRaster(
    template,
    template_file,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW"), datatype = "FLT4S")
  )
}

predictor_specs <- list(
  hm = paths$hm,
  sdm = paths$sdm,
  pop_density = paths$pop_density,
  dis_to_pa = paths$dis_to_pa,
  dem = paths$dem,
  tri = paths$tri,
  gdp = paths$gdp,
  poaching = paths$poaching
)

aligned_files <- character(0)
aligned_needed <- file.path(aligned_dir, paste0(names(predictor_specs), "_aligned.tif"))
aligned_ready <- reuse_aligned_predictors && file.exists(template_file) && all(file.exists(aligned_needed))
if (aligned_ready) {
  message("检测到已对齐的变量栅格，跳过重投影/重采样。")
  aligned_files <- aligned_needed
} else {
  message("开始对齐变量栅格。")
  for (nm in names(predictor_specs)) {
    message(sprintf("对齐变量: %s", nm))
    r <- load_raster(predictor_specs[[nm]])
    r_aligned <- terra::project(r, template, method = "bilinear")
    r_aligned <- apply_mask_and_fill(r_aligned, china_mask, fill_window = 3, fill_iter = 3)
    out_file <- file.path(aligned_dir, paste0(nm, "_aligned.tif"))
    terra::writeRaster(
      r_aligned,
      out_file,
      overwrite = TRUE,
      wopt = list(gdal = c("COMPRESS=LZW"), datatype = "FLT4S")
    )
    aligned_files <- c(aligned_files, out_file)
  }
  writeLines(processing_signature, con = signature_file, useBytes = TRUE)
}

predictor_stack <- raster::stack(aligned_files)
names(predictor_stack) <- names(predictor_specs)

species_proj <- terra::project(species_raw, target_crs_txt)
species_xy <- terra::crds(species_proj)
species_xy <- species_xy[complete.cases(species_xy), , drop = FALSE]

template_raster <- raster::raster(template_file)
cell_id <- raster::cellFromXY(template_raster, species_xy)
keep_idx <- !is.na(cell_id)
species_xy <- species_xy[keep_idx, , drop = FALSE]
cell_id <- cell_id[keep_idx]
unique_idx <- !duplicated(cell_id)
species_xy <- species_xy[unique_idx, , drop = FALSE]

presence_vals <- raster::extract(predictor_stack, species_xy)
presence_keep <- complete.cases(presence_vals)
species_xy <- species_xy[presence_keep, , drop = FALSE]
presence_vals <- presence_vals[presence_keep, , drop = FALSE]

if (nrow(species_xy) < 20) {
  stop("有效冲突点过少，无法稳定进行 MaxEnt 拟合。")
}

presence_df <- data.frame(x = species_xy[, 1], y = species_xy[, 2])
write.csv(presence_df, file.path(report_dir, "01_presences_used.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(presence_vals, file.path(report_dir, "01_presence_predictor_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")

jar_candidates <- c(
  file.path(system.file(package = "dismo"), "java", "maxent.jar"),
  file.path(output_dir, "maxent.jar"),
  Sys.getenv("MAXENT_JAR")
)
jar_candidates <- jar_candidates[nzchar(jar_candidates)]
maxent_jar <- jar_candidates[file.exists(jar_candidates)][1]
if (is.na(maxent_jar) || !nzchar(maxent_jar)) {
  stop(paste0(
    "未找到 maxent.jar。请将 maxent.jar 放到以下任一路径之一后重试：\n",
    paste(unique(jar_candidates), collapse = "\n")
  ))
}

set.seed(20260525)
maxent_args <- c(
  "linear=true",
  "quadratic=true",
  "product=true",
  "hinge=true",
  "threshold=true",
  "autofeature=false",
  "responsecurves=true",
  "jackknife=true",
  "randomseed=true",
  "outputformat=logistic",
  "pictures=true"
)

model_name <- "All_species_mainland_1km_MaxEnt_main"
maxent_fit <- dismo::maxent(
  x = predictor_stack,
  p = presence_df,
  path = model_dir,
  args = maxent_args
)

saveRDS(maxent_fit, file.path(model_dir, "maxent_model.rds"))
results_df <- as.data.frame(maxent_fit@results)
results_df$metric <- rownames(results_df)
rownames(results_df) <- NULL
write.csv(results_df, file.path(model_dir, "maxent_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")

export_default_jackknife <- function(model_path, out_file, scale = 0.72) {
  default_plot <- file.path(model_path, "plots", "species_jacknife.png")
  if (!file.exists(default_plot)) {
    message("未找到 MaxEnt 默认 Jackknife 图，跳过导出。")
    return(invisible(FALSE))
  }

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  img <- png::readPNG(default_plot)
  h <- dim(img)[1]
  w <- dim(img)[2]
  out_w <- max(900L, as.integer(round(w * scale)))
  out_h <- max(600L, as.integer(round(h * scale)))

  png(out_file, width = out_w, height = out_h, bg = "white")
  op <- par(mar = c(0, 0, 0, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  plot.new()
  rasterImage(as.raster(img), 0, 0, 1, 1)
  invisible(TRUE)
}

jackknife_plot_file <- file.path(report_dir, "03_jackknife_default.png")
export_default_jackknife(model_dir, jackknife_plot_file)

prediction_file <- file.path(prediction_dir, "maxent_logistic_risk.tif")
ensure_writable_output(prediction_file)
logistic_map <- raster::predict(
  predictor_stack,
  maxent_fit,
  args = c("outputformat=logistic"),
  progress = "text",
  filename = prediction_file,
  overwrite = TRUE
)

risk_cols <- colorRampPalette(c("#8ecae6", "#fff3a6", "#d7191c"))(120)
legend_ticks <- seq(0, 1, by = 0.2)
legend_labels <- format(legend_ticks, nsmall = 1, trim = TRUE)

risk_with_points_file <- file.path(prediction_dir, "maxent_logistic_risk.png")
ensure_writable_output(risk_with_points_file)
png(risk_with_points_file, width = 1800, height = 1600, res = 180)
par(mar = c(4, 4, 2, 7))
plot(
  logistic_map,
  main = "MaxEnt logistic conflict risk",
  col = risk_cols,
  zlim = c(0, 1),
  axes = FALSE,
  box = FALSE,
  axis.args = list(at = legend_ticks, labels = legend_labels)
)
lines(china_border_line, col = "#2b2b2b", lwd = 1)
points(species_xy, pch = 20, cex = 0.4, col = rgb(0, 0, 0, 0.35))
dev.off()

risk_no_points_file <- file.path(prediction_dir, "risk_sin_points.png")
ensure_writable_output(risk_no_points_file)
png(risk_no_points_file, width = 1800, height = 1600, res = 180)
par(mar = c(4, 4, 2, 7))
plot(
  logistic_map,
  main = "MaxEnt logistic conflict risk",
  col = risk_cols,
  zlim = c(0, 1),
  axes = FALSE,
  box = FALSE,
  axis.args = list(at = legend_ticks, labels = legend_labels)
)
lines(china_border_line, col = "#2b2b2b", lwd = 1)
dev.off()

summary_lines <- c(
  sprintf("Model: %s", model_name),
  sprintf("Output dir: %s", output_dir),
  sprintf("MaxEnt jar: %s", maxent_jar),
  sprintf("Presence count used: %d", nrow(species_xy)),
  sprintf("Predictor files aligned: %d", length(aligned_files)),
  "",
  "Data check:",
  capture.output(print(data_check)),
  "",
  "MaxEnt results:",
  capture.output(print(head(results_df, 30)))
)
writeLines(summary_lines, con = file.path(report_dir, "02_model_summary.txt"), useBytes = TRUE)

message("MaxEnt workflow completed successfully.")
message(sprintf("Results are in: %s", output_dir))
