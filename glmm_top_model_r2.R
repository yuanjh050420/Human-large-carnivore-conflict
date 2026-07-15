suppressPackageStartupMessages({
  library(data.table)
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

ensure_package("glmmTMB")
ensure_package("performance")

input_file <- "GLMM/GLMM_main_9vars_standardized_full_coords_sdm_binary.csv"
output_dir <- "GLMM/glmm_outputs_stage3_binary_all_subsets"
best_model_file <- file.path(output_dir, "best_model_ic_stage3_binary_all_subsets.csv")

if (!file.exists(best_model_file)) {
  stop(sprintf("缺少最优模型信息文件: %s", best_model_file))
}

required_cols <- c(
  "Join_Count", "X_aea_m", "Y_aea_m",
  "gdp", "hm", "sdm_binary", "poaching_pressure",
  "pop_density", "dis_to_pa", "dem_250m", "dem_250m_tri"
)

block_size_m <- 100000

make_cond_formula <- function(cond_terms) {
  rhs <- if (length(cond_terms) > 0) paste(cond_terms, collapse = " + ") else "1"
  as.formula(paste("Join_Count ~", rhs, "+ (1 | space_block)"))
}

make_zi_formula <- function(zi_terms) {
  rhs <- if (length(zi_terms) > 0) paste(zi_terms, collapse = " + ") else "1"
  as.formula(paste("~", rhs))
}

parse_terms <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(x) == "" || trimws(x) == "1") {
    return(character(0))
  }
  trimws(strsplit(x, " \\+ ", perl = TRUE)[[1]])
}

safe_write_csv <- function(obj, filename, row.names = FALSE) {
  tryCatch(
    {
      write.csv(obj, file.path(output_dir, filename), row.names = row.names)
      TRUE
    },
    error = function(e) {
      message(sprintf("写入失败 [%s]: %s", filename, e$message))
      FALSE
    }
  )
}

safe_write_lines <- function(text, filename) {
  tryCatch(
    {
      writeLines(text, con = file.path(output_dir, filename), useBytes = TRUE)
      TRUE
    },
    error = function(e) {
      message(sprintf("写入失败 [%s]: %s", filename, e$message))
      FALSE
    }
  )
}

dat <- fread(input_file, encoding = "UTF-8")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop(sprintf("缺少必要字段: %s", paste(missing_cols, collapse = ", ")))
}

model_dat <- as.data.frame(dat[, ..required_cols])
for (nm in names(model_dat)) {
  model_dat[[nm]] <- suppressWarnings(as.numeric(model_dat[[nm]]))
}
model_dat <- model_dat[complete.cases(model_dat), ]
if (nrow(model_dat) < 100) stop("有效样本量过小，无法稳定建模。")
if (any(model_dat$Join_Count < 0, na.rm = TRUE)) stop("Join_Count 存在负值。")

model_dat$sdm <- model_dat$sdm_binary
model_dat$X_block <- floor(model_dat$X_aea_m / block_size_m)
model_dat$Y_block <- floor(model_dat$Y_aea_m / block_size_m)
model_dat$space_block <- factor(paste(model_dat$X_block, model_dat$Y_block, sep = "_"))

best_model_info <- read.csv(best_model_file, stringsAsFactors = FALSE)
if (nrow(best_model_info) < 1) {
  stop("最优模型信息文件为空，无法计算 R²。")
}

best_row <- best_model_info[1, , drop = FALSE]
cond_terms <- parse_terms(best_row$cond_terms)
zi_terms <- parse_terms(best_row$zi_terms)

best_fit <- glmmTMB::glmmTMB(
  formula = make_cond_formula(cond_terms),
  ziformula = make_zi_formula(zi_terms),
  family = glmmTMB::nbinom2,
  data = model_dat
)

r2_obj <- performance::r2_nakagawa(best_fit)

r2_df <- data.frame(
  model = "top_model_best_by_AICc",
  cond_terms = if (length(cond_terms) > 0) paste(cond_terms, collapse = " + ") else "1",
  zi_terms = if (length(zi_terms) > 0) paste(zi_terms, collapse = " + ") else "1",
  R2_marginal = unname(r2_obj$R2_marginal),
  R2_conditional = unname(r2_obj$R2_conditional),
  R2_random_effect_gap = unname(r2_obj$R2_conditional - r2_obj$R2_marginal),
  stringsAsFactors = FALSE
)

safe_write_csv(r2_df, "top_model_r2_nakagawa_stage3_binary.csv", row.names = FALSE)

summary_lines <- c(
  "Stage3 binary top model R2 计算完成",
  sprintf("cond_terms: %s", r2_df$cond_terms[1]),
  sprintf("zi_terms: %s", r2_df$zi_terms[1]),
  sprintf("marginal R2: %.6f", r2_df$R2_marginal[1]),
  sprintf("conditional R2: %.6f", r2_df$R2_conditional[1]),
  sprintf("conditional - marginal: %.6f", r2_df$R2_random_effect_gap[1])
)

safe_write_lines(summary_lines, "top_model_r2_nakagawa_stage3_binary.txt")

print(r2_df)
