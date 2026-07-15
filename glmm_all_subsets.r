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
ensure_package("parallel")

input_file <- "GLMM/GLMM_main_9vars_standardized_full_coords_sdm_binary.csv"
output_dir <- "GLMM/glmm_outputs_stage3_binary_all_subsets"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

block_size_m <- 100000
parallel_workers <- max(1, parallel::detectCores(logical = TRUE) - 1)
aicc_delta_cutoff <- 2
cum_weight_cutoff <- 0.95

use_cache <- TRUE
cache_file <- file.path(output_dir, "all_subsets_fit_cache_stage3_binary.rds")
cache_meta_file <- file.path(output_dir, "all_subsets_fit_cache_stage3_binary_meta.rds")
force_rebuild_cache <- TRUE
posthoc_from_cache_only <- FALSE

dat <- fread(input_file, encoding = "UTF-8")
required_cols <- c(
  "Join_Count", "X_aea_m", "Y_aea_m",
  "gdp", "hm", "sdm_binary", "poaching_pressure",
  "pop_density", "dis_to_pa", "dem_250m", "dem_250m_tri"
)
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

zero_prop <- mean(model_dat$Join_Count == 0)
cat(sprintf("样本量: %d\n", nrow(model_dat)))
cat(sprintf("零值比例: %.4f\n", zero_prop))
cat(sprintf("空间 block 数量: %d\n", nlevels(model_dat$space_block)))

cond_main_terms <- c(
  "gdp", "hm", "sdm", "poaching_pressure",
  "pop_density", "dis_to_pa", "dem_250m", "dem_250m_tri"
)
cond_inter_terms <- c("sdm:hm")
cond_quad_terms <- c("I(hm^2)")
cond_all_terms <- c(cond_main_terms, cond_inter_terms, cond_quad_terms)

zi_candidate_terms <- c("pop_density", "sdm")

ensure_output_dir <- function() {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(output_dir)
}

out_file <- function(name) {
  ensure_output_dir()
  file.path(output_dir, name)
}

safe_write_csv <- function(obj, filename, row.names = FALSE) {
  tryCatch(
    {
      write.csv(obj, out_file(filename), row.names = row.names)
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
      writeLines(text, con = out_file(filename), useBytes = TRUE)
      TRUE
    },
    error = function(e) {
      message(sprintf("写入失败 [%s]: %s", filename, e$message))
      FALSE
    }
  )
}

safe_capture_output <- function(expr, filename) {
  tryCatch(
    {
      capture.output(expr, file = out_file(filename))
      TRUE
    },
    error = function(e) {
      message(sprintf("写入失败 [%s]: %s", filename, e$message))
      FALSE
    }
  )
}

make_cond_formula <- function(cond_terms) {
  rhs <- if (length(cond_terms) > 0) paste(cond_terms, collapse = " + ") else "1"
  as.formula(paste("Join_Count ~", rhs, "+ (1 | space_block)"))
}

make_zi_formula <- function(zi_terms) {
  rhs <- if (length(zi_terms) > 0) paste(zi_terms, collapse = " + ") else "1"
  as.formula(paste("~", rhs))
}

calc_aicc <- function(fit) {
  aic <- as.numeric(AIC(fit))
  ll <- logLik(fit)
  k <- attr(ll, "df")
  n <- nobs(fit)
  denom <- n - k - 1
  if (is.na(denom) || denom <= 0) return(Inf)
  aic + (2 * k * (k + 1)) / denom
}

is_fit_ok <- function(fit_obj, warns) {
  if (is.null(fit_obj)) return(FALSE)

  pd_hess_ok <- !is.null(fit_obj$sdr) && isTRUE(fit_obj$sdr$pdHess)
  conv_code_ok <- !is.null(fit_obj$fit$convergence) && fit_obj$fit$convergence == 0
  singular_warn <- any(grepl(
    "singular convergence|non-positive-definite|false convergence",
    warns,
    ignore.case = TRUE
  ))

  pd_hess_ok && conv_code_ok && !singular_warn
}

make_model_key <- function(cond_terms, zi_terms) {
  cond_key <- if (length(cond_terms) == 0) "1" else paste(sort(cond_terms), collapse = " + ")
  zi_key <- if (length(zi_terms) == 0) "1" else paste(sort(zi_terms), collapse = " + ")
  paste0("cond:", cond_key, " || zi:", zi_key)
}

make_cache_signature <- function() {
  fi <- file.info(input_file)
  if (nrow(fi) == 0 || is.na(fi$size[1])) {
    stop(sprintf("输入文件不存在或不可读: %s", input_file))
  }

  paste(
    normalizePath(input_file, winslash = "/", mustWork = FALSE),
    as.character(fi$size[1]),
    as.character(fi$mtime[1]),
    as.character(block_size_m),
    paste(required_cols, collapse = ","),
    paste(cond_main_terms, collapse = ","),
    paste(cond_inter_terms, collapse = ","),
    paste(cond_quad_terms, collapse = ","),
    paste(zi_candidate_terms, collapse = ","),
    sep = "|"
  )
}

is_cond_subset_valid <- function(term_subset) {
  if ("sdm:hm" %in% term_subset && !all(c("sdm", "hm") %in% term_subset)) return(FALSE)
  if ("I(hm^2)" %in% term_subset && !("hm" %in% term_subset)) return(FALSE)
  TRUE
}

all_subsets <- function(terms) {
  n <- length(terms)
  out <- vector("list", 2^n)
  idx <- 1L
  for (k in 0:n) {
    cmb <- combn(terms, k, simplify = FALSE)
    for (s in cmb) {
      out[[idx]] <- s
      idx <- idx + 1L
    }
  }
  out
}

extract_coef_table <- function(fit, part = c("cond", "zi")) {
  part <- match.arg(part)
  sm <- summary(fit)$coefficients
  if (!part %in% names(sm)) return(data.frame())
  tab <- as.data.frame(sm[[part]])
  if (nrow(tab) == 0) return(data.frame())
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab$part <- part
  tab
}

# full averaging + unconditional SE
model_average_part <- function(selected_results, part = c("cond", "zi")) {
  part <- match.arg(part)

  term_union <- character(0)
  for (res in selected_results) {
    tb <- res$coef_tables[[part]]
    if (!is.null(tb) && nrow(tb) > 0) {
      term_union <- union(term_union, tb$term)
    }
  }
  term_union <- sort(unique(term_union))

  if (length(term_union) == 0) return(data.frame())

  w <- sapply(selected_results, function(x) x$weight)
  avg_beta <- numeric(length(term_union))
  uncond_se <- numeric(length(term_union))

  names(avg_beta) <- term_union
  names(uncond_se) <- term_union

  for (tm in term_union) {
    beta_i <- numeric(length(selected_results))
    se_i <- numeric(length(selected_results))

    for (i in seq_along(selected_results)) {
      tb <- selected_results[[i]]$coef_tables[[part]]
      if (!is.null(tb) && nrow(tb) > 0 && any(tb$term == tm)) {
        rw <- tb[tb$term == tm, , drop = FALSE][1, ]
        beta_i[i] <- as.numeric(rw$Estimate)
        se_i[i] <- as.numeric(rw$`Std. Error`)
      } else {
        beta_i[i] <- 0
        se_i[i] <- 0
      }
    }

    bbar <- sum(w * beta_i)
    var_u <- sum(w * (se_i^2 + (beta_i - bbar)^2))

    avg_beta[tm] <- bbar
    uncond_se[tm] <- sqrt(var_u)
  }

  data.frame(
    term = term_union,
    Estimate = as.numeric(avg_beta),
    Uncond_SE = as.numeric(uncond_se),
    part = part,
    stringsAsFactors = FALSE
  )
}

cond_subsets_raw <- all_subsets(cond_all_terms)
cond_subsets <- Filter(is_cond_subset_valid, cond_subsets_raw)
zi_subsets <- all_subsets(zi_candidate_terms)

cat(sprintf("条件部分候选子集: %d（层级约束后）\n", length(cond_subsets)))
cat(sprintf("零膨胀部分候选子集: %d\n", length(zi_subsets)))
cat(sprintf("理论总组合数: %d\n", length(cond_subsets) * length(zi_subsets)))

model_tasks <- list()
for (cset in cond_subsets) {
  for (zset in zi_subsets) {
    key <- make_model_key(cset, zset)
    model_tasks[[length(model_tasks) + 1L]] <- list(
      key = key,
      cond_terms = cset,
      zi_terms = zset
    )
  }
}

task_keys <- vapply(model_tasks, function(x) x$key, character(1))
model_tasks <- model_tasks[!duplicated(task_keys)]
cat(sprintf("去重后待拟合模型数: %d\n", length(model_tasks)))

cached_results <- list()
cache_signature <- make_cache_signature()

if (force_rebuild_cache) {
  cat("force_rebuild_cache=TRUE：本次忽略历史缓存并全量重跑。\n")
} else if (use_cache && file.exists(cache_file)) {
  cache_meta <- if (file.exists(cache_meta_file)) {
    tryCatch(readRDS(cache_meta_file), error = function(e) NULL)
  } else {
    NULL
  }

  meta_ok <- !is.null(cache_meta) && is.list(cache_meta) &&
    !is.null(cache_meta$signature) && identical(cache_meta$signature, cache_signature)

  if (!meta_ok) {
    cat("检测到缓存签名不匹配或缺失元数据：将忽略旧缓存并重跑。\n")
  } else {
    tmp <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.list(tmp) && length(tmp) > 0) {
      cached_results <- tmp
      cat(sprintf("已加载缓存结果: %d\n", length(cached_results)))
    }
  }
}

cached_keys <- names(cached_results)
if (is.null(cached_keys)) cached_keys <- character(0)

pending_tasks <- Filter(function(x) !(x$key %in% cached_keys), model_tasks)
cat(sprintf("本次需新拟合模型数: %d\n", length(pending_tasks)))

if (posthoc_from_cache_only && length(pending_tasks) > 0) {
  stop(sprintf(
    paste0(
      "当前设置为 posthoc_from_cache_only=TRUE（仅后处理，不重跑模型），",
      "但缓存缺少 %d 个模型结果。请先关闭该开关补齐缓存，或继续使用现有 ΔAICc<=2 结果。"
    ),
    length(pending_tasks)
  ))
}

fit_one_task <- function(task, dat_local) {
  key <- task$key
  cond_terms <- task$cond_terms
  zi_terms <- task$zi_terms

  cond_formula <- make_cond_formula(cond_terms)
  zi_formula <- make_zi_formula(zi_terms)

  fit_once <- function(ctrl = NULL) {
    warns <- character(0)
    fit_obj <- withCallingHandlers(
      tryCatch(
        glmmTMB::glmmTMB(
          formula = cond_formula,
          ziformula = zi_formula,
          family = glmmTMB::nbinom2,
          data = dat_local,
          control = ctrl
        ),
        error = function(e) {
          structure(list(message = e$message), class = "fit_error")
        }
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(obj = fit_obj, warns = unique(warns))
  }

  res1 <- fit_once(glmmTMB::glmmTMBControl())
  fit <- if (inherits(res1$obj, "fit_error")) NULL else res1$obj
  warns <- res1$warns
  fail_reason <- NULL

  if (is.null(fit) || !is_fit_ok(fit, warns)) {
    res2 <- fit_once(glmmTMB::glmmTMBControl(
      optimizer = optim,
      optArgs = list(method = "BFGS")
    ))
    fit2 <- if (inherits(res2$obj, "fit_error")) NULL else res2$obj

    if (!is.null(fit2) && is_fit_ok(fit2, res2$warns)) {
      fit <- fit2
      warns <- res2$warns
    } else {
      fit <- NULL
      fail_reason <- if (inherits(res2$obj, "fit_error")) {
        paste0("error: ", res2$obj$message)
      } else {
        "convergence_or_hessian_failed"
      }
      warns <- unique(c(warns, res2$warns))
    }
  }

  if (is.null(fit)) {
    return(list(
      key = key,
      status = "failed",
      cond_terms = cond_terms,
      zi_terms = zi_terms,
      fail_reason = fail_reason,
      warnings = paste(warns, collapse = " | ")
    ))
  }

  ll <- as.numeric(logLik(fit))
  aic <- as.numeric(AIC(fit))
  bic <- as.numeric(BIC(fit))
  aicc <- calc_aicc(fit)

  coef_cond <- extract_coef_table(fit, "cond")
  coef_zi <- extract_coef_table(fit, "zi")

  list(
    key = key,
    status = "ok",
    cond_terms = cond_terms,
    zi_terms = zi_terms,
    logLik = ll,
    AIC = aic,
    BIC = bic,
    AICc = aicc,
    warnings = paste(warns, collapse = " | "),
    coef_tables = list(cond = coef_cond, zi = coef_zi)
  )
}

new_results <- list()
if (length(pending_tasks) > 0) {
  cat(sprintf("开始并行拟合，worker 数: %d\n", parallel_workers))

  cl <- parallel::makeCluster(parallel_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(glmmTMB))
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      "model_dat", "pending_tasks", "fit_one_task", "make_cond_formula", "make_zi_formula",
      "calc_aicc", "is_fit_ok", "extract_coef_table"
    ),
    envir = environment()
  )

  par_out <- parallel::parLapply(
    cl,
    pending_tasks,
    function(task) fit_one_task(task, dat_local = model_dat)
  )

  new_results <- par_out
  names(new_results) <- vapply(new_results, function(x) x$key, character(1))

  cat(sprintf("并行拟合完成：新增结果 %d\n", length(new_results)))
}

all_results <- c(cached_results, new_results)
if (length(all_results) == 0) {
  stop("没有可用拟合结果（缓存为空且新拟合为空）。")
}

if (use_cache) {
  saveRDS(all_results, cache_file)
  saveRDS(
    list(
      signature = cache_signature,
      created_at = as.character(Sys.time()),
      input_file = normalizePath(input_file, winslash = "/", mustWork = FALSE)
    ),
    cache_meta_file
  )
}

res_df <- data.frame(
  key = vapply(all_results, function(x) x$key, character(1)),
  status = vapply(all_results, function(x) x$status, character(1)),
  cond_terms = vapply(all_results, function(x) {
    if (length(x$cond_terms) == 0) "1" else paste(x$cond_terms, collapse = " + ")
  }, character(1)),
  zi_terms = vapply(all_results, function(x) {
    if (length(x$zi_terms) == 0) "1" else paste(x$zi_terms, collapse = " + ")
  }, character(1)),
  AIC = vapply(all_results, function(x) if (!is.null(x$AIC)) x$AIC else NA_real_, numeric(1)),
  AICc = vapply(all_results, function(x) if (!is.null(x$AICc)) x$AICc else NA_real_, numeric(1)),
  BIC = vapply(all_results, function(x) if (!is.null(x$BIC)) x$BIC else NA_real_, numeric(1)),
  logLik = vapply(all_results, function(x) if (!is.null(x$logLik)) x$logLik else NA_real_, numeric(1)),
  fail_reason = vapply(all_results, function(x) if (!is.null(x$fail_reason)) x$fail_reason else "", character(1)),
  warnings = vapply(all_results, function(x) if (!is.null(x$warnings)) x$warnings else "", character(1)),
  stringsAsFactors = FALSE
)

safe_write_csv(res_df, "all_subsets_model_results_stage3_binary.csv", row.names = FALSE)

failed_df <- res_df[res_df$status == "failed", , drop = FALSE]
if (nrow(failed_df) > 0) {
  safe_write_csv(failed_df, "all_subsets_failed_models_stage3_binary.csv", row.names = FALSE)
}

ok_df <- res_df[res_df$status == "ok" & is.finite(res_df$AICc), , drop = FALSE]
if (nrow(ok_df) == 0) {
  stop("没有稳定收敛且可计算AICc的模型。")
}

ok_df <- ok_df[order(ok_df$AICc), , drop = FALSE]
best_aicc <- ok_df$AICc[1]
ok_df$delta_AICc <- ok_df$AICc - best_aicc
ok_df$rel_lik <- exp(-0.5 * ok_df$delta_AICc)
ok_df$akaike_weight <- ok_df$rel_lik / sum(ok_df$rel_lik)

safe_write_csv(ok_df, "all_subsets_ok_models_ranked_by_AICc_stage3_binary.csv", row.names = FALSE)

selected_df <- ok_df[ok_df$delta_AICc <= aicc_delta_cutoff, , drop = FALSE]
if (nrow(selected_df) == 0) {
  selected_df <- ok_df[1, , drop = FALSE]
}

selected_df$selected_weight <- selected_df$akaike_weight / sum(selected_df$akaike_weight)
safe_write_csv(selected_df, "selected_models_deltaAICc_le_2_stage3_binary.csv", row.names = FALSE)

ok_by_weight <- ok_df[order(-ok_df$akaike_weight, ok_df$AICc), , drop = FALSE]
ok_by_weight$cum_akaike_weight <- cumsum(ok_by_weight$akaike_weight)
cut_idx <- which(ok_by_weight$cum_akaike_weight >= cum_weight_cutoff)[1]
if (is.na(cut_idx)) cut_idx <- nrow(ok_by_weight)

selected95_df <- ok_by_weight[seq_len(cut_idx), , drop = FALSE]
selected95_df$selected_weight95 <- selected95_df$akaike_weight / sum(selected95_df$akaike_weight)
safe_write_csv(selected95_df, "selected_models_cumweight95_stage3_binary.csv", row.names = FALSE)

selected_keys <- selected95_df$key
selected_results <- all_results[selected_keys]

for (i in seq_along(selected_results)) {
  selected_results[[i]]$weight <- selected95_df$selected_weight95[i]
}

avg_cond <- model_average_part(selected_results, "cond")
avg_zi <- model_average_part(selected_results, "zi")

avg_all <- data.table::rbindlist(list(avg_cond, avg_zi), fill = TRUE)
if (nrow(avg_all) > 0) {
  safe_write_csv(avg_all, "model_averaged_coefficients_stage3_binary.csv", row.names = FALSE)
}

term_importance <- function(ok_models, part = c("cond", "zi")) {
  part <- match.arg(part)
  terms_col <- if (part == "cond") "cond_terms" else "zi_terms"
  term_list <- strsplit(ok_models[[terms_col]], " \\+ ", perl = TRUE)
  term_list <- lapply(term_list, function(x) {
    if (length(x) == 1 && x[1] == "1") character(0) else x
  })

  all_terms <- sort(unique(unlist(term_list)))
  if (length(all_terms) == 0) {
    return(data.frame(
      term = character(0),
      importance_weight = numeric(0),
      part = character(0),
      stringsAsFactors = FALSE
    ))
  }

  imp <- vapply(all_terms, function(tm) {
    has_tm <- vapply(term_list, function(v) tm %in% v, logical(1))
    sum(ok_models$akaike_weight[has_tm])
  }, numeric(1))

  out <- data.frame(
    term = all_terms,
    importance_weight = as.numeric(imp),
    part = part,
    stringsAsFactors = FALSE
  )
  out[order(-out$importance_weight, out$term), , drop = FALSE]
}

imp_cond <- term_importance(ok_df, "cond")
imp_zi <- term_importance(ok_df, "zi")
imp_all <- data.table::rbindlist(list(imp_cond, imp_zi), fill = TRUE)
if (nrow(imp_all) > 0) {
  safe_write_csv(imp_all, "variable_importance_akaike_weight_stage3_binary.csv", row.names = FALSE)
}

best_row <- ok_df[1, , drop = FALSE]
best_cond_terms <- strsplit(best_row$cond_terms, " \\+ ", perl = TRUE)[[1]]
best_zi_terms <- strsplit(best_row$zi_terms, " \\+ ", perl = TRUE)[[1]]
if (length(best_cond_terms) == 1 && best_cond_terms == "1") best_cond_terms <- character(0)
if (length(best_zi_terms) == 1 && best_zi_terms == "1") best_zi_terms <- character(0)

best_fit <- glmmTMB::glmmTMB(
  formula = make_cond_formula(best_cond_terms),
  ziformula = make_zi_formula(best_zi_terms),
  family = glmmTMB::nbinom2,
  data = model_dat
)

safe_capture_output(summary(best_fit), "summary_best_model_all_subsets_binary.txt")

zi_link_pred <- as.numeric(predict(best_fit, type = "zlink"))
zi_prob_pred <- as.numeric(predict(best_fit, type = "zprob"))
mu_cond_pred <- as.numeric(predict(best_fit, type = "conditional"))
y_pred_mean <- as.numeric(predict(best_fit, type = "response"))

pred_df <- cbind(
  model_dat,
  data.frame(
    zi_link_pred = zi_link_pred,
    zi_prob_pred = zi_prob_pred,
    mu_cond_pred = mu_cond_pred,
    y_pred_mean = y_pred_mean,
    stringsAsFactors = FALSE
  )
)
safe_write_csv(pred_df, "best_model_predictions_with_zi_prob_stage3_binary.csv", row.names = FALSE)

ic_tab <- data.frame(
  model = "best_by_AICc",
  AIC = as.numeric(AIC(best_fit)),
  AICc = calc_aicc(best_fit),
  BIC = as.numeric(BIC(best_fit)),
  logLik = as.numeric(logLik(best_fit)),
  cond_terms = best_row$cond_terms,
  zi_terms = best_row$zi_terms,
  stringsAsFactors = FALSE
)
safe_write_csv(ic_tab, "best_model_ic_stage3_binary_all_subsets.csv", row.names = FALSE)

if (requireNamespace("DHARMa", quietly = TRUE)) {
  sim_res <- DHARMa::simulateResiduals(best_fit, n = 250)
  tryCatch(
    {
      png(out_file("diagnostic_DHARMa_best_model_binary_all_subsets.png"), width = 1600, height = 1200, res = 180)
      plot(sim_res)
      dev.off()
    },
    error = function(e) {
      message(sprintf("写入失败 [diagnostic_DHARMa_best_model_binary_all_subsets.png]: %s", e$message))
    }
  )

  safe_capture_output(DHARMa::testDispersion(sim_res), "diagnostic_dispersion_test_binary_all_subsets.txt")
  safe_capture_output(DHARMa::testZeroInflation(sim_res), "diagnostic_zeroinflation_test_binary_all_subsets.txt")
}

summary_lines <- c(
  "Stage3 binary 全子集 + 模型平均 完成",
  sprintf("输出目录: %s", output_dir),
  sprintf("理论组合数: %d", length(cond_subsets) * length(zi_subsets)),
  sprintf("成功模型数: %d", nrow(ok_df)),
  sprintf("失败模型数: %d", nrow(failed_df)),
  sprintf("ΔAICc <= %.1f 入选模型数: %d", aicc_delta_cutoff, nrow(selected_df)),
  sprintf("累计Akaike权重达到 %.0f%% 的入选模型数: %d", cum_weight_cutoff * 100, nrow(selected95_df)),
  sprintf("95%%模型集累计原始权重: %.6f", sum(selected95_df$akaike_weight)),
  sprintf("最优模型 cond: %s", best_row$cond_terms),
  sprintf("最优模型 zi: %s", best_row$zi_terms),
  sprintf("最优模型 AICc: %.6f", best_row$AICc)
)

safe_write_lines(summary_lines, "run_summary_stage3_binary_all_subsets.txt")

cat("--------------------------------------------------\n")
cat("Stage3 binary 全子集 + 模型平均 完成。\n")
cat(sprintf("输出目录: %s\n", output_dir))
cat(sprintf("理论组合数: %d\n", length(cond_subsets) * length(zi_subsets)))
cat(sprintf("成功模型数: %d; 失败模型数: %d\n", nrow(ok_df), nrow(failed_df)))
cat(sprintf("ΔAICc <= %.1f 入选模型数: %d\n", aicc_delta_cutoff, nrow(selected_df)))
cat(sprintf("累计Akaike权重达到 %.0f%% 的入选模型数: %d\n", cum_weight_cutoff * 100, nrow(selected95_df)))
cat(sprintf("95%%模型集累计原始权重: %.6f\n", sum(selected95_df$akaike_weight)))
cat("最优模型 IC:\n")
print(ic_tab)
cat("--------------------------------------------------\n")
