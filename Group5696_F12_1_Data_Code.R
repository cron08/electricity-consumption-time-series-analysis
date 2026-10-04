# MATH5845 Group Project Electricity Group: T200

#NOTE ABOUT REPRODUCIBILITY: line 10 should be changed to the appropriate path name.

# Packages and file paths
options(error = NULL)
options(stringsAsFactors = FALSE)
set.seed(5845)

projectdir <- "CHANGE THIS TO OUR GROUP'S PROJECT FOLDER PATH NAME"
data_path <- file.path(projectdir, "Group5696_F12_1_Data.tsf")

# Output folder

runstamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
rundir <- file.path(projectdir, paste0("Run_", runstamp))


if (dir.exists(rundir)) {
  suffix <- 1L
  repeat {
    candidate <- paste0(rundir, "_", suffix)
    if (!dir.exists(candidate)) {
      rundir <- candidate
      break
    }
    suffix <- suffix + 1L
  }
}

seriesid <- "T200"
seasonalperiod <- 52L
start_date <- as.Date("2012-01-01")
lambda_feedback <- 1.999927

# Range of lambda values used for the joint likelihood


joint_lambda_bounds <- c(-1.5, 2.5)
joint_lambda_grid_points <- 81L
joint_lambda_ci_level <- 0.95

# Output folders
figuredir <- file.path(rundir, "Figures")
topmodelfiguredir <- file.path(figuredir, "Top_Models")
tabledir <- file.path(rundir, "Tables")
textdir <- file.path(rundir, "Model_Output")

for (d in c(
  figuredir,
  topmodelfiguredir,
  tabledir,
  textdir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Packages
required_packages <- c("astsa", "forecast", "tseries")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(astsa)
library(forecast)
library(tseries)

cat("Project folder:", projectdir, "\n")
cat("This run will be saved to:", rundir, "\n")
cat("Figures will be saved to:", figuredir, "\n")


# Read T200 series
if (!file.exists(data_path)) {
  stop("Dataset not found at: ", data_path)
}

read_tsf_base <- function(path) {
  lines <- readLines(path, warn = FALSE)
  data_marker <- which(trimws(lines) == "@data")

  if (length(data_marker) != 1L) {
    stop("Could not identify a unique @data line in the .tsf file.")
  }

  data_lines <- lines[(data_marker + 1L):length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]

  parsed <- lapply(data_lines, function(z) {
    pieces <- strsplit(z, ":", fixed = TRUE)[[1]]
    if (length(pieces) < 3L) return(NULL)

    values <- as.numeric(
      strsplit(pieces[length(pieces)], ",", fixed = TRUE)[[1]]
    )

    data.frame(
      series_name = pieces[1],
      start_timestamp = paste(
        pieces[2:(length(pieces) - 1L)],
        collapse = ":"
      ),
      t = seq_along(values),
      value = values,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, parsed)
}

electricity <- read_tsf_base(data_path)
df <- subset(electricity, series_name == seriesid)

if (nrow(df) == 0L) stop("T200 was not found in the dataset.")
if (anyNA(df$value)) stop("T200 contains missing values.")

df$date <- start_date + 7L * (df$t - 1L)
y <- ts(df$value, start = c(2012, 1), frequency = seasonalperiod)

cat("Series:", seriesid, "\n")
cat("Observations:", length(y), "\n")
cat("Start date:", as.character(start_date), "\n")
cat("Frequency:", frequency(y), "\n")

write.csv(
  df,
  file.path(tabledir, "T200_extracted.csv"),
  row.names = FALSE
)


# Functions

fitarima <- function(
    x,
    order,
    seasonal = list(order = c(0, 0, 0), period = seasonalperiod),
    include.mean = FALSE,
    method = "ML"
) {
  tryCatch(
    stats::arima(
      x,
      order = order,
      seasonal = seasonal,
      include.mean = include.mean,
      method = method
    ),
    error = function(e) {
      structure(
        list(error = conditionMessage(e)),
        class = "fit_error"
      )
    }
  )
}

fitfailed <- function(x) inherits(x, "fit_error")

# For D = 1, remove the first 52 residuals before the
# residual diagnostics, following the lecturer feedback.
clean_residuals <- function(fit, D, s = seasonalperiod) {
  if (fitfailed(fit)) return(numeric(0))

  r <- as.numeric(residuals(fit))

  if (D == 1L && length(r) > s) {
    r <- r[-seq_len(s)]
  }

  r[is.finite(r)]
}


nafterdifferencing <- function(n, d, D, s = seasonalperiod) {
  n - d - D * s
}

information_criteria <- function(fit, n_eff) {
  if (fitfailed(fit)) {
    return(c(AIC = NA_real_, AICc = NA_real_, BIC = NA_real_))
  }

  k <- length(coef(fit)) + 1L  # fitted coefficients + white noise variance
  ll <- as.numeric(logLik(fit))
  aic <- -2 * ll + 2 * k

  aicc <- if (n_eff > k + 1L) {
    aic + (2 * k * (k + 1)) / (n_eff - k - 1)
  } else {
    NA_real_
  }

  bic <- -2 * ll + k * log(n_eff)

  c(AIC = aic, AICc = aicc, BIC = bic)
}

parameter_table <- function(model_name, fit) {
  if (fitfailed(fit)) return(NULL)

  est <- coef(fit)
  if (length(est) == 0L) return(NULL)

  se <- if (is.null(fit$var.coef) || length(fit$var.coef) == 0L) {
    rep(NA_real_, length(est))
  } else {
    sqrt(diag(fit$var.coef))
  }

  z <- est / se

  data.frame(
    Model = model_name,
    Parameter = names(est),
    Estimate = as.numeric(est),
    SE = as.numeric(se),
    z_value = as.numeric(z),
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

allarmacoefficientssignificant <- function(fit) {
  if (fitfailed(fit)) return(NA)

  pars <- parameter_table("temp", fit)
  if (is.null(pars)) return(NA)

  armacoefficients <- pars[
    grepl("^(ar|ma|sar|sma)[0-9]+$", pars$Parameter),
    ,
    drop = FALSE
  ]

  if (nrow(armacoefficients) == 0L) return(NA)
  all(armacoefficients$p_value < 0.05, na.rm = TRUE)
}

# Roots of the AR and MA polynomials

root_minima <- function(fit) {
  if (fitfailed(fit)) {
    return(c(min_AR_root = NA_real_, min_MA_root = NA_real_))
  }

  cf <- coef(fit)

  ar_cf <- cf[grep("^ar[0-9]+$", names(cf))]
  ma_cf <- cf[grep("^ma[0-9]+$", names(cf))]
  sar_cf <- cf[grep("^sar[0-9]+$", names(cf))]
  sma_cf <- cf[grep("^sma[0-9]+$", names(cf))]

  ar_roots <- if (length(ar_cf)) {
    Mod(polyroot(c(1, -ar_cf)))
  } else {
    numeric(0)
  }

  ma_roots <- if (length(ma_cf)) {
    Mod(polyroot(c(1, ma_cf)))
  } else {
    numeric(0)
  }

  
  sar_roots <- if (length(sar_cf)) {
    Mod(polyroot(c(1, -sar_cf)))^(1 / seasonalperiod)
  } else {
    numeric(0)
  }

  sma_roots <- if (length(sma_cf)) {
    Mod(polyroot(c(1, sma_cf)))^(1 / seasonalperiod)
  } else {
    numeric(0)
  }

  all_ar <- c(ar_roots, sar_roots)
  all_ma <- c(ma_roots, sma_roots)

  c(
    min_AR_root = if (length(all_ar)) min(all_ar) else Inf,
    min_MA_root = if (length(all_ma)) min(all_ma) else Inf
  )
}

rolling_mean_sd <- function(x, window = 26L) {
  x <- as.numeric(x)
  if (length(x) < window) stop("Window is longer than the series.")

  endpoints <- seq.int(window, length(x))

  means <- vapply(endpoints, function(i) {
    mean(x[(i - window + 1L):i])
  }, numeric(1))

  sds <- vapply(endpoints, function(i) {
    sd(x[(i - window + 1L):i])
  }, numeric(1))

  data.frame(mean = means, sd = sds)
}


# Joint Box-Cox likelihood


boxcox_log_jacobian <- function(x, lambda, order, seasonal) {
  x <- as.numeric(x)

  if (any(!is.finite(x)) || any(x <= 0)) {
    stop("Box-Cox transformation requires positive finite observations.")
  }

  d <- as.integer(order[2])
  D <- as.integer(seasonal$order[2])
  s <- as.integer(seasonal$period)
  n_drop <- d + D * s

  if (n_drop >= length(x)) {
    stop("Differencing removes all observations.")
  }

  x_used <- x[(n_drop + 1L):length(x)]
  (lambda - 1) * sum(log(x_used))
}


joint_loglik_for_spec <- function(
    x,
    lambda,
    order,
    seasonal,
    include.mean = FALSE,
    method = "ML"
) {
  z <- forecast::BoxCox(x, lambda = lambda)

  fit <- fitarima(
    x = z,
    order = order,
    seasonal = seasonal,
    include.mean = include.mean,
    method = method
  )

  if (fitfailed(fit)) {
    return(list(
      fit = fit,
      transformed_logLik = NA_real_,
      log_jacobian = NA_real_,
      joint_logLik = -Inf
    ))
  }

  ll_transformed <- as.numeric(logLik(fit))
  log_jacobian <- boxcox_log_jacobian(x, lambda, order, seasonal)

  list(
    fit = fit,
    transformed_logLik = ll_transformed,
    log_jacobian = log_jacobian,
    joint_logLik = ll_transformed + log_jacobian
  )
}

profile_likelihood_ci <- function(
    profile_df,
    lambda_hat,
    ll_hat,
    level = 0.95
) {
  cutoff <- ll_hat - 0.5 * qchisq(level, df = 1)

  prof <- profile_df[is.finite(profile_df$JointLogLik), , drop = FALSE]
  prof <- prof[order(prof$Lambda), , drop = FALSE]

  
  prof <- rbind(
    prof,
    data.frame(Lambda = lambda_hat, JointLogLik = ll_hat)
  )
  prof <- prof[order(prof$Lambda), , drop = FALSE]

  
  
  left <- prof[prof$Lambda <= lambda_hat, , drop = FALSE]
  left <- left[order(left$Lambda, decreasing = TRUE), , drop = FALSE]

  lower <- NA_real_
  if (nrow(left) >= 2L) {
    below <- which(left$JointLogLik < cutoff)
    if (length(below) > 0L) {
      j <- below[1L]
      if (j > 1L) {
        x1 <- left$JointLogLik[j]
        x2 <- left$JointLogLik[j - 1L]
        y1 <- left$Lambda[j]
        y2 <- left$Lambda[j - 1L]
        lower <- y1 + (cutoff - x1) * (y2 - y1) / (x2 - x1)
      }
    }
  }

  
  right <- prof[prof$Lambda >= lambda_hat, , drop = FALSE]
  right <- right[order(right$Lambda), , drop = FALSE]

  upper <- NA_real_
  if (nrow(right) >= 2L) {
    below <- which(right$JointLogLik < cutoff)
    if (length(below) > 0L) {
      j <- below[1L]
      if (j > 1L) {
        x1 <- right$JointLogLik[j - 1L]
        x2 <- right$JointLogLik[j]
        y1 <- right$Lambda[j - 1L]
        y2 <- right$Lambda[j]
        upper <- y1 + (cutoff - x1) * (y2 - y1) / (x2 - x1)
      }
    }
  }

  c(
    lower = lower,
    upper = upper,
    cutoff = cutoff
  )
}

fit_joint_boxcox_sarima <- function(
    x,
    order,
    seasonal,
    include.mean = FALSE,
    lambda_bounds = joint_lambda_bounds,
    grid_points = joint_lambda_grid_points,
    ci_level = joint_lambda_ci_level,
    method = "ML"
) {
  if (length(lambda_bounds) != 2L || lambda_bounds[1] >= lambda_bounds[2]) {
    stop("lambda_bounds must contain an increasing lower and upper bound.")
  }

  lambda_grid <- seq(
    lambda_bounds[1],
    lambda_bounds[2],
    length.out = grid_points
  )

  profile_value <- function(lambda) {
    joint_loglik_for_spec(
      x = x,
      lambda = lambda,
      order = order,
      seasonal = seasonal,
      include.mean = include.mean,
      method = method
    )$joint_logLik
  }

  ll_grid <- vapply(lambda_grid, profile_value, numeric(1))

  if (!any(is.finite(ll_grid))) {
    stop("All joint Box-Cox/SARIMA likelihood evaluations failed.")
  }

  best_grid_index <- which.max(ll_grid)

  
  
  if (best_grid_index == 1L || best_grid_index == length(lambda_grid)) {
    lambda_hat <- lambda_grid[best_grid_index]
    boundary_solution <- TRUE
  } else {
    local_interval <- c(
      lambda_grid[best_grid_index - 1L],
      lambda_grid[best_grid_index + 1L]
    )

    refined <- optimize(
      f = function(lambda) {
        ll <- profile_value(lambda)
        if (!is.finite(ll)) return(.Machine$double.xmax / 100)
        -ll
      },
      interval = local_interval,
      tol = 1e-5
    )

    lambda_hat <- refined$minimum
    boundary_solution <- FALSE
  }

  optimum <- joint_loglik_for_spec(
    x = x,
    lambda = lambda_hat,
    order = order,
    seasonal = seasonal,
    include.mean = include.mean,
    method = method
  )

  profile_df <- data.frame(
    Lambda = lambda_grid,
    JointLogLik = ll_grid
  )

  ci <- profile_likelihood_ci(
    profile_df = profile_df,
    lambda_hat = lambda_hat,
    ll_hat = optimum$joint_logLik,
    level = ci_level
  )

  list(
    lambda_hat = lambda_hat,
    lambda_CI_lower = unname(ci["lower"]),
    lambda_CI_upper = unname(ci["upper"]),
    profile_cutoff = unname(ci["cutoff"]),
    fit = optimum$fit,
    transformed_logLik = optimum$transformed_logLik,
    log_jacobian = optimum$log_jacobian,
    joint_logLik = optimum$joint_logLik,
    profile = profile_df,
    boundary_solution = boundary_solution,
    lambda_bounds = lambda_bounds
  )
}

joint_information_criteria <- function(
    joint_logLik,
    fit,
    n_eff,
    lambda_estimated = TRUE
) {
  if (!is.finite(joint_logLik) || fitfailed(fit)) {
    return(c(
      k = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      BIC = NA_real_
    ))
  }

  
  k <- length(coef(fit)) + 1L + as.integer(lambda_estimated)

  aic <- -2 * joint_logLik + 2 * k
  aicc <- if (n_eff > k + 1L) {
    aic + (2 * k * (k + 1)) / (n_eff - k - 1)
  } else {
    NA_real_
  }
  bic <- -2 * joint_logLik + k * log(n_eff)

  c(k = k, AIC = aic, AICc = aicc, BIC = bic)
}


# Exploratory analysis

# Original series

png(
  file.path(figuredir, "01_original_series.png"),
  width = 1800,
  height = 950,
  res = 180,
  bg = "white"
)
plot(
  y,
  type = "l",
  lwd = 1.3,
  main = "T200 Weekly Electricity Consumption",
  xlab = "Year",
  ylab = "Electricity consumption"
)
dev.off()


# Series by week of year

annual_matrix <- matrix(
  as.numeric(y),
  nrow = seasonalperiod,
  ncol = length(y) / seasonalperiod
)

png(
  file.path(figuredir, "02_annual_seasonal_overlay.png"),
  width = 1800,
  height = 1000,
  res = 180,
  bg = "white"
)
matplot(
  1:seasonalperiod,
  annual_matrix,
  type = "l",
  lty = 1:ncol(annual_matrix),
  lwd = 1.3,
  xlab = "Week of year",
  ylab = "Electricity consumption",
  main = "T200: Annual Pattern by Week of Year"
)
legend(
  "topright",
  legend = paste0("Year ", seq_len(ncol(annual_matrix))),
  lty = 1:ncol(annual_matrix),
  bty = "n"
)
dev.off()


# ACF and PACF

png(
  file.path(figuredir, "03_original_acf_pacf.png"),
  width = 1700,
  height = 1000,
  res = 180,
  bg = "white"
)
astsa::acf2(y, max.lag = 104, main = "T200: Original Series")
dev.off()


# Transformation and differencing
# Compare the Guerrero estimate, the likelihood estimate and the lecturer's value.

lambda_guerrero <- forecast::BoxCox.lambda(y, method = "guerrero")
lambda_loglik <- forecast::BoxCox.lambda(y, method = "loglik")

lambda_table <- data.frame(
  Method = c(
    "Guerrero (preliminary)",
    "forecast profile log-likelihood (preliminary)",
    "Lecturer feedback"
  ),
  Lambda = c(lambda_guerrero, lambda_loglik, lambda_feedback)
)

write.csv(
  lambda_table,
  file.path(tabledir, "boxcox_lambda_values.csv"),
  row.names = FALSE
)

print(lambda_table)

y_bc_feedback <- forecast::BoxCox(y, lambda_feedback)

# Box-Cox transformation


raw_mv <- rolling_mean_sd(y, window = 26L)
bc_mv <- rolling_mean_sd(y_bc_feedback, window = 26L)

png(
  file.path(figuredir, "04_boxcox_variance_check.png"),
  width = 1800,
  height = 1500,
  res = 180,
  bg = "white"
)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(
  y,
  type = "l",
  main = "Original series",
  xlab = "Year",
  ylab = "Electricity consumption"
)

plot(
  y_bc_feedback,
  type = "l",
  main = paste0("Preliminary Box-Cox transform (lecturer lambda = ", round(lambda_feedback, 4), ")"),
  xlab = "Year",
  ylab = "Transformed value"
)

plot(
  raw_mv$mean,
  raw_mv$sd,
  pch = 1,
  main = "Original scale: local mean vs local SD",
  xlab = "26-week local mean",
  ylab = "26-week local SD"
)
abline(lm(sd ~ mean, data = raw_mv), lty = 2)

plot(
  bc_mv$mean,
  bc_mv$sd,
  pch = 1,
  main = "Box-Cox scale: local mean vs local SD",
  xlab = "26-week local mean",
  ylab = "26-week local SD"
)
abline(lm(sd ~ mean, data = bc_mv), lty = 2)

par(mfrow = c(1, 1))
dev.off()


# Differenced series
y_diff <- diff(y, lag = 1)
y_sdiff <- diff(y, lag = seasonalperiod)
y_both <- diff(y_sdiff, lag = 1)

# Differenced series plots

png(
  file.path(figuredir, "05_candidate_differences.png"),
  width = 1800,
  height = 1450,
  res = 180,
  bg = "white"
)
par(mfrow = c(3, 1), mar = c(3.5, 4, 3, 1))

plot(
  y_diff,
  main = "Ordinary first difference",
  xlab = "Year",
  ylab = expression((1-B)*Y[t])
)
abline(h = 0, lty = 2)

plot(
  y_sdiff,
  main = "Seasonal difference at lag 52",
  xlab = "Year",
  ylab = expression((1-B^52)*Y[t])
)
abline(h = 0, lty = 2)

plot(
  y_both,
  main = "Ordinary and seasonal differences",
  xlab = "Year",
  ylab = expression((1-B)*(1-B^52)*Y[t])
)
abline(h = 0, lty = 2)

par(mfrow = c(1, 1))
dev.off()


# ACF and PACF after differencing

png(
  file.path(figuredir, "06_ordinary_difference_acf_pacf.png"),
  width = 1700,
  height = 1000,
  res = 180,
  bg = "white"
)
astsa::acf2(y_diff, max.lag = 104, main = "Ordinary Difference")
dev.off()

png(
  file.path(figuredir, "07_seasonal_difference_acf_pacf.png"),
  width = 1700,
  height = 1000,
  res = 180,
  bg = "white"
)
astsa::acf2(y_sdiff, max.lag = 104, main = "Seasonal Difference at Lag 52")
dev.off()

png(
  file.path(figuredir, "08_double_difference_acf_pacf.png"),
  width = 1700,
  height = 1000,
  res = 180,
  bg = "white"
)
astsa::acf2(y_both, max.lag = 80, main = "Ordinary + Seasonal Difference")
dev.off()


# ADF and KPSS tests

d_rec <- forecast::ndiffs(y)
D_rec <- forecast::nsdiffs(y)

stationarity_tests <- data.frame(
  Series = c(
    "Original",
    "Ordinary difference",
    "Seasonal difference",
    "Both differences"
  ),
  ADF_p = c(
    tryCatch(tseries::adf.test(y)$p.value, error = function(e) NA_real_),
    tryCatch(tseries::adf.test(y_diff)$p.value, error = function(e) NA_real_),
    tryCatch(tseries::adf.test(y_sdiff)$p.value, error = function(e) NA_real_),
    tryCatch(tseries::adf.test(y_both)$p.value, error = function(e) NA_real_)
  ),
  KPSS_p = c(
    tryCatch(tseries::kpss.test(y, null = "Level")$p.value, error = function(e) NA_real_),
    tryCatch(tseries::kpss.test(y_diff, null = "Level")$p.value, error = function(e) NA_real_),
    tryCatch(tseries::kpss.test(y_sdiff, null = "Level")$p.value, error = function(e) NA_real_),
    tryCatch(tseries::kpss.test(y_both, null = "Level")$p.value, error = function(e) NA_real_)
  )
)

write.csv(
  stationarity_tests,
  file.path(tabledir, "stationarity_tests.csv"),
  row.names = FALSE
)

cat("ndiffs(y) =", d_rec, "\n")
cat("nsdiffs(y) =", D_rec, "\n")
print(stationarity_tests)


# ARMA, ARIMA and SARIMA models


modelcandidates <- data.frame(
  Model = c(
    # ARMA models
    "AR(1)",
    "AR(2)",
    "MA(1)",
    "MA(2)",
    "ARMA(1,1)",
    "ARMA(2,1)",
    "ARMA(1,2)",

    # ARIMA models
    "ARIMA(0,1,0)",
    "ARIMA(1,1,0)",
    "ARIMA(0,1,1)",
    "ARIMA(2,1,0)",
    "ARIMA(0,1,2)",
    "ARIMA(1,1,1)",
    "ARIMA(2,1,1)",
    "ARIMA(0,1,3)",

    # SARIMA models with d = 0, D = 1
    "SARIMA(1,0,0)(0,1,0)[52]",
    "SARIMA(0,0,1)(0,1,0)[52]",
    "SARIMA(1,0,1)(0,1,0)[52]",
    "SARIMA(1,0,0)(0,1,1)[52]",
    "SARIMA(0,0,1)(0,1,1)[52]",
    "SARIMA(1,0,1)(0,1,1)[52]",
    "SARIMA(0,0,3)(1,1,0)[52]",
    "SARIMA(1,0,0)(1,1,0)[52]",

    # SARIMA models with d = 1, D = 1
    "SARIMA(0,1,1)(0,1,1)[52]",
    "SARIMA(1,1,0)(0,1,1)[52]",
    "SARIMA(1,1,1)(0,1,1)[52]",
    "SARIMA(0,1,3)(1,1,0)[52]",
    "SARIMA(2,1,1)(0,1,0)[52]",
    "SARIMA(2,1,1)(1,1,0)[52]"
  ),

  Family = c(
    rep("ARMA on original", 7),
    rep("Nonseasonal ARIMA", 8),
    rep("SARIMA d=0,D=1", 8),
    rep("SARIMA d=1,D=1", 6)
  ),

  p = c(
    1,2,0,0,1,2,1,
    0,1,0,2,0,1,2,0,
    1,0,1,1,0,1,0,1,
    0,1,1,0,2,2
  ),

  d = c(
    rep(0, 7),
    rep(1, 8),
    rep(0, 8),
    rep(1, 6)
  ),

  q = c(
    0,0,1,2,1,1,2,
    0,0,1,0,2,1,1,3,
    0,1,1,0,1,1,3,0,
    1,0,1,3,1,1
  ),

  P = c(
    rep(0, 15),
    0,0,0,0,0,0,1,1,
    0,0,0,1,0,1
  ),

  D = c(
    rep(0, 15),
    rep(1, 14)
  ),

  Q = c(
    rep(0, 15),
    0,0,0,1,1,1,0,0,
    1,1,1,0,0,0
  ),

  stringsAsFactors = FALSE
)

fit_spec <- function(spec_row, x = y) {
  include_mean <- (spec_row$d == 0L && spec_row$D == 0L)

  fitarima(
    x = x,
    order = c(spec_row$p, spec_row$d, spec_row$q),
    seasonal = list(
      order = c(spec_row$P, spec_row$D, spec_row$Q),
      period = seasonalperiod
    ),
    include.mean = include_mean,
    method = "ML"
  )
}

candidate_fits <- setNames(
  lapply(seq_len(nrow(modelcandidates)), function(i) {
    fit_spec(modelcandidates[i, , drop = FALSE])
  }),
  modelcandidates$Model
)

candidate_summary_one <- function(i) {
  modelcandidate <- modelcandidates[i, , drop = FALSE]
  fit <- candidate_fits[[modelcandidate$Model]]

  if (fitfailed(fit)) {
    return(data.frame(
      Model = modelcandidate$Model,
      Family = modelcandidate$Family,
      AIC = NA_real_,
      AICc = NA_real_,
      BIC = NA_real_,
      Fit_status = fit$error,
      stringsAsFactors = FALSE
    ))
  }

  n_eff <- nafterdifferencing(length(y), modelcandidate$d, modelcandidate$D)
  ic <- information_criteria(fit, n_eff)

  data.frame(
    Model = modelcandidate$Model,
    Family = modelcandidate$Family,
    AIC = unname(ic["AIC"]),
    AICc = unname(ic["AICc"]),
    BIC = unname(ic["BIC"]),
    Fit_status = "OK",
    stringsAsFactors = FALSE
  )
}

candidate_comparison <- do.call(
  rbind,
  lapply(seq_len(nrow(modelcandidates)), candidate_summary_one)
)

candidate_comparison <- candidate_comparison[
  order(candidate_comparison$AICc),
]

write.csv(
  candidate_comparison,
  file.path(tabledir, "candidate_model_comparison.csv"),
  row.names = FALSE
)

print(
  candidate_comparison[, c("Model", "Family", "AIC", "AICc", "BIC")],
  row.names = FALSE
)


# auto.arima

cat("\nrunning auto.arima models for comparison...\n")

auto_raw <- forecast::auto.arima(
  y,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE,
  trace = FALSE
)

auto_boxcox <- forecast::auto.arima(
  y,
  seasonal = TRUE,
  lambda = lambda_feedback,
  stepwise = FALSE,
  approximation = FALSE,
  trace = FALSE
)

capture.output(
  list(
    raw = auto_raw,
    boxcox_lambda_1_999927 = auto_boxcox
  ),
  file = file.path(textdir, "auto_arima_guides.txt")
)

cat("Raw auto.arima guide:\n")
print(auto_raw)
cat("\nbox-Cox auto.arima guide:\n")
print(auto_boxcox)


# Joint Box-Cox likelihood


seasonal_010 <- list(order = c(0, 1, 0), period = seasonalperiod)
seasonal_110 <- list(order = c(1, 1, 0), period = seasonalperiod)

cat("\njoint Joint Box-Cox likelihood...\n")

joint_bc_211_010 <- fit_joint_boxcox_sarima(
  x = y,
  order = c(2, 1, 1),
  seasonal = seasonal_010,
  include.mean = FALSE
)

joint_bc_013_110 <- fit_joint_boxcox_sarima(
  x = y,
  order = c(0, 1, 3),
  seasonal = seasonal_110,
  include.mean = FALSE
)

joint_bc_211_110 <- fit_joint_boxcox_sarima(
  x = y,
  order = c(2, 1, 1),
  seasonal = seasonal_110,
  include.mean = FALSE
)


raw_joint_013_110 <- joint_loglik_for_spec(
  y,
  lambda = 1,
  order = c(0, 1, 3),
  seasonal = seasonal_110,
  include.mean = FALSE
)

raw_joint_211_110 <- joint_loglik_for_spec(
  y,
  lambda = 1,
  order = c(2, 1, 1),
  seasonal = seasonal_110,
  include.mean = FALSE
)

fit_raw_013_110 <- raw_joint_013_110$fit
fit_raw_211_110 <- raw_joint_211_110$fit


fit_joint_211_010 <- joint_bc_211_010$fit
fit_joint_013_110 <- joint_bc_013_110$fit
fit_joint_211_110 <- joint_bc_211_110$fit


joint_boxcox_estimates <- data.frame(
  Structure = c(
    "SARIMA(2,1,1)(0,1,0)[52]",
    "SARIMA(0,1,3)(1,1,0)[52]",
    "SARIMA(2,1,1)(1,1,0)[52]"
  ),
  Lambda_joint_MLE = c(
    joint_bc_211_010$lambda_hat,
    joint_bc_013_110$lambda_hat,
    joint_bc_211_110$lambda_hat
  ),
  Lambda_CI_lower = c(
    joint_bc_211_010$lambda_CI_lower,
    joint_bc_013_110$lambda_CI_lower,
    joint_bc_211_110$lambda_CI_lower
  ),
  Lambda_CI_upper = c(
    joint_bc_211_010$lambda_CI_upper,
    joint_bc_013_110$lambda_CI_upper,
    joint_bc_211_110$lambda_CI_upper
  ),
  Joint_logLik_max = c(
    joint_bc_211_010$joint_logLik,
    joint_bc_013_110$joint_logLik,
    joint_bc_211_110$joint_logLik
  ),
  Joint_logLik_at_lambda_1 = c(
    joint_loglik_for_spec(y, 1, c(2,1,1), seasonal_010)$joint_logLik,
    joint_loglik_for_spec(y, 1, c(0,1,3), seasonal_110)$joint_logLik,
    joint_loglik_for_spec(y, 1, c(2,1,1), seasonal_110)$joint_logLik
  ),
  Joint_logLik_at_lecturer_lambda = c(
    joint_loglik_for_spec(y, lambda_feedback, c(2,1,1), seasonal_010)$joint_logLik,
    joint_loglik_for_spec(y, lambda_feedback, c(0,1,3), seasonal_110)$joint_logLik,
    joint_loglik_for_spec(y, lambda_feedback, c(2,1,1), seasonal_110)$joint_logLik
  ),
  Boundary_solution = c(
    joint_bc_211_010$boundary_solution,
    joint_bc_013_110$boundary_solution,
    joint_bc_211_110$boundary_solution
  ),
  stringsAsFactors = FALSE
)

joint_boxcox_estimates$Improvement_over_lecturer_lambda <-
  joint_boxcox_estimates$Joint_logLik_max -
  joint_boxcox_estimates$Joint_logLik_at_lecturer_lambda

write.csv(
  joint_boxcox_estimates,
  file.path(tabledir, "joint_boxcox_estimates.csv"),
  row.names = FALSE
)

print(joint_boxcox_estimates, digits = 5, row.names = FALSE)


if (any(joint_boxcox_estimates$Boundary_solution)) {
  warning(
    "At least one joint Box-Cox estimate is on the lambda search boundary. ",
    "Expand joint_lambda_bounds and rerun before interpreting that estimate."
  )
}

profile_211_010 <- transform(
  joint_bc_211_010$profile,
  Structure = "SARIMA(2,1,1)(0,1,0)[52]"
)
profile_013_110 <- transform(
  joint_bc_013_110$profile,
  Structure = "SARIMA(0,1,3)(1,1,0)[52]"
)
profile_211_110 <- transform(
  joint_bc_211_110$profile,
  Structure = "SARIMA(2,1,1)(1,1,0)[52]"
)

joint_profile_table <- rbind(
  profile_211_010,
  profile_013_110,
  profile_211_110
)

write.csv(
  joint_profile_table,
  file.path(tabledir, "joint_boxcox_profile_likelihoods.csv"),
  row.names = FALSE
)

# Box-Cox likelihood over lambda

profile_list <- list(
  "SARIMA(2,1,1)(0,1,0)[52]" = joint_bc_211_010,
  "SARIMA(0,1,3)(1,1,0)[52]" = joint_bc_013_110,
  "SARIMA(2,1,1)(1,1,0)[52]" = joint_bc_211_110
)

all_profile_ll <- unlist(lapply(profile_list, function(obj) {
  obj$profile$JointLogLik[is.finite(obj$profile$JointLogLik)]
}))

png(
  file.path(figuredir, "09_joint_boxcox_profile_likelihood.png"),
  width = 1900,
  height = 1100,
  res = 180,
  bg = "white"
)

plot(
  x = joint_lambda_bounds,
  y = range(all_profile_ll),
  type = "n",
  xlab = expression(lambda),
  ylab = "Joint log-likelihood (original-data scale)",
  main = "Joint Box-Cox + SARIMA Likelihood over Lambda"
)

for (i in seq_along(profile_list)) {
  obj <- profile_list[[i]]
  lines(
    obj$profile$Lambda,
    obj$profile$JointLogLik,
    lty = i,
    lwd = 1.4
  )
  points(
    obj$lambda_hat,
    obj$joint_logLik,
    pch = i
  )
}

abline(v = 1, lty = 4)
abline(v = lambda_feedback, lty = 5)

legend(
  "bottomright",
  legend = c(
    names(profile_list),
    "lambda = 1 (raw/no transform)",
    paste0("lecturer lambda = ", round(lambda_feedback, 4))
  ),
  lty = c(1, 2, 3, 4, 5),
  bty = "n",
  cex = 0.85
)

dev.off()



# Guerrero Box-Cox models

lambda_guerrero <- BoxCox.lambda(y, method = "guerrero")
y_guerrero <- BoxCox(y, lambda_guerrero)

guerrero_013_110 <- Arima(
  y_guerrero,
  order = c(0, 1, 3),
  seasonal = list(order = c(1, 1, 0), period = seasonalperiod),
  include.constant = FALSE,
  method = "ML"
)

guerrero_211_110 <- Arima(
  y_guerrero,
  order = c(2, 1, 1),
  seasonal = list(order = c(1, 1, 0), period = seasonalperiod),
  include.constant = FALSE,
  method = "ML"
)

guerrero_211_010 <- Arima(
  y_guerrero,
  order = c(2, 1, 1),
  seasonal = list(order = c(0, 1, 0), period = seasonalperiod),
  include.constant = FALSE,
  method = "ML"
)

# Top models

topmodels <- list(
  "SARIMA(0,1,3)(1,1,0)[52], raw" = fit_raw_013_110,
  "SARIMA(2,1,1)(1,1,0)[52], raw" = fit_raw_211_110,
  "SARIMA(0,1,3)(1,1,0)[52], Guerrero Box-Cox" = guerrero_013_110,
  "SARIMA(2,1,1)(1,1,0)[52], Guerrero Box-Cox" = guerrero_211_110,
  "SARIMA(2,1,1)(0,1,0)[52], Guerrero Box-Cox" = guerrero_211_010
)

topmodelnames <- names(topmodels)

topmodeltransform <- c(
  "Raw",
  "Raw",
  "Guerrero Box-Cox",
  "Guerrero Box-Cox",
  "Guerrero Box-Cox"
)

topmodellambda <- c(
  1,
  1,
  lambda_guerrero,
  lambda_guerrero,
  lambda_guerrero
)

topmodellambdaestimated <- c(
  FALSE,
  FALSE,
  TRUE,
  TRUE,
  TRUE
)

topmodelorder <- list(
  c(0, 1, 3),
  c(2, 1, 1),
  c(0, 1, 3),
  c(2, 1, 1),
  c(2, 1, 1)
)

topmodelseasonal <- list(
  list(order = c(1, 1, 0), period = seasonalperiod),
  list(order = c(1, 1, 0), period = seasonalperiod),
  list(order = c(1, 1, 0), period = seasonalperiod),
  list(order = c(1, 1, 0), period = seasonalperiod),
  list(order = c(0, 1, 0), period = seasonalperiod)
)

topmodelfitdf <- c(4, 4, 4, 4, 3)

topmodelsummary <- function(
    fit,
    modelname,
    transformation,
    lambda,
    lambdaestimated,
    order,
    seasonal,
    fitdf
) {

  D <- seasonal$order[2]
  r <- clean_residuals(fit, D)

  stdres <- (r - mean(r)) / sd(r)

  d <- order[2]
  s <- seasonal$period
  n <- length(y) - d - D * s

  loglik <- as.numeric(logLik(fit))

  if (lambda != 1) {
    loglik <- loglik + boxcox_log_jacobian(
      y,
      lambda,
      order,
      seasonal
    )
  }

  k <- length(coef(fit)) + 1 + as.integer(lambdaestimated)

  aic <- -2 * loglik + 2 * k
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  bic <- -2 * loglik + k * log(n)

  lb20 <- Box.test(
    stdres,
    lag = min(20, length(stdres) - 1),
    type = "Ljung-Box",
    fitdf = fitdf
  )

  lb31 <- Box.test(
    stdres,
    lag = min(31, length(stdres) - 1),
    type = "Ljung-Box",
    fitdf = fitdf
  )

  sw <- shapiro.test(stdres)

  estimates <- coef(fit)
  ses <- sqrt(diag(fit$var.coef))
  pvalues <- 2 * pnorm(
    abs(estimates / ses),
    lower.tail = FALSE
  )

  roots <- root_minima(fit)

  data.frame(
    Model = modelname,
    Transformation = transformation,
    Lambda = lambda,
    Lambda_estimated = lambdaestimated,
    n_eff = n,
    k = k,
    adjusted_logLik = loglik,
    adjusted_AIC = aic,
    adjusted_AICc = aicc,
    adjusted_BIC = bic,
    LjungBox20_p = lb20$p.value,
    LjungBox31_p = lb31$p.value,
    Shapiro_W = unname(sw$statistic),
    Shapiro_p = sw$p.value,
    min_AR_root = unname(roots["min_AR_root"]),
    min_MA_root = unname(roots["min_MA_root"]),
    All_AR_MA_coefficients_significant_5pct =
      all(pvalues < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

topmodeltable <- do.call(
  rbind,
  Map(
    topmodelsummary,
    topmodels,
    topmodelnames,
    topmodeltransform,
    topmodellambda,
    topmodellambdaestimated,
    topmodelorder,
    topmodelseasonal,
    topmodelfitdf
  )
)

topmodeltable$Delta_AICc <-
  topmodeltable$adjusted_AICc - min(topmodeltable$adjusted_AICc)

topmodeltable$Delta_BIC <-
  topmodeltable$adjusted_BIC - min(topmodeltable$adjusted_BIC)

topmodeltable <- topmodeltable[
  order(topmodeltable$adjusted_AICc),
]

print(topmodeltable, digits = 4, row.names = FALSE)

write.csv(
  topmodeltable,
  file.path(tabledir, "top_model_comparison.csv"),
  row.names = FALSE
)

topmodelparameters <- do.call(
  rbind,
  lapply(seq_along(topmodels), function(i) {

    fit <- topmodels[[i]]
    estimates <- coef(fit)
    ses <- sqrt(diag(fit$var.coef))

    data.frame(
      Model = topmodelnames[i],
      Transformation = topmodeltransform[i],
      Lambda = topmodellambda[i],
      Coefficient = names(estimates),
      Estimate = as.numeric(estimates),
      SE = as.numeric(ses),
      z = as.numeric(estimates / ses),
      p_value = 2 * pnorm(
        abs(estimates / ses),
        lower.tail = FALSE
      ),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  topmodelparameters,
  file.path(tabledir, "top_model_parameters.csv"),
  row.names = FALSE
)

# Residual diagnostics


topmodeldiagnosticplot <- function(
    fit,
    model_name,
    D,
    fitdf,
    max_lag = 52L,
    qq_simulations = 5000L,
    seed = 5845L
) {
  if (fitfailed(fit)) stop("Cannot diagnose failed fit: ", model_name)

  r <- clean_residuals(fit, D)
  std_r <- (r - mean(r)) / sd(r)
  n <- length(std_r)

  
  sw <- shapiro.test(std_r)

  
  last_lag <- min(max_lag, n - 1L)
  first_lag <- max(fitdf + 1L, 1L)
  lags <- seq.int(first_lag, last_lag)

  lb_p <- vapply(lags, function(h) {
    Box.test(
      std_r,
      lag = h,
      type = "Ljung-Box",
      fitdf = fitdf
    )$p.value
  }, numeric(1))

  
  theoretical <- qnorm(ppoints(n))
  observed <- sort(std_r)

  
  
  set.seed(seed)
  qq_sim <- replicate(qq_simulations, {
    z <- rnorm(n)
    z <- (z - mean(z)) / sd(z)
    sort(z)
  })

  qq_lower <- apply(qq_sim, 1, quantile, probs = 0.025)
  qq_upper <- apply(qq_sim, 1, quantile, probs = 0.975)

  
  
  layout(
    matrix(
      c(
        1, 1, 1,
        2, 3, 4,
        5, 5, 5
      ),
      nrow = 3,
      byrow = TRUE
    ),
    heights = c(1, 1.2, 1)
  )

  par(mar = c(4, 4, 3, 1))

  
  plot(
    std_r,
    type = "l",
    main = paste("Standardised Residuals:", model_name),
    xlab = "Residual index",
    ylab = "Standardized residual"
  )
  abline(h = 0, lty = 2)

  
  acf(
    std_r,
    lag.max = last_lag,
    main = "ACF of Residuals"
  )

  
  pacf(
    std_r,
    lag.max = last_lag,
    main = "PACF of Residuals"
  )

  
  plot(
    theoretical,
    observed,
    pch = 1,
    xlab = "Theoretical quantiles",
    ylab = "Sample quantiles",
    main = "Normal Q-Q Plot"
  )
  lines(theoretical, qq_lower, lty = 2)
  lines(theoretical, qq_upper, lty = 2)
  qqline(std_r)

  mtext(
    paste0(
      "95% pointwise envelope; Shapiro-Wilk p = ",
      format.pval(sw$p.value, digits = 3, eps = 0.001)
    ),
    side = 3,
    line = 0.15,
    cex = 0.72
  )

  
  plot(
    lags,
    lb_p,
    type = "o",
    pch = 1,
    ylim = c(0, 1),
    xlab = "Lag (H)",
    ylab = "p-value",
    main = "p-values for Ljung-Box Statistic"
  )
  abline(h = 0.05, lty = 2)
  abline(v = seasonalperiod, lty = 3)

  invisible(
    list(
      residuals = std_r,
      Shapiro_W = unname(sw$statistic),
      Shapiro_p = sw$p.value,
      lags = lags,
      LjungBox_p = lb_p
    )
  )
}

topmodelfiles <- c(
  "01_raw_SARIMA_013_110_diagnostics.png",
  "02_raw_SARIMA_211_110_diagnostics.png",
  "03_guerrero_SARIMA_013_110_diagnostics.png",
  "04_guerrero_SARIMA_211_110_diagnostics.png",
  "05_guerrero_SARIMA_211_010_diagnostics.png"
)

topmodeldiagnostics <- vector("list", length(topmodels))
names(topmodeldiagnostics) <- topmodelnames

for (i in seq_along(topmodels)) {

  png(
    file.path(topmodelfiguredir, topmodelfiles[i]),
    width = 1800,
    height = 1400,
    res = 180,
    bg = "white"
  )

  topmodeldiagnostics[[i]] <- topmodeldiagnosticplot(
    topmodels[[i]],
    topmodelnames[i],
    D = topmodelseasonal[[i]]$order[2],
    fitdf = topmodelfitdf[i]
  )

  dev.off()
}

# Save results

figureguide <- c(
  "Report figures",
  "",
  "Exploratory plots",
  "  Figures/01_original_series.png",
  "  Figures/02_annual_seasonal_overlay.png",
  "  Figures/03_original_acf_pacf.png",
  "  Figures/04_boxcox_variance_check.png",
  "  Figures/05_candidate_differences.png",
  "  Figures/06_ordinary_difference_acf_pacf.png",
  "  Figures/07_seasonal_difference_acf_pacf.png",
  "  Figures/08_double_difference_acf_pacf.png",
  "",
  "Top model diagnostics",
  "  Figures/Top_Models/01_raw_SARIMA_013_110_diagnostics.png",
  "  Figures/Top_Models/02_raw_SARIMA_211_110_diagnostics.png",
  "  Figures/Top_Models/03_guerrero_SARIMA_013_110_diagnostics.png",
  "  Figures/Top_Models/04_guerrero_SARIMA_211_110_diagnostics.png",
  "  Figures/Top_Models/05_guerrero_SARIMA_211_010_diagnostics.png",
  "",
  "Joint Box-Cox likelihood figure",
  "  Figures/09_joint_boxcox_profile_likelihood.png",
  "",
  "Tables",
  "  Tables/boxcox_lambda_values.csv",
  "  Tables/joint_boxcox_estimates.csv",
  "  Tables/joint_boxcox_profile_likelihoods.csv",
  "  Tables/stationarity_tests.csv",
  "  Tables/candidate_model_comparison.csv",
  "  Tables/top_model_comparison.csv",
  "  Tables/top_model_parameters.csv"
)

writeLines(
  figureguide,
  file.path(rundir, "T200_report_figure_guide.txt")
)

sink(file.path(rundir, "T200_analysis_summary.txt"))

cat("MATH5845 T200 analysis summary\n\n")
cat("Observations:", length(y), "\n")
cat("Seasonal period:", seasonalperiod, "weeks\n")
cat("ndiffs(y):", d_rec, "\n")
cat("nsdiffs(y):", D_rec, "\n\n")

cat("Box-Cox values\n")
print(lambda_table, row.names = FALSE)
cat("\nGuerrero lambda used for top models:", lambda_guerrero, "\n")

cat("\njoint Joint Box-Cox likelihood\n")
print(joint_boxcox_estimates, digits = 5, row.names = FALSE)

cat("\ntop models\n")
print(topmodeltable, digits = 4, row.names = FALSE)

cat("\nauto.arima - raw series\n")
print(auto_raw)

cat("\nauto.arima - lecturer Box-Cox value\n")
print(auto_boxcox)

sink()

cat("\nanalysis complete\n")
cat("Run folder:", rundir, "\n")
cat("Figures:", figuredir, "\n")
cat("Tables:", tabledir, "\n")
cat("Model output:", textdir, "\n")

# Frequency-domain analysis
# Uses the same T200 series and the same run folder created above.

frequencyfiguredir <- file.path(figuredir, "Frequency_Domain")
dir.create(frequencyfiguredir, recursive = TRUE, showWarnings = FALSE)

# Raw periodogram
# spec.pgram frequencies are already in cycles per year because y has frequency 52.
raw_spec <- spec.pgram(
  y,
  spans = NULL,
  taper = 0,
  demean = TRUE,
  detrend = FALSE,
  log = "no",
  plot = FALSE
)

raw_freq <- raw_spec$freq
raw_power <- raw_spec$spec

raw_local_peak_id <- which(
  raw_power[2:(length(raw_power) - 1L)] > raw_power[1:(length(raw_power) - 2L)] &
    raw_power[2:(length(raw_power) - 1L)] > raw_power[3:length(raw_power)]
) + 1L

raw_periodogram_peaks <- data.frame(
  Frequency_cycles_per_year = raw_freq[raw_local_peak_id],
  Spectral_power = raw_power[raw_local_peak_id]
)
raw_periodogram_peaks$Period_years <- 1 / raw_periodogram_peaks$Frequency_cycles_per_year
raw_periodogram_peaks$Period_weeks <- seasonalperiod / raw_periodogram_peaks$Frequency_cycles_per_year
raw_periodogram_peaks <- raw_periodogram_peaks[
  order(raw_periodogram_peaks$Spectral_power, decreasing = TRUE),
  , drop = FALSE
]

write.csv(
  raw_periodogram_peaks,
  file.path(tabledir, "frequency_domain_raw_periodogram_peaks.csv"),
  row.names = FALSE
)

raw_top_peaks <- head(raw_periodogram_peaks, 5L)

dominant_frequency <- raw_periodogram_peaks$Frequency_cycles_per_year[1L]
dominant_period_weeks <- raw_periodogram_peaks$Period_weeks[1L]
second_frequency <- raw_periodogram_peaks$Frequency_cycles_per_year[2L]
second_period_weeks <- raw_periodogram_peaks$Period_weeks[2L]
harmonic_ratio <- second_frequency / dominant_frequency

png(
  file.path(frequencyfiguredir, "01_raw_periodogram.png"),
  width = 1800,
  height = 1000,
  res = 180,
  bg = "white"
)
plot(
  raw_freq,
  raw_power,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral power",
  main = "T200: Raw Periodogram"
)
points(
  raw_top_peaks$Frequency_cycles_per_year,
  raw_top_peaks$Spectral_power,
  pch = 19
)
text(
  raw_top_peaks$Frequency_cycles_per_year,
  raw_top_peaks$Spectral_power,
  labels = round(raw_top_peaks$Frequency_cycles_per_year, 3),
  pos = 3,
  cex = 0.8
)
abline(v = c(1, 2), lty = 2)
dev.off()


# Remove a linear trend before recalculating the periodogram.
time_index <- as.numeric(time(y))
trend_model <- lm(as.numeric(y) ~ time_index)

y_detrended <- ts(
  residuals(trend_model),
  start = start(y),
  frequency = frequency(y)
)

png(
  file.path(frequencyfiguredir, "02_linear_trend_and_detrended_series.png"),
  width = 1800,
  height = 1450,
  res = 180,
  bg = "white"
)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(
  y,
  type = "l",
  lwd = 1.3,
  xlab = "Year",
  ylab = "Electricity consumption",
  main = "T200: Series with Fitted Linear Trend"
)
lines(
  time_index,
  as.numeric(fitted(trend_model)),
  lwd = 2
)
plot(
  y_detrended,
  type = "l",
  lwd = 1.3,
  xlab = "Year",
  ylab = "Detrended electricity consumption",
  main = "T200: Detrended Series"
)
abline(h = 0, lty = 2)
par(mfrow = c(1, 1))
dev.off()


detrended_spec <- spec.pgram(
  y_detrended,
  spans = NULL,
  taper = 0,
  demean = TRUE,
  detrend = FALSE,
  log = "no",
  plot = FALSE
)

png(
  file.path(frequencyfiguredir, "03_detrended_periodogram.png"),
  width = 1800,
  height = 1000,
  res = 180,
  bg = "white"
)
plot(
  detrended_spec$freq,
  detrended_spec$spec,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral power",
  main = "T200: Detrended Periodogram"
)
abline(v = c(1, 2), lty = 2)
dev.off()


# Daniell-smoothed periodogram of the detrended series.
smoothed_spec <- spec.pgram(
  y_detrended,
  spans = c(3, 3),
  taper = 0.1,
  demean = TRUE,
  detrend = FALSE,
  log = "no",
  plot = FALSE
)

smooth_freq <- smoothed_spec$freq
smooth_power <- smoothed_spec$spec
smooth_local_peak_id <- which(
  smooth_power[2:(length(smooth_power) - 1L)] > smooth_power[1:(length(smooth_power) - 2L)] &
    smooth_power[2:(length(smooth_power) - 1L)] > smooth_power[3:length(smooth_power)]
) + 1L

smoothed_periodogram_peaks <- data.frame(
  Frequency_cycles_per_year = smooth_freq[smooth_local_peak_id],
  Spectral_density = smooth_power[smooth_local_peak_id]
)
smoothed_periodogram_peaks$Period_years <- 1 / smoothed_periodogram_peaks$Frequency_cycles_per_year
smoothed_periodogram_peaks$Period_weeks <- seasonalperiod / smoothed_periodogram_peaks$Frequency_cycles_per_year
smoothed_periodogram_peaks <- smoothed_periodogram_peaks[
  order(smoothed_periodogram_peaks$Spectral_density, decreasing = TRUE),
  , drop = FALSE
]

write.csv(
  smoothed_periodogram_peaks,
  file.path(tabledir, "frequency_domain_smoothed_periodogram_peaks.csv"),
  row.names = FALSE
)

png(
  file.path(frequencyfiguredir, "04_smoothed_periodogram.png"),
  width = 1800,
  height = 1000,
  res = 180,
  bg = "white"
)
plot(
  smooth_freq,
  smooth_power,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral density",
  main = "T200: Smoothed Periodogram"
)
abline(v = c(1, 2), lty = 2)
dev.off()


# Parametric spectral estimate.
# Fit AR(p) to the detrended series by Yule-Walker and let AIC choose p,
# then compare its spectral density with the smoothed periodogram.
ar_spectral_fit <- ar(
  y_detrended,
  method = "yw",
  aic = TRUE,
  order.max = 15
)

ar_spectrum <- spec.ar(
  ar_spectral_fit,
  n.freq = 500,
  log = "no",
  plot = FALSE
)

png(
  file.path(frequencyfiguredir, "05_smoothed_periodogram_vs_ar_spectrum.png"),
  width = 1800,
  height = 1000,
  res = 180,
  bg = "white"
)
plot(
  smooth_freq,
  smooth_power,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral density",
  main = paste0(
    "T200: Smoothed Periodogram vs AR(",
    ar_spectral_fit$order,
    ") Spectral Density"
  )
)
lines(
  ar_spectrum$freq,
  ar_spectrum$spec,
  lty = 2,
  lwd = 2
)
legend(
  "topright",
  legend = c(
    "Smoothed periodogram",
    paste0("AR(", ar_spectral_fit$order, ") spectral density")
  ),
  lty = c(1, 2),
  lwd = c(1.3, 2),
  bty = "n"
)
dev.off()


# One report-ready comparison figure.
png(
  file.path(frequencyfiguredir, "06_frequency_domain_summary.png"),
  width = 1800,
  height = 1800,
  res = 180,
  bg = "white"
)
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
plot(
  raw_freq,
  raw_power,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral power",
  main = "Raw Periodogram"
)
abline(v = c(1, 2), lty = 2)
plot(
  detrended_spec$freq,
  detrended_spec$spec,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral power",
  main = "Detrended Periodogram"
)
abline(v = c(1, 2), lty = 2)
plot(
  smooth_freq,
  smooth_power,
  type = "l",
  lwd = 1.3,
  xlab = "Frequency (cycles per year)",
  ylab = "Spectral density",
  main = "Smoothed Periodogram"
)
abline(v = c(1, 2), lty = 2)
par(mfrow = c(1, 1))
dev.off()


frequency_domain_summary <- data.frame(
  Quantity = c(
    "Dominant frequency (cycles per year)",
    "Dominant period (weeks)",
    "Second frequency (cycles per year)",
    "Second period (weeks)",
    "Second-to-first frequency ratio",
    "AR order selected by AIC"
  ),
  Value = c(
    dominant_frequency,
    dominant_period_weeks,
    second_frequency,
    second_period_weeks,
    harmonic_ratio,
    ar_spectral_fit$order
  )
)

write.csv(
  frequency_domain_summary,
  file.path(tabledir, "frequency_domain_summary.csv"),
  row.names = FALSE
)

cat("\nFrequency-domain analysis complete.\n")
cat("Frequency-domain figures saved to:", frequencyfiguredir, "\n")
print(frequency_domain_summary)
