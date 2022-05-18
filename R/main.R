#' naive
#'
#' @param df A data frame with time features on columns. In case of missing values, automatic missing imputation through kalman filter will be performed.
#' @param seq_len Positive integer. Time-step number of the forecasting sequence. Default: NULL (random selection within boundaries).
#' @param ci Confidence interval for prediction. Default: 0.8
#' @param smoother Logical. Flag to TRUE for loess smoothing. Default: FALSE.
#' @param cover Positive numeric. The quantile cover around the location parameter (between 0 and 1). Default: NULL (random selection within boundaries).
#' @param stride Positive integer. Shift between subsequent sequences. Default: NULL (random selection within boundaries).
#' @param method String. Distance method using during the comparison of time sequences. Possible options are: "euclidean", "manhattan", "canberra1", "minimum", "maximum", "minkowski", "bhattacharyya", "kullback_leibler", "jensen_shannon". Default: NULL (random selection within boundaries).
#' @param location String. Statistic used to center the cover parameter. Possible options are: "mean", "mode" (parzen method), "median". Default: NULL (random selection within boundaries).
#' @param n_windows Positive integer. Number of validation windows to test prediction error. Default: 10.
#' @param n_samp Positive integer. Number of sample selected during random search. Default: 30.
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: collection of all the models explored with random search
#' \item history: a table with the explored models' hyper-parameters and validation errors
#' \item best_model: best combination resulting from the average prediction score across different ranks and features, including:
#' \itemize{
#' \item quant_preds: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, median range ratio, upside probability and divergence for each point fo predicted sequences
#' \item errors: testing errors for each time feature averaged across validation windows
#' \item plots: standard plot with confidence interval for each time feature
#' }
#' \item time_log
#' }
#'
#' @export
#'
#' @import purrr
#' @import tictoc
#' @import greybox
#' @import ggplot2
#' @importFrom scales number
#' @importFrom readr parse_number
#' @importFrom lubridate seconds_to_period is.Date as.duration
#' @importFrom modeest mlv1
#' @importFrom moments kurtosis skewness
#' @import stats
#' @importFrom imputeTS na_kalman
#' @importFrom fANCOVA loess.as
#' @importFrom utils combn head tail
#' @importFrom Rfast Dist

#'@examples
#'naive(time_features[,c(2, 3)], seq_len = 100, n_samp = 1, n_windows = 3)
#'


###
naive <- function(df, seq_len = NULL, ci = 0.8, smoother = FALSE, cover = NULL, stride = NULL, method = NULL, location = NULL, n_windows = 10, n_samp = 30, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  tic.clearlog()
  tic("time")

  set.seed(seed)
  n_length <- nrow(df)
  n_feats <- ncol(df)
  feat_names <- colnames(df)

  if(anyNA(df)){df <- as.data.frame(map(df, ~ na_kalman(.x))); message("kalman imputation on time features\n")}
  if(smoother==TRUE){df <- as.data.frame(map(df, ~ suppressWarnings(loess.as(x=1:n_length, y=.x)$fitted))); message("performing optimal smoothing\n")}
  if(any(seq_len <= 0) | any(stride <= 0) | any(cover >= 1) | any(cover <= 0)){stop("at least one parameter out of boundaries")}
  if(length(seq_len) == 1 & length(cover) == 1 & length(stride) == 1 & length(method) == 1 & length(location) == 1){n_samp <- 1}

  deriv <- map_dbl(df, ~ best_deriv(.x))
  max_stride <- function(n_length, seq_len, n_row = 5) {floor(((n_length - n_length%%seq_len - max(deriv)) - seq_len)/(n_row - 1))}

  sqln_set <- sampler(seq_len, n_samp, range = c(1, round(sqrt(n_length))), integer = TRUE)
  cvr_set <- sampler(cover, n_samp, range = c(0.1, 0.9), integer = FALSE)
  max_set <- map_dbl(sqln_set, ~ floor(max_stride(n_length/(n_windows + 1), .x)))
  fix_set <- map2_dbl(max_set, sqln_set, ~ sample(min(c(.x, .y)), 1))
  if(any(fix_set < 1) & is.null(stride)){stop("not enough data for the validation windows")}
  strd_set <- sampler(stride, n_samp, range = NULL, integer = TRUE, fun = fix_set)
  mthd_set <- sampler(method, n_samp, range = c("euclidean", "manhattan", "canberra1", "minimum", "maximum", "minkowski", "bhattacharyya", "kullback_leibler", "jensen_shannon"), integer = FALSE)
  lctn_set <- sampler(location, n_samp, range = c("mean", "median", "mode"), integer = FALSE)

  exploration <- pmap(list(sqln_set, cvr_set, strd_set, mthd_set, lctn_set), ~ tryCatch(windower(df, ..1, ci, ..2, ..3, ..4, ..5, n_windows, error_scale, error_benchmark, dates), error = function(e) NA))

  collected <- flatten(map(exploration, ~ .x["errors"]))
  avg_error <- as.data.frame(Reduce(rbind, map(collected, ~ apply(.x, 2, mean))))
  if(n_samp == 1){avg_error <- t(as.data.frame(avg_error))}
  history <- data.frame(seq_len = sqln_set, cover = round(cvr_set, 4), stride = strd_set, method = mthd_set, location = lctn_set, round(avg_error, 4))
  rownames(history) <- NULL
  if(n_samp > 1){
    ranking_order <- order(- avg_error$pred_score, abs(avg_error$me), avg_error$mae,  avg_error$mse, avg_error$rmsse, abs(avg_error$mpe), avg_error$mape, avg_error$rmae, avg_error$rrmse, avg_error$rame, avg_error$mase, avg_error$smse, abs(avg_error$sce), avg_error$gmrae)
    history <- history[ranking_order,]
    best_idx <- as.integer(rownames(history)[1])
    best_model <- exploration[[best_idx]]}
  else {best_model <- exploration[[1]]}

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)

  return(outcome)
}

###
windower <- function(df, seq_len, ci = 0.8, cover = 0.5, stride = 1, method = "euclidean", location = "mode", n_windows = 10, error_scale, error_benchmark, dates = NULL)
{
  feat_names <- colnames(df)
  n_length <- nrow(df)
  n_feats <- ncol(df)

  idx <- c(rep(1, n_length%%(n_windows + 1)), rep(1:(n_windows + 1), each = n_length/(n_windows + 1)))

  window_results <- map(1:(n_windows + 1), ~ engine(df[idx <= .x,, drop = FALSE], seq_len, ci, cover, stride, method, location, dates))
  ground_truth <- map(1:n_windows, ~ head(df[idx == (.x + 1),, drop = FALSE], seq_len))
  window_hist <- map(1:n_windows, ~ df[idx <= .x,, drop = FALSE])
  window_quant_preds <- map_depth(map(window_results, ~ .x$quant_preds), 2, ~ .x[,"mean"])
  window_raw_preds <- map(window_results, ~ .x$raw_preds)

  pred_score <- map_dbl(transpose(map2(head(window_raw_preds, -1), ground_truth, ~ map2(.x, .y, ~ prediction_score(.x, .y)))), ~ mean(unlist(.x)))

  reducer <- function(vector_list)
  {
    if(length(vector_list) > 1){out <- colMeans(Reduce(rbind, vector_list))}
    if(length(vector_list) == 1){out <- as.matrix(unlist(vector_list))}
    return(out)
  }

  errors <- as.data.frame(round(cbind(pred_score = pred_score, t(as.data.frame(map(transpose(pmap(list(ground_truth, head(window_quant_preds, -1), window_hist),  ~ pmap(list(..1, ..2, ..3), ~ my_metrics(..1, ..2, ..3, error_scale, error_benchmark)))), ~ reducer(.x))))), 4))
  final <- flatten(tail(window_results, 1))
  model <- list(quant_preds = final$quant_preds, plots = final$plots, errors = errors)

  return(model)
}

###
engine <- function(df, seq_len, ci = 0.8, cover = 0.5, stride = 1, method = "euclidean", location = "mode", dates = NULL)
{
  feat_names <- colnames(df)
  n_length <- nrow(df)
  deriv <- map_dbl(df, ~ best_deriv(.x))

  diff_models <- map2(df, deriv, ~ recursive_diff(.x, .y))
  ddf <- as.data.frame(map(diff_models, ~ smart_tail(.x$vector, n_length - max(deriv))))
  segmented_sets <- map(ddf, ~ smart_reframer(.x, seq_len, stride))
  dmat_sets <- map(segmented_sets, ~ as.matrix(Dist(.x, method, p = 3)))

  selector <- function(mat, cover, location)
  {
    diag_out <- sort(mat[mat != 0])
    min_value <- min(diag_out)
    max_value <- max(diag_out)
    if(location == "mode"){loc <- suppressWarnings(mlv1(diag_out, method = "parzen"))}
    if(location == "mean"){loc <- suppressWarnings(mean(diag_out))}
    if(location == "median"){loc <- suppressWarnings(median(diag_out))}
    min_fth <- tail(diag_out[diag_out <= loc], 1)###min_feasible_threshold
    max_fth <- head(diag_out[diag_out >= loc], 1)###max_feasible_threshold
    mode_quant <- ecdf(diag_out)(loc)
    min_cover <- mode_quant - cover/2
    max_cover <- mode_quant + cover/2
    min_cover <- min_cover * (min_cover >= 0)
    max_cover <- max_cover * (max_cover <= 1) + 1 * (max_cover > 1)
    th <- quantile(diag_out, probs = c(min_cover, max_cover))
    mask <- (mat != 0) & (mat <= max(th)) & (mat >= min(th))
    if(all(mask == FALSE)){mask <- (mat != 0) & (mat <= max_fth) & (mat >= min_fth)}
    idx <- unique(as.vector(which(mask == TRUE, arr.ind = TRUE)))

    return(idx)
  }

  idx_sets <- map(dmat_sets, ~ selector(.x, cover, location))
  core_sets <- map2(segmented_sets, idx_sets,  ~ as.data.frame(.x[.y, ]))

  empirical <- function(vect, n){sample(vect, size = n, replace = TRUE)}

  if(seq_len > 1){
    diff_preds <- map(map_depth(core_sets, 2, ~ empirical(.x, 1000)), ~ Reduce(cbind, .x))
    raw_preds <- map2(diff_models, diff_preds, ~ t(apply(.y, 1, function(x) invdiff(x, .x$tail_value))))
  }

  if(seq_len == 1){
    diff_preds <- map(map_depth(core_sets, 2, ~ empirical(.x, 1000)), ~ t(Reduce(cbind, .x)))
    raw_preds <- map2(diff_models, diff_preds, ~ apply(.y, 1, function(x) invdiff(x, .x$tail_value)))
  }

  raw_preds <-  map2(df, raw_preds, ~ doxa_filter(.x, .y))

  quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))
  p_stats <- function(x){stats <- c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = tryCatch(suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA), skewness = tryCatch(suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA)); return(stats)}
  quant_preds <- map(raw_preds, ~ as.data.frame(t(as.data.frame(apply(.x, 2, p_stats)))))

  iqr_to_range <- map(quant_preds, ~ tryCatch((.x[, "75%"] - .x[, "25%"])/(.x[, "max"] - .x[, "min"]), error = function(e) NA))
  median_range_ratio <- map(quant_preds, ~ tryCatch((.x[, "max"] - .x[, "50%"])/(.x[, "50%"] - .x[, "min"]), error = function(e) NA))
  growths <- map(quant_preds, ~ mapply(function(m, s) rnorm(1000, m, s), m = .x[, "mean"], s = .x[, "sd"]))
  upside_prob <- map(growths, ~ tryCatch(c(NA, colMeans(apply(.x[,-1, drop = FALSE]/.x[,-ncol(.x), drop = FALSE], 2, function(x) x > 1))), error = function(e) NA))
  pvalues <- map(quant_preds, ~ mapply(function(m, s) pnorm(seq(min(.x[, "min"]), max(.x[, "max"]), length.out = 1000), m, s), m = .x[, "mean"], s = .x[, "sd"]))
  divergence <- map(pvalues, ~ tryCatch(c(NA, apply(.x[,-1, drop = FALSE] - .x[,-ncol(.x), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE)))), error = function(e) NA))
  quant_preds <- pmap(list(quant_preds, iqr_to_range, median_range_ratio, upside_prob, divergence), ~ round(cbind(..1, iqr_to_range = ..2, median_range_ratio = ..3, upside_prob = ..4, divergence = ..5), 4))
  names(quant_preds) <- feat_names

  if(is.null(dates)){hist_dates <- 1:n_length; forcat_dates <- (n_length + 1):(n_length + seq_len); quant_preds <- map(quant_preds, ~ {rownames(.x) <- paste0("t", 1:seq_len); return(.x)})}
  if(!is.null(dates)){hist_dates <- tail(dates, n_length); forcat_dates <- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len); quant_preds <- map(quant_preds, ~ {rownames(.x) <- as.character(forcat_dates); return(.x)})}
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_labs <- paste0("Forecasting Values for ", feat_names)

  plots <- pmap(list(df, quant_preds, y_labs), ~ ts_graph(x_hist = hist_dates, y_hist = ..1, x_forcat = forcat_dates, y_forcat = ..2[,"50%"], lower = ..2[,2], upper = ..2[,6], label_y = ..3))

  outcome <- list(quant_preds = quant_preds, raw_preds = raw_preds, plots = plots)
  return(outcome)
}

###
smart_reframer <- function(ts, seq_len, stride)
{
  n_length <- length(ts)
  if(seq_len > n_length | stride > n_length){stop("vector too short for sequence length or stride")}
  if(n_length%%seq_len > 0){ts <- tail(ts, - (n_length%%seq_len))}
  n_length <- length(ts)
  idx <- base::seq(from = 1, to = (n_length - seq_len + 1), by = 1)
  reframed <- t(sapply(idx, function(x) ts[x:(x+seq_len-1)]))
  if(seq_len == 1){reframed <- t(reframed)}
  idx <- rev(base::seq(nrow(reframed), 1, - stride))
  reframed <- reframed[idx,,drop = FALSE]
  colnames(reframed) <- paste0("t", 1:seq_len)
  return(reframed)
}

###
best_deriv <- function(ts, max_diff = 3, thresh = 0.001)
{
  pvalues <- vector(mode = "double", length = as.integer(max_diff))

  for(d in 1:(max_diff + 1))
  {
    model <- lm(ts ~ t, data.frame(ts, t = 1:length(ts)))
    pvalues[d] <- with(summary(model), pf(fstatistic[1], fstatistic[2], fstatistic[3],lower.tail=FALSE))
    ts <- diff(ts)
  }

  best <- tail(cumsum(pvalues < thresh), 1)

  return(best)
}

###
recursive_diff <- function(vector, deriv)
{
  vector <- unlist(vector)
  head_value <- vector("numeric", deriv)
  tail_value <- vector("numeric", deriv)
  if(deriv==0){head_value = NULL; tail_value = NULL}
  if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
  outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
  return(outcome)
}

###
invdiff <- function(vector, heads, add = FALSE)
{
  vector <- unlist(vector)
  if(is.null(heads)){return(vector)}
  for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
  if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
}

###
smart_head <- function(x, n)
{
  if(n != 0){return(head(x, n))}
  if(n == 0){return(x)}
}

###
smart_tail <- function(x, n)
{
  if(n != 0){return(tail(x, n))}
  if(n == 0){return(x)}
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "seagreen2", forcat_line = "seagreen4", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = c(y_hist, y_forcat))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- lower; forcat_data$upper <- upper}

  plot <- ggplot()+ geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}

###
my_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive")
{
  scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
  benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))

  me <- ME(holdout, forecast, na.rm = TRUE)
  mae <- MAE(holdout, forecast, na.rm = TRUE)
  mse <- MSE(holdout, forecast, na.rm = TRUE)
  rmsse <- RMSSE(holdout, forecast, scale, na.rm = TRUE)
  mre <- MRE(holdout, forecast, na.rm = TRUE)
  mpe <- MPE(holdout, forecast, na.rm = TRUE)
  mape <- MAPE(holdout, forecast, na.rm = TRUE)
  rmae <- rMAE(holdout, forecast, benchmark, na.rm = TRUE)
  rrmse <- rRMSE(holdout, forecast, benchmark, na.rm = TRUE)
  rame <- rAME(holdout, forecast, benchmark, na.rm = TRUE)
  mase <- MASE(holdout, forecast, scale, na.rm = TRUE)
  smse <- sMSE(holdout, forecast, scale, na.rm = TRUE)
  sce <- sCE(holdout, forecast, scale, na.rm = TRUE)
  gmrae <- GMRAE(holdout, forecast, benchmark, na.rm = TRUE)

  out <- round(c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae), 3)
  return(out)
}

###
doxa_filter <- function(orig, mat, n_class = NULL)
{
  orig <- dim(mat)
  discrete_check <- all(orig%%1 == 0)
  all_positive_check <- all(orig >= 0)
  all_negative_check <- all(orig <= 0)
  monotonic_increase_check <- all(diff(orig) >= 0)
  monotonic_decrease_check <- all(diff(orig) <= 0)
  class_check <- FALSE
  if(is.integer(n_class)){class_check <- length(unique(orig)) <= n_class}

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(discrete_check){mat <- floor(mat)}
  if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
  if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
  if(class_check){mat[!(mat %in% unique(orig))] <- ((mat[!(mat %in% unique(orig))] > max(unique(orig))) * max(unique(orig))) + ((mat[!(mat %in% unique(orig))] < min(unique(orig))) * min(unique(orig)))}
  if(all_positive_check){mat[mat < 0] <- 0}
  if(all_negative_check){mat[mat > 0] <- 0}

  dim(mat) <- orig
  return(mat)
}

###
sampler <- function(vect, n_samp, range = NULL, integer = FALSE, fun = NULL)
{
  if(is.null(vect) & is.null(fun))
  {
    if(!is.character(range)){if(integer){set <- min(range):max(range)} else {set <- seq(min(range), max(range), length.out = 1000)}} else {set <- range}
    samp <- sample(set, n_samp, replace = TRUE)
  }

  if(is.null(vect) & !is.null(fun)){samp <- fun}

  if(length(vect)==1){samp <- rep(vect, n_samp)}
  if(length(vect) > 1){samp <- sample(vect, n_samp, replace = TRUE)}
  return(samp)
}

###
prediction_score <- function(integrated_preds, ground_truth)
{
  pfuns <- apply(integrated_preds, 2, ecdf)
  pvalues <- map2_dbl(pfuns, ground_truth, ~ .x(.y))
  scores <- mean(1 - 2 * abs(pvalues - 0.5))
  return(scores)
}
