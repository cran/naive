#' naive
#'
#' @param df A data frame with time features on columns (all numerics or all categories, but not both). In case of missing values, automatic missing imputation through kalman filter will be performed.
#' @param seq_len Positive integer. Time-step number of the forecasting sequence. Default: NULL (random selection within boundaries).
#' @param ci Confidence interval for prediction. Default: 0.8
#' @param smoother Logical. Flag to TRUE for loess smoothing (only for numeric series). Default: FALSE.
#' @param cover Positive numeric. The quantile cover around the location parameter (between 0 and 1). Default: NULL (random selection within boundaries).
#' @param stride Positive integer. Shift between subsequent sequences. Default: NULL (random selection within boundaries).
#' @param method String. Distance method using during the comparison of time sequences. Possible options are: "euclidean", "manhattan", "minkowski". Default: NULL (random selection).
#' @param location String. Statistic used to center the cover parameter. Possible options are: "mean", "mode" (parzen method), "median". Default: NULL (random selection).
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
#' \item quant_preds: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, above to below median range, upside probability and divergence for each point fo predicted sequences
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
#' @import fastDummies
#' @import entropy
#' @import philentropy

#'@examples
#'\dontrun{
#'naive(time_features[, 2:3, drop = F], seq_len = 30, n_samp = 1, n_windows = 5)
#'naive(time_features[, 3, drop = F], seq_len = 20, n_samp = 10, n_windows = 5)
#'naive(time_features[, -1, drop = F], seq_len = 60, n_samp = 15, n_windows = 5)
#'}


###
naive <- function(df, seq_len = NULL, ci = 0.8, smoother = FALSE, cover = NULL, stride = NULL, method = NULL, location = NULL, n_windows = 10, n_samp = 30, dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  tic.clearlog()
  tic("time")

  set.seed(seed)

  if(!is.data.frame(df)){stop("time features must be in dataframe format")}
  if(n_windows < 2){n_windows <- 2; message("setting validation to the minimum number of windows (2)")}

  n_length <- nrow(df)

  class_index <- any(map_lgl(df, ~ is.factor(.x) | is.character(.x)))
  all_classes <- all(class_index)
  numeric_index <- map_lgl(df, ~ is.integer(.x) | is.numeric(.x))
  all_numerics <- all(numeric_index)
  if(!(all_classes | all_numerics)){stop("only all numerics or all classes, not both")}

  if(all_classes){df <- dummy_cols(df, select_columns = NULL, remove_first_dummy = FALSE, remove_most_frequent_dummy = TRUE, ignore_na = FALSE, split = NULL, remove_selected_columns = TRUE); binary_class <- rep(T, ncol(df))}
  if(all_numerics){binary_class <- rep(F, ncol(df))}

  if(anyNA(df) & all_numerics){df <- as.data.frame(na_kalman(df)); message("kalman imputation on time features\n")}
  if(anyNA(df) & all_classes){df <- floor(as.data.frame(na_kalman(df))); message("kalman imputation on time features\n")}
  if(smoother == TRUE & all_numerics){df <- as.data.frame(purrr::map(df, ~ suppressWarnings(loess.as(x=1:n_length, y=.x)$fitted))); message("performing optimal smoothing\n")}


  n_feats <- ncol(df)
  feat_names <- colnames(df)

  if(any(seq_len <= 0) | any(stride <= 0) | any(cover >= 1) | any(cover <= 0)){stop("at least one parameter out of boundaries")}
  if(length(seq_len) == 1 & length(cover) == 1 & length(stride) == 1 & length(method) == 1 & length(location) == 1){n_samp <- 1}

  deriv <- map_dbl(df, ~ best_deriv(.x))
  max_stride <- function(n_length, seq_len, n_row = 5) {floor(((n_length - n_length%%seq_len - max(deriv)) - seq_len)/(n_row - 1))}

  min_limit <- 2 + max(deriv)
  max_limit <- round(n_length/(2*(n_windows + 1)))
  if(max_limit < min_limit){stop("not enough data for validation windows")}

  if(is.numeric(seq_len) && any(seq_len < min_limit)){seq_len[seq_len < min_limit] <- min_limit; message("fixing seq_len for min limit")}
  if(is.numeric(seq_len) && any(seq_len > max_limit)){seq_len[seq_len > max_limit] <- max_limit; message("fixing seq_len for max limit")}
  sqln_set <- sampler(seq_len, n_samp, range = c(min_limit, max_limit), integer = TRUE)

  cvr_set <- sampler(cover, n_samp, range = c(0.1, 0.9), integer = FALSE)
  strd_set <- sampler(stride, n_samp, range = NULL, integer = TRUE, fun = map_dbl(sqln_set, ~ sample(ceiling(sqrt(.x)), 1)))###STRIDE MUST DEPENDS ON SEQ_LEN, NEED TO FIX
  mthd_set <- sampler(method, n_samp, range = c("euclidean", "manhattan", "minkowski"), integer = FALSE)
  lctn_set <- sampler(location, n_samp, range = c("mean", "median", "mode"), integer = FALSE)

  hyper_params <- list(sqln_set, cvr_set, strd_set, mthd_set, lctn_set)

  exploration <- pmap(hyper_params, ~ tryCatch(windower(df, ..1, ci, ..2, ..3, ..4, ..5, n_windows, error_scale, error_benchmark, binary_class, dates, seed), error = function(e) NA))

  collected <- flatten(map(exploration, ~ .x["errors"]))
  collected_dim <- length(dim(collected[[1]]))
  if(collected_dim == 2){avg_error <- as.data.frame(Reduce(rbind, map(collected, ~ apply(.x, 2, mean))))}
  if(collected_dim == 1){avg_error <- as.data.frame(Reduce(rbind, collected))}
  if(n_samp == 1){avg_error <- t(as.data.frame(avg_error))}
  history <- data.frame(seq_len = sqln_set, cover = round(cvr_set, 4), stride = strd_set, method = mthd_set, location = lctn_set, round(avg_error, 4))
  rownames(history) <- NULL
  if(n_samp > 1){
    if(all_numerics == TRUE){history <- ranker(history, focus = -c(1:5), inverse = NULL, absolute = c("me", "mpe", "sce"), reverse = FALSE)}
    if(all_classes == TRUE){history <- ranker(history, focus = -c(1:5), inverse = NULL, absolute = NULL, reverse = FALSE)}
    best_idx <- as.integer(rownames(history)[1])
    best_model <- exploration[[best_idx]]
    }

  else {best_model <- exploration[[1]]}

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)

  return(outcome)
}

###
windower <- function(df, seq_len, ci = 0.8, cover = 0.5, stride = 1, method = "euclidean", location = "mode", n_windows = 10, error_scale, error_benchmark, binary_class, dates = NULL, seed)
{
  feat_names <- colnames(df)
  n_length <- nrow(df)
  n_feats <- ncol(df)

  idx <- c(rep(1, n_length%%(n_windows + 1)), rep(1:(n_windows + 1), each = n_length/(n_windows + 1)))

  window_results <- map(1:(n_windows + 1), ~ engine(df[idx <= .x,, drop = FALSE], seq_len, ci, cover, stride, method, location, binary_class, dates[idx <= .x], seed))
  ground_truth <- map(1:n_windows, ~ head(df[idx == (.x + 1),, drop = FALSE], seq_len))
  window_hist <- map(1:n_windows, ~ df[idx <= .x,, drop = FALSE])
  if(!all(binary_class)){window_quant_preds <- map_depth(map(window_results, ~ .x$quant_preds), 2, ~ .x[,"mean"])}
  if(all(binary_class)){window_quant_preds <- map_depth(map(window_results, ~ .x$quant_preds), 2, ~ .x[,"prop"])}
  window_raw_preds <- map(window_results, ~ .x$raw_preds)


  pred_score <- map(transpose(map2(head(window_raw_preds, -1), ground_truth, ~ map2(.x, .y, ~ prediction_score(.x, .y)))), ~ round(colMeans(Reduce(rbind, .x)), 4))

  reducer <- function(vector_list)
  {
    if(length(vector_list) > 1){out <- colMeans(Reduce(rbind, vector_list))}
    if(length(vector_list) == 1){out <- as.matrix(unlist(vector_list))}
    return(out)
  }

  errors <- as.data.frame(round(cbind(t(as.data.frame(map(transpose(pmap(list(ground_truth, head(window_quant_preds, -1), window_hist),  ~ pmap(list(..1, ..2, ..3, binary_class), ~ custom_metrics(..1, ..2, ..3, error_scale, error_benchmark, binary_class = ..4)))), ~ reducer(.x))))), 4))
  final <- flatten(tail(window_results, 1))
  quant_preds <- map2(final$quant_preds, pred_score, ~ cbind(.x, pred_score = .y))

  plots <- pmap(list(df, quant_preds, feat_names), ~ plotter(quant_pred = ..2, ci, ts = ..1, dates, feat_name = ..3))

  model <- list(quant_preds = quant_preds, plots = plots, errors = errors)

  return(model)
}

###
engine <- function(df, seq_len, ci = 0.8, cover = 0.5, stride = 1, method = "euclidean", location = "mode", binary_class, dates = NULL, seed)
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
    diag_out <- sort(mat[!(mat %in% diag(mat))])
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

  raw_preds <-  pmap(list(df, raw_preds, binary_class), ~ doxa_filter(..1, ..2, ..3))
  quant_preds <-  pmap(list(df, raw_preds, binary_class), ~ fast_qpred(raw_pred = ..2, ts = ..1, ci, error_scale, error_benchmark, binary_class = ..3, dates, seed))

  outcome <- list(quant_preds = quant_preds, raw_preds = raw_preds)
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
  scores <- 1 - 2 * abs(pvalues - 0.5)
  return(scores)
}

###
fast_qpred <- function(raw_pred, ts, ci, error_scale = "naive", error_benchmark = "naive", binary_class = F, dates, seed = 42)
{
  set.seed(seed)

  quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

  if(binary_class == F)
  {
    p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), kurtosis = suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), skewness = suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)))}
    quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
    p_value <- apply(raw_pred, 2, function(x) ecdf(x)(seq(min(raw_pred), max(raw_pred), length.out = 1000)))
    divergence <- c(max(p_value[,1] - seq(0, 1, length.out = 1000)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
    upside_prob <- c(mean((raw_pred[,1]/tail(ts, 1)) > 1, na.rm = T), apply(apply(raw_pred[,-1, drop = FALSE]/raw_pred[,-ncol(raw_pred), drop = FALSE], 2, function(x) x > 1), 2, mean, na.rm = T))
    iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
    above_to_below_range <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
    quant_pred <- round(cbind(quant_pred, iqr_to_range, above_to_below_range, upside_prob, divergence), 4)
  }

  if(binary_class == T)
  {
    p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), prop = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), entropy = entropy(x))}
    quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
    p_value <- apply(raw_pred, 2, function(x) ecdf(x)(c(0, 1)))
    divergence <- c(max(p_value[,1] - c(0, 1)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
    upgrade_prob <- c(mean(((raw_pred[,1] + 1)/tail(ts + 1, 1)) > 1, na.rm = T), apply(apply((raw_pred[,-1, drop = FALSE] + 1)/(raw_pred[,-ncol(raw_pred), drop = FALSE] + 1), 2, function(x) x > 1), 2, mean, na.rm = T))
    quant_pred <- round(cbind(quant_pred, upgrade_prob = upgrade_prob, divergence = divergence), 4)
  }


  if(is.Date(dates))
  {
    new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + nrow(quant_pred) * mean(diff(dates)), length.out = nrow(quant_pred))
    rownames(quant_pred) <- as.character(new_dates)
  }
  else
  {
    rownames(quant_pred) <- paste0("t", 1:nrow(quant_pred))
  }

  return(quant_pred)
}

###
plotter <- function(quant_pred, ci, ts, dates = NULL, feat_name)
{
  seq_len <- nrow(quant_pred)
  n_ts <- length(ts)

  if(is.Date(dates))
  {
    new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len)
    x_hist <- dates
    x_forcat <- new_dates
    rownames(quant_pred) <- as.character(new_dates)
  }
  else
  {
    x_hist <- 1:n_ts
    x_forcat <- (n_ts + 1):(n_ts + seq_len)
    rownames(quant_pred) <- paste0("t", 1:seq_len)
  }

  quant_pred <- as.data.frame(quant_pred)
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_lab <- paste0("Forecasting Values for ", feat_name)

  lower_b <- paste0((1-ci)/2 * 100, "%")
  upper_b <- paste0((ci+(1-ci)/2) * 100, "%")

  plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = quant_pred[, "50%"], lower = quant_pred[, lower_b], upper = quant_pred[, upper_b], label_x = x_lab, label_y = y_lab)
  return(plot)
}

###
doxa_filter <- function(ts, mat, binary_class = F)
{
  discrete_check <- all(ts%%1 == 0)
  all_positive_check <- all(ts >= 0)
  all_negative_check <- all(ts <= 0)
  monotonic_increase_check <- all(diff(ts) >= 0)
  monotonic_decrease_check <- all(diff(ts) <= 0)

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(all_positive_check){mat[mat < 0] <- 0}
  if(all_negative_check){mat[mat > 0] <- 0}
  if(discrete_check){mat <- round(mat)}
  if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
  if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}

  if(binary_class == T){mat[mat > 1] <- 1; mat[mat < 1] <- 0}
  mat <- na.omit(mat)

  return(mat)
}


###
custom_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive", binary_class = F)
{

  if(binary_class == F)
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
  }

  if(binary_class == T)
  {
    dice <- suppressMessages(distance(rbind(holdout, forecast), method = "dice"))
    jaccard <- suppressMessages(distance(rbind(holdout, forecast), method = "jaccard"))
    cosine <- suppressMessages(distance(rbind(holdout, forecast), method = "cosine"))
    canberra <- suppressMessages(distance(rbind(holdout, forecast), method = "canberra"))
    gower <- suppressMessages(distance(rbind(holdout, forecast), method = "gower"))
    tanimoto <- suppressMessages(distance(rbind(holdout, forecast), method = "tanimoto"))
    hassebrook <- 1 - suppressMessages(distance(rbind(holdout, forecast), method = "hassebrook"))
    taneja <- suppressMessages(distance(rbind(holdout, forecast), method = "taneja"))
    lorentzian <- suppressMessages(distance(rbind(holdout, forecast), method = "lorentzian"))
    clark <- suppressMessages(distance(rbind(holdout, forecast), method = "clark"))
    sorensen <- suppressMessages(distance(rbind(holdout, forecast), method = "sorensen"))
    harmonic_mean <- suppressMessages(distance(rbind(holdout, forecast), method = "harmonic_mean"))
    avg <- suppressMessages(distance(rbind(holdout, forecast), method = "avg"))

    out <- round(c(dice, jaccard, cosine, canberra, gower, tanimoto, hassebrook, taneja, lorentzian, clark, sorensen, harmonic_mean, avg), 4)
  }

  return(out)
}

###
ranker <- function(df, focus, inverse = NULL, absolute = NULL, reverse = FALSE)
{
  rank_set <- df[, focus, drop = FALSE]
  if(!is.null(inverse)){rank_set[, inverse] <- - rank_set[, inverse]}###INVERSION BY COL NAMES
  if(!is.null(absolute)){rank_set[, absolute] <- abs(rank_set[, absolute])}###ABS BY COL NAMES
  index <- apply(scale(rank_set), 1, mean, na.rm = TRUE)
  if(reverse == FALSE){df <- df[order(index),]}
  if(reverse == TRUE){df <- df[order(-index),]}
  return(df)
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "seagreen2", forcat_line = "seagreen4", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  if(is.character(y_hist)){y_hist <- as.factor(y_hist)}
  if(is.character(y_forcat)){y_forcat <- factor(y_forcat, levels = levels(y_hist))}
  if(is.character(lower)){lower <- factor(lower, levels = levels(y_hist))}
  if(is.character(upper)){upper <- factor(upper, levels = levels(y_hist))}

  n_class <- NULL
  if(is.factor(y_hist)){class_levels <- levels(y_hist); n_class <- length(class_levels)}

  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = as.numeric(c(y_hist, y_forcat)))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = as.numeric(y_forcat))

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- as.numeric(lower); forcat_data$upper <- as.numeric(upper)}

  plot <- ggplot()+ geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  if(is.null(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)}
  if(is.numeric(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), breaks = 1:n_class, labels = class_levels)}
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}

