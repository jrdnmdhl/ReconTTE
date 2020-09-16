
###FUNCTION INPUTS
#' @export
guyot <- function(survival, n_risk, n_events = 'NA') {

  survival <- clean_digitized_survival(survival, n_risk)
  n_risk <- clean_n_risk(n_risk, survival)
  n_points <- nrow(survival)
  res <- tibble(
    n_censors = rep(0, n_points),
    n_events = rep(0, n_points),
    n_risk_start = rep(n_risk$n_start[1], n_points),
    n_risk_end = rep(n_risk$n_start[1], n_points),
    survival = rep(1, n_points)
  )
  n_risk %>%
    rowwise() %>%
    group_split() %>%
    reduce(solve_interval, .init = res, survival = survival)
}

solve_interval <- function(results_so_far, interval, survival) {
  #if (interval$upper == nrow(survival)) return(results_so_far)
  target_nrisk <- interval$n_end
  estimated_nrisk <- -1
  censor_guess <- 0
  guesses <- c()
  while(target_nrisk != estimated_nrisk) {
    if (censor_guess %in% guesses) {
      stop('guess already tried')
    }
    guesses <- append(guesses, censor_guess)
    guess_result <- guess_interval(results_so_far, interval, survival, censor_guess)
    estimated_nrisk <- guess_result$n_risk_end[interval$upper]
    censor_guess <- max(censor_guess + (estimated_nrisk - target_nrisk), 0)
  }
  guess_result
}

guess_interval <- function(results_so_far, interval, survival, n_censors) {
  lower <- interval$lower
  upper <- interval$upper
  last_interval <- upper == nrow(results_so_far)
  interval_start_time <- survival$t[lower]
  if (last_interval) {
    interval_end_time <- survival$t[upper]
    interval_breaks <- survival$t[lower:upper]
    interval_selector <- lower:(upper - 1)
  } else {
    interval_end_time <- survival$t[upper + 1]
    interval_breaks <- survival$t[lower:(upper+1)]
    interval_selector <- lower:upper
  }
  censor_times <- survival$t[lower] + seq_len(n_censors) * (interval_end_time - interval_start_time) / (n_censors + 1)
  if (n_censors == 0) {
    results_so_far$n_censors[interval_selector] <- 0
  } else {
    results_so_far$n_censors[interval_selector] <- hist(censor_times, breaks = interval_breaks, plot = F)$counts
  }
  n_points <- nrow(survival)
  for (i in lower:upper) {
    results_so_far$n_events[i] <- round(results_so_far$n_risk_start[i] * survival$p_death[i])
    estimated_p_death <- results_so_far$n_events[i] / results_so_far$n_risk_start[i]
    if (i == 1) {
      prior_survival <- 1
    } else {
      prior_survival <- results_so_far$survival[i - 1]
    }
    results_so_far$survival[i] <- prior_survival * (1 - estimated_p_death)
    results_so_far$n_risk_end[i] <- results_so_far$n_risk_start[i] - results_so_far$n_censors[i] - results_so_far$n_events[i]
    if (i < n_points) {
      results_so_far$n_risk_start[i + 1] <- results_so_far$n_risk_end[i]
    }
  }
  results_so_far
}
