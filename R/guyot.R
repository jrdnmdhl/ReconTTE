#' @export
prep_data <- function(surv, nrisk) {
  cleaned_surv <- clean_digitized_survival(surv, nrisk)
  cleaned_n_risk <- clean_n_risk(nrisk, cleaned_surv)
  list(survival = cleaned_surv, n_risk = cleaned_n_risk)
}

# Clean up digitized survival data to ensure that survival does not increase
# and time does not decrease.
clean_digitized_survival <- function(surv, nrisk) {

  # Make sure that dataset includes time zero, sort it, remove duplicates,
  # and number it.
  temp_surv <- surv %>%
    select(time, survival) %>%
    rbind(data.frame(time = nrisk$time, survival = 1)) %>%
    distinct(time, .keep_all = T) %>%
    arrange(time) %>%
    mutate(t = time, i = seq_len(n()), s = 1)

  # Iterate over each row and ensure survival does not increase
  n_row <- nrow(temp_surv)
  for (i in seq_len(n_row - 1)) {
    temp_surv$s[i+1] <- min(temp_surv$s[i], temp_surv$survival[i+1])
  }

  # Seelct only necessary columns
  outSurv <- select(temp_surv, i, t, s) %>%
    mutate(p_death = c(0, (1 - s / lag(s))[-1])) %>%
    as_tibble()

  return(outSurv)
}

# Clean up numbers at risk
clean_n_risk <- function(nrisk, surv) {

  index_n_risk <- nrisk %>%
    rowwise() %>%
    group_split() %>%
    map(function(d){
      prior_surv = surv %>% filter(surv$t >= d$time)
      nRows = nrow(prior_surv)
      d$lower = head(prior_surv, 1)$i
      return(d)
    }) %>%
    bind_rows() %>%
    mutate(
      upper = ifelse(
        seq_len(n()) == max(n()),
        max(surv$i),
        lead(lower) - 1
      ),
      t = time,
      i = seq_len(n())
    ) %>%
    mutate(
      n_start = n,
      n_end = c(n[-1], 0)
    ) %>%
    select(i, t, lower, upper, n_start, n_end)

  return(index_n_risk)
}
