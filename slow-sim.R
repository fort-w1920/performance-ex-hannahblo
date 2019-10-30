simulate <- function(reps, seed, data, true_coef = 0:ncol(data), df = 4) {
  set.seed(seed)
  coefs <- NULL
  for (rep in seq_len(reps)) {
    coefs <- cbind(coefs, simulate_once(data, true_coef, df))
  }
  return(structure(coefs, seed = seed))
}

simulate_once <- function(data, true_coef, df) {
  data <- simulate_response(data, true_coef, df)
  estimate_coef(data)
}

simulate_response <- function(data, true_coef, df) {
  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef
  data$y <- expected + rt(nrow(data), df = df)
  data
}

estimate_coef <- function(data) {
  model <- lm(y ~ ., data = data)
  unname(coef(model))
}
