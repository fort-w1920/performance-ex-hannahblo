
############### a #############################################################
# Ausführung von:
# library(profvis)
# profvis(test <- simulate(reps = 100, seed = 20141028, data = testdata))

# Die meiste Zeit wird in Zeile 5 gebraucht bei der Berechnung des Linearen
# Models verbraucht.





############## b ##############################################################
library(parallel)
simulate_opt <- function(reps, seed, data, true_coef = 0:ncol(data), df = 4) {
  
  # @ reps (count): number of repetition
  # @ seed (count): seed of the simulation
  # @ data (dataframe): data
  # @ true_coef (numeric): ture underlying coefficients
  # @ df (count): degree of freedem (caluclation of t-distribution)
  
  # This function calcuates a simulation study
  
  # Return: Simulated data 
  
  # Input checking
  checkmate::check_count(reps)
  checkmate::check_count(seed)
  checkmate::check_data_frame(data, any.missing = FALSE, 
                              types = c("numeric"),
                              ncols = length(true_coef) - 1)
  checkmate::check_numeric(true_coef)
  checkmate::check_count(df)
  
  
  # set seed
  set.seed(seed)
  
  # Design matrix and the expected value are fixed for every simulation loops
  design <- model.matrix(~., data = data)
  expected <- design %*% true_coef
  solve_matrix <- tcrossprod(solve(crossprod(design, design)), design)
  
  # Calculate the number of cores ('-1' to prevent freezing)
  no_cores <- detectCores() - 1
  
  if (no_cores > 1) {
    # inital cluster
    cl <- makeCluster(no_cores)
    # calculate coeficients
    coefs <- parallel::parLapply(X = seq_len(reps), fun = simulate_once_opt,
                                 data_row = nrow(data),
                                 solve_matrix = solve_matrix,
                                 expected = expected, df = df, cl = cl)
    # stop cluster
    parallel::stopCluster(cl)
  }
  if (!(no_cores > 1)) {
    coefs <- lapply(X = seq_len(reps), fun = simulate_once_opt,
                    data_row = nrow(data), design = design,
                    expected = expected, df = df)
  }
  
  return(matrix(unlist(coefs), ncol = reps))
}

simulate_once_opt <- function(rep, data_row, solve_matrix, expected, df = 4) {
  # @rep (count): for lapply
  # @data_row (count): number of rows
  # @solve_matrix (matrix): solved system of equation using design matrix 
  # @expected (vector): expected values
  # @df (count): degree of freedem (caluclation of t-distribution)
  
  # Return: estimated values
  
  y <- rt(data_row, df = df) + expected
  return(solve_matrix %*% y)
}




#### Comparison of both implementations
source("slow-sim.R")
set.seed <- 232323
observations <- 5000
covariates <- 10
testdata <- as.data.frame(matrix(rnorm(observations * covariates),
                                 nrow = observations
))

# Possibility 1.
system.time(test <- simulate(reps = 100, seed = 20141028, data = testdata))
system.time(test_opt <- simulate_opt(reps = 100, seed = 20141028, 
                                     data = testdata))

# Possibility 2.
rbenchmark::benchmark(simulate(reps = 100, seed = 20141028, data = testdata),
                      simulate_opt(reps = 100, seed = 20141028, 
                                   data = testdata),
                      replications = 10, 
                      columns = c("test", "elapsed", "relative"),
                      order = "elapsed"
)

# Possibility 3.
library(profvis)
profvis(simulate(reps = 100, seed = 20141028, data = testdata))
profvis(simulate_opt(reps = 100, seed = 20141028, data = testdata))
