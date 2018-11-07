#' Print to STDOUT.
#' 
#' \code{sum} print using C-style printf syntax.
printf <- function(...) invisible(cat(sprintf(...)))

#' Run benchmark on function fun n times.
#'
#' @param n repeat n times
#' @param fun function to repeat
#' @param verbose displays more info when TRUE
bench = function(n, fun, verbose=FALSE) {
  digits = as.integer(log10(n)) + 1
  tests = rep(NA, length=n)
  progress_str = sprintf("\rExecuting test %%%dd/%%%dd... average: %%%d.1f", digits, digits, digits)
  
  start.time = Sys.time()
  
  average = 0
  for (i in 1:n) {
    capture.output({
      tests[i] = fun()
    })
    
    if (i %% 10 == 0) {
      average = mean(tests, na.rm=TRUE)
    }
    
    printf(progress_str, i, n, average)
  }
  
  total.time = Sys.time() - start.time
  average.time = total.time / n
  
  printf(" done!\n-------\n")
  printf("Total time: %s\n", format(total.time))
  printf("Average time per test: %s\n", format(average.time))
  if (any(is.na(tests)) > 0)
    printf("Failed runs: %d\n", sum(is.na(tests)))
  if (verbose) {
    printf("Test results (%d): ", n); cat(tests, sep=', ')
  }
  
  plotDensity(tests)
  
  tests
}

#' Compare execution times between basicDM and aStarDM.
#' 
#' Note: Both BasicDM and AStarPQueueDM will use the same seed.
#'
#' @param n Number of tries.
#' @param del Number of deliveries.
#' @param dim 
#' @param seed Seed value to use.
compareDMs <- function(n=10, del=5, dim=10, seed=NULL, verbose=TRUE) {
  if (is.null(seed)) {
    seed = runif(1, min=-2e9, max=2e9)
  }
  
  seed = as.integer(seed)
  
  printRun <- function(title, carReady) {
    set.seed(seed)
    printf("---- %s ----\n", title)
    bench(n, function() runDeliveryMan(carReady=carReady, del=del, dim=dim, pause=0, doPlot=0), verbose)
  }
  
  printf("==================\n")
  printf("n:    %d\ndel:  %d\ndim:  %d\nseed: %d\n", n, del, dim, seed)
  printRun("BasicDM", basicDM)
  printRun("AStarPQueueDM", aStarPQueueDM)
}

calculateMinMaxInterval = function(n, confidence_level, probs_fun) {
  rows = nrow(probs_fun()$salinity)
  min_values = vector(mode="numeric", length=n)
  max_values = vector(mode="numeric", length=n)
  
  for (i in 1:n) {
    p = probs_fun()
    min_values[i] = min(apply(p$salinity, 1, function(row) row[1] - qnorm(confidence_level) * (row[2] / sqrt(rows))),
                        apply(p$phosphate, 1, function(row) row[1] - qnorm(confidence_level) * (row[2] / sqrt(rows))),
                        apply(p$nitrogen, 1, function(row) row[1] - qnorm(confidence_level) * (row[2] / sqrt(rows))))
    max_values[i] = max(apply(p$salinity, 1, function(row) row[1] + qnorm(confidence_level) * (row[2] / sqrt(rows))),
                        apply(p$phosphate, 1, function(row) row[1] + qnorm(confidence_level) * (row[2] / sqrt(rows))),
                        apply(p$nitrogen, 1, function(row) row[1] + qnorm(confidence_level) * (row[2] / sqrt(rows))))
    if (i %% 1000 == 0) { printf("%d/%d\n", i, n) }
  }
  
  list(min=list(min=min(min_values), mean=mean(min_values)),
       max=list(max=max(max_values), mean=mean(max_values)))
}

averageTest <- function(tests, fun=function() runDeliveryMan(carRead=aStarPQueueDM, dim=10, turns=2000, doPlot=F, pause=0, del=5)){
  start.time = Sys.time()
  sum = 0
  for (i in 1:tests) {
    sum = sum + fun()
    if (i %% 10 == 0){
      print(i)
      print(sum / i)
    }
  }
  print(sum / i)
  total = Sys.time() - start.time
  printf("Total time: %s\n", format(total))
  return(0)
}

plotDensity = function(tests) {
  require(ggplot2)
  
  N = length(tests)
  
  d = data.frame(count = factor(rep(c("A", "B"), each=N)),
                 moves = c(tests))

  ggplot(d, aes(x=moves)) +
    geom_histogram(binwidth=.5,
                   colour="black",
                   fill="white") +
    geom_vline(aes(xintercept=mean(moves, na.rm=T)),   # Ignore NA values for mean
               color="red", linetype="dashed", size=1)
}