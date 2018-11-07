#' Print to STDOUT.
#' 
#' \code{sum} print using C-style printf syntax.
printf <- function(...) invisible(cat(sprintf(...)))

#' Run a function n times and measure the total and average execution times.
#'
#' @param n repeat n times
#' @param fun function to repeat
timer <- function(n, fun, verbose=TRUE) {
  runs = NULL
  failed = NULL
  
  start.time = Sys.time()
  
  capture.output({
    for (i in 1:n) {
      t = fun()
      if (is.na(t)) failed = c(failed, i)
      else runs = c(runs, t)
    }
  })
  
  end.time = Sys.time()
  
  total = end.time - start.time
  average = total / n
  
  if (verbose) {
    printf("Total time: %s\n", format(total))
    if (length(failed) > 0)
      printf("Failed runs: %d\n", length(failed))
    printf("Runs (%d): ", n); cat(runs, sep=', ')
    printf("\nAverage turns: %g\n", round(mean(runs)))
    printf("Average time per run: %s\n", format(average))
  } else {
    printf("Avg. turns: %d  (%s / run)\n", round(mean(runs)), format(average))
  }
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
    timer(n, function() runDeliveryMan(carReady=carReady, del=del, dim=dim, pause=0, doPlot=0), verbose)
  }
  
  printf("==================\n")
  printf("n:    %d\ndel:  %d\ndim:  %d\nseed: %d\n", n, del, dim, seed)
  printRun("BasicDM", basicDM)
  printRun("AStarPQueueDM", aStarPQueueDM)
}

averageTest <- function(tests){
  start.time = Sys.time()
  sum = 0
  for (i in 1:tests) {
    sum=sum+runDeliveryMan(carReady = aStarPQueueDM, dim = 10, turns = 2000, doPlot = F, pause = 0, del = 5)
    if(i%%10==0){
      print(i)
      print(sum/i)
    }
  }
  print(sum/i)
  total = Sys.time() - start.time
  printf("Total time: %s\n", format(total))
  return(0)
}