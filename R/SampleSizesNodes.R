#' Sample sizes in each node
#'
#' @description The function computes the releavant sample sizes in each node of
#' the given partition.
#'
#' @param stratsd a vector specifying the standard deviations of each stratum of the
#' partition.
#'
#' @param n the desired sample size.
#'
#' @param stratsize a vector specifying the number of observations in the population
#' in each stratum of the partition.
#'
#' @param thhold a number indicating the minimum number of observations to be selecting 
#' from each stratum.
#'
#' @return A named vector with the same size as the \code{stratsd} and \code{stratsize} with
#' the number of observations to be selected from each stratum in the sample of size \code{n}.
#'

#'
#' @keywords internal


SampleSizesNodes<- function(stratsd, n, stratsize, thhold){
  
  nstrata = length(stratsize)
  popsize = sum(stratsize)
  #With the ingredients above, we like a sample of size n so that 
  # (1) n_i units are selected from stratum i
  # (2) n_i's sum to n
  # (3) n_i's are (approximately) proportional to sqrt(stratsize)*stratsd
  # (4) except that the n_i's should at least be equal to a certain threshold
  # (5) and that n_i cannot exceed the i-th stratum size
  
  # Let's insure that the n_i's are at least equal to a certain threshold unless
  # the threshold exceeds the stratum size
  
  # The target for the stratum sample sizes (although not integers)
  target = n*sqrt(stratsize)*stratsd/sum(sqrt(stratsize)*stratsd) 
  
  samplesize = round(target) # Now the numbers are integers, but may not satisfy 
  # (2), (4) and (5)
  
  # We start by fixing (4) and (5)
  # thhold = 1 # One could also make the threshold a function of n
  for (i in 1:nstrata){
    if (samplesize[i] < thhold) {samplesize[i] = min(thhold, stratsize[i])}
    if (samplesize[i] > stratsize[i]) {samplesize[i] = stratsize[i]}
  }
  
  # The sum of the stratum sample sizes could now be more or less than n
  # We need to make some (typically) minor adjustments by increasing or
  # decreasing some of the sample sizes
  # We cannot increase those that are already at the stratum size, and we
  # cannot decrease those at the threshold size
  # We will decide which sizes to increase or decrease based on the
  # discrepancies with the target values
  
  discr = (samplesize - target)/target
  # print(target)
  #print("samplesize")
  #print(samplesize)
  if (sum(samplesize) - n > 0) {decrease = sum(samplesize) - n
  for (i in 1:decrease) {
    label = which(samplesize > thhold)
    samplesize[which(discr == max(discr[label]))]=samplesize[which(discr == max(discr[label]))]-1
    discr = (samplesize - target)/target}
  }
  if (sum(samplesize) - n < 0) {increase = n - sum(samplesize)
  for (i in 1:increase) {
    label = which(samplesize < stratsize)
    samplesize[which(discr == min(discr[label]))]=samplesize[which(discr == min(discr[label]))]+1
    discr = (samplesize - target)/target
  }
  }
  samplesize=as.numeric(samplesize)
  names(samplesize) <- names(stratsize)
  return(samplesize)
} 