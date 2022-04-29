############################
# JOINT ESTIMATION FUNCTIONS
############################

#### M=1 RECONSTRUCTION
ancestor_pair_estimation <- function(L, m, pairs_est, ones_est){
  # define matrix of needed pair estimates
  pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
  for (i in 1:(L-m)){
    pairs_needed[,i] <- c(rep(pairs_est[1,i], times = (1/2^(m+i))*2^L), 
                          rep(pairs_est[2,i], times = (1/2^(m+i))*2^L), 
                          rep(pairs_est[3,i], times = (1/2^(m+i))*2^L), 
                          rep(pairs_est[4,i], times = (1/2^(m+i))*2^L))
  }
  
  # define matrix of needed onewise estimates
  ones_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
  for (j in 1:(L-m-1)){
    ones_needed[,j] <- rep(c(rep(ones_est[1,(j+1)], times = (1/2^(m+j))*2^L),
                             rep(ones_est[2,(j+1)], times = (1/2^(m+j))*2^L)), times = 2)
  }
  
  # concatenate into one matrix
  products <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
  products[,1:(L-m)] <- pairs_needed
  products[,(L-m+1):(2*L - 2*m - 1)] <- (1/ones_needed)
  
  # calculate estimates and replace NA's with 0
  ancestor_estimates <- matrix(NA, nrow = (2^L), ncol = 1)
  ancestor_estimates <- apply(products, 1, prod)
  ancestor_estimates[is.na(ancestor_estimates)] <- 0
  ancestor_estimates <- as.matrix(ancestor_estimates)
  
  # return matrix of estimates
  return(ancestor_estimates)
}


##### M=2 RECONSTRUCTION
ancestor_three_estimation <- function(L, m, three_est, pairs_est){
  
  # define matrix of needed threewise estimates
  threes_needed <- matrix(NA, nrow = (2^L) , ncol = (L-m))
  for (i in 1:(L-m)){
    threes_needed[,i] <- c(rep(three_est[1,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[2,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[3,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[4,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[5,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[6,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[7,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[8,i], times = (1/2^(m+i))*2^L))
  }
  
  # define matrix of needed pairwise estimates
  pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
  for(j in 1:(L-m-1)){
    pairs_needed[,j] <- c(rep(pairs_est[1,(j+1)], times = (1/2^(m+j))*2^L),
                          rep(pairs_est[2,(j+1)], times = (1/2^(m+j))*2^L),
                          rep(pairs_est[3,(j+1)], times = (1/2^(m+j))*2^L),
                          rep(pairs_est[4,(j+1)], times = (1/2^(m+j))*2^L))
  }
  
  # concatenate into one matrix and invert pairs
  to_mult <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m -1))
  to_mult[,(1:(L-m))] <- threes_needed
  to_mult[,((L-m+1):(2*L - 2*m - 1))] <- 1/(pairs_needed)
  
  # calculate estimates
  an_est <- matrix(NA, nrow = (2^L), ncol = 1)
  an_est <- apply(to_mult, 1, prod)
  
  # replace NA's with 0
  an_est[is.na(an_est)] <- 0
  
  # return matrix of elements
  an_est <- as.matrix(an_est)
  return(an_est)
}


##### M=3 RECONSTRUCTION
ancestor_four_estimation <- function(L, m, four_est, three_est){
  # initialize matrix for needed fourwise estimates
  fours_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
  
  # fill matrix with needed fourwise estimates
  for (i in 1:(L-m)){
    fours_needed[,i] <- c(rep(four_est[1,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[2,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[3,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[4,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[5,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[6,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[7,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[8,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[9,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[10,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[11,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[12,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[13,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[14,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[15,i], times = (1/2^(m+i))*2^L),
                          rep(four_est[16,i], times = (1/2^(m+i))*2^L))
  }
  
  # initialize matrix for needed threewise estimates
  threes_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
  
  # fill matrix with needed threewise estimates
  for (j in 1:(L-m-1)){
    threes_needed[,j] <- c(rep(three_est[1, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[2, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[3, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[4, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[5, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[6, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[7, (j+1)], times = (1/2^(m+j))*2^L),
                           rep(three_est[8, (j+1)], times = (1/2^(m+j))*2^L))
  }
  
  # concatenate into one matrix & take inverse of threewise
  to_mult <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
  to_mult[,(1:(L-m))] <- fours_needed
  to_mult[,((L-m+1):(2*L - 2*m - 1))] <- 1/(threes_needed)
  
  # calculate estimates
  an_ests <- matrix(NA, nrow = (2^L), ncol = 1)
  an_ests <- apply(to_mult, 1, prod)
  
  # replace NA's with zeros
  an_ests[is.na(an_ests)] <- 0
  
  # return matrix of estimates
  an_ests <- as.matrix(an_ests)
  return(an_ests)
}