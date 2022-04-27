##############################
# MARGINAL ESTIMATION FUNCTIONS
##############################

################
# MAIN FUNCTIONS
################

### M=0 MARGINAL
########
estimates_m0 <- function(descendents, L, n){
  onewise_margins <- matrix(NA, nrow = 2, ncol = L)
  
  
  rownames(onewise_margins) <- paste("pi_",c(0,1), sep="")
  colnames(onewise_margins) <- paste("site_", (1:(L)), sep="")
  for (i in 1:(dim(descendents)[2])){
    onewise_margins[1,i] <- (n - sum(descendents[,i]))/n
    onewise_margins[2,i] <- sum(descendents[,i])/n
  }
  return(onewise_margins)
}
########

### M=1 MARGINAL
#######
estimates_m1 <- function(L, q, n, d){
  # fix lower order marginals
  pi_m0 <- estimates_m0(d, L, n) 
  
  # initialize matrix for pairwise estimates; this will be filled in and returned
  pairwise_margins <- matrix(NA, nrow = 4, ncol = (L-1))
  rownames(pairwise_margins) <- paste("pi_", c("0,0", "0,1" , "1,0" , "1,1") )
  colnames(pairwise_margins) <- paste( "site_" , seq(1:(L-1)))
  
  # iterate over the pairs
  for (j in 1:(L-1)){
    # select relevant onewise marginals, pi_s(0)
    pi_zeroes <- pi_m0[1, j:(j+1)]
    
    
    # calculate upper and lower bounds
    l_bound <- max(0, (sum(pi_zeroes)-1))
    u_bound <- min(pi_zeroes)
    
    #check if there is space bw bounds; if not set estimate to be avg of the bounds
    if (abs(l_bound - u_bound) <= 10^(-10)){
      phi <- (l_bound + u_bound)/2
    }
    if (abs(l_bound - u_bound) > 10^(-10)){
      # calculate unconstrained MLE
      n_00 <- sum(d[,j]==0 & d[,j+1]==0)
      mle <- (1 / (1-q))*((n_00/n) - q*pi_zeroes[1]*pi_zeroes[2])
      
      # find constrained estimate
      phi <- ((l_bound <= mle) & (mle <= u_bound))*mle + l_bound*(mle < l_bound) + u_bound*(mle > u_bound)
    }
    
    # calculate pi's from phi
    pi_00 <- phi
    pi_01 <- pi_zeroes[1] - phi
    pi_10 <- pi_zeroes[2] - phi
    pi_11 <- 1 - pi_zeroes[1] - pi_zeroes[2] + phi
    
    #attach to matrix
    pairwise_margins[,j] <- round(c(pi_00, pi_01, pi_10, pi_11), 10)
  }
  
  return (pairwise_margins)
}
####### 

### M=2 MARGINAL
#######
# write a function to calculate the threewise marginal estimates
estimates_m2 <- function(L, q, n, descend){
  # fix onewise estimates
  pi_one <- estimates_m0(descend, L,n)
  # fix pairwise estimates
  pi_two <- estimates_m1(L, q, n, descend)
  
  # initialize matrix to put the threewise margin estimates into
  threewise_margins <- matrix(NA, nrow = 8, ncol = (L-2))
  
  ## all m+1-wise ancestral sequences
  an_seq <- matrix(c(0:(2^(2+1)-1)), ncol=1)
  an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=2+1) ) 
  an_seq <- matrix( as.character(an_seq), nrow=2^{2+1}, ncol=2+1)
  an_seq <- apply(an_seq, 1, paste, collapse="")
  
  rownames(threewise_margins) <- paste("pi_", an_seq ,sep="")
  colnames(threewise_margins) <- paste("site_",seq(1:(L-2)),sep="")
  
  #### ESTIMATE PHIS ####
  # iterate over groups of 3 sites
  for (j in 1:(L-2)){
    pi_ones <- pi_one[,j:(j+2)]
    pi_pairs <- pi_two[,j:(j+1)]
    
    ### first estimate phi(0) ###
    # calculate n(xs, 0, xs+2)
    n_000 <- sum(descend[,j]==0 & descend[,j+1]==0 & descend[,j+2]==0)
    n_001 <- sum(descend[,j]==0 & descend[,j+1]==0 & descend[,j+2]==1)
    n_100 <- sum(descend[,j]==1 & descend[,j+1]==0 & descend[,j+2]==0)
    n_101 <- sum(descend[,j]==1 & descend[,j+1]==0 & descend[,j+2]==1)
    n_0 <- c(n_000, n_001, n_100, n_101)
    
    # find bounds for phi(0)
    lower_phi0 <- max(0, (pi_pairs[1,1] + pi_pairs[1,2] - pi_ones[1,2]))
    upper_phi0 <- min(pi_pairs[1,1], pi_pairs[1,2])
    
    # estimate phi(0)
    if (abs(lower_phi0 - upper_phi0) <= 10^(-10)){
      hat_phi0 <- (lower_phi0 + upper_phi0)/2
    }
    if (abs(lower_phi0 - upper_phi0) > 10^(-10)){
      hat_phi0 <- roots_threewise_0(q, n_0, pi_ones, pi_pairs, lower_phi0, upper_phi0)
    }
    
    # calculate other pi's based on estimated phi
    pi_000 <- hat_phi0
    pi_001 <- round(pi_pairs[1,1] - hat_phi0, 10)
    pi_100 <- round(pi_pairs[1,2] - hat_phi0, 10)
    pi_101 <- round(pi_ones[1,2] - pi_pairs[1,2] - pi_pairs[1,1] + hat_phi0, 10)
    
    ### now estimate phi(1) ###
    # calculate the n(xs, 1, xs+2)
    n_010 <- sum(descend[,j]==0 & descend[,j+1]==1 & descend[,j+2]==0)
    n_011 <- sum(descend[,j]==0 & descend[,j+1]==1 & descend[,j+2]==1)
    n_110 <- sum(descend[,j]==1 & descend[,j+1]==1 & descend[,j+2]==0)
    n_111 <- sum(descend[,j]==1 & descend[,j+1]==1 & descend[,j+2]==1)
    n_1 <- c(n_010, n_011, n_110, n_111)
    
    # find bounds for phi(1)
    lower_phi1 <- max(0, (pi_pairs[2,1] + pi_pairs[3,2] - pi_ones[2,2]))
    upper_phi1 <- min(pi_pairs[2,1], pi_pairs[3,2])
    
    #estimate phi(1)
    if (abs(lower_phi1 - upper_phi1) <= 10^(-10)){
      hat_phi1 <- (lower_phi1 + upper_phi1)/2
    }
    if (abs(lower_phi1 - upper_phi1) > 10^(-10)){
      hat_phi1 <- roots_threewise_1(q, n_1, pi_ones, pi_pairs, lower_phi1, upper_phi1)
    }
    
    # calculate other pi's
    pi_010 <- hat_phi1
    pi_011 <- round(pi_pairs[2,1] - hat_phi1, 10)
    pi_110 <- round(pi_pairs[3,2] - hat_phi1, 10)
    pi_111 <- round(pi_ones[2,2] - pi_pairs[3,2] - pi_pairs[2,1] + hat_phi1, 10)
    
    ### append estimates into matrix 
    threewise_margins[,j] <- c(pi_000,pi_001,pi_010,pi_011,pi_100,pi_101,pi_110,pi_111)
    
  }
  return(threewise_margins)
}
#######

##M=3 MARGINAL
#######
estimates_m3 <- function(L, q, n, d){
  # fix onewise
  onewise <- estimates_m0(d, L, n)
  # fix pairwise
  pairwise <- estimates_m1(L, q, n , d)
  # fix threewise
  threewise <- estimates_m2(L, q, n, d)
  
  # initialize matrix to hold fourwise estimates
  fourwise_marginals <- matrix(NA, nrow = 16, ncol = L-3)
  
  ## all m+1-wise ancestral sequences
  an_seq <- matrix(c(0:(2^(3+1)-1)), ncol=1)
  an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=3+1) ) 
  an_seq <- matrix( as.character(an_seq), nrow=2^{3+1}, ncol=3+1)
  an_seq <- apply(an_seq, 1, paste, collapse="")
  
  rownames(fourwise_marginals) <- paste("pi_", an_seq ,sep="")
  colnames(fourwise_marginals) <- paste("site_",seq(1:(L-3)),sep="")
  
  # iterate over the fourwise combos
  for (s in 1:(L-3)){
    # find relevant onewise
    pi_one <- onewise[,s:(s+3)]
    
    # find relevant pairwise
    pi_two <- pairwise[,s:(s+2)]
    
    # find relevant threewise
    pi_three <- threewise[,s:(s+1)]
    
    ### FIND PHI(0,0)=PI_{S,...,S+3}(0,0,0,0) ###
    
    # calculate num observations for ea sequence
    n_0000 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==0) #pi(0,0,0,0)
    n_0001 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==1) #pi(0,0,0,1)
    n_1000 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==0) #pi(1,0,0,0)
    n_1001 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==0 & d[,s+3]==1) #pi(1,0,0,1)
    n_gr1 <- c(n_0000, n_0001, n_1000, n_1001)
    
    # calculate bounds on phi(0,0)
    lower_gr1 <- max(0, (pi_three[1,2]+ pi_three[1,1] - pi_two[1,2]))
    upper_gr1 <- min(pi_three[1,2], pi_three[1,1])
    
    # estimate phi(0,0)
    if (abs(lower_gr1 - upper_gr1) <= 10^(-10)){
      hat_phi_gr1 <- (lower_gr1 + upper_gr1)/2
    }
    if (abs(lower_gr1 - upper_gr1) > 10^(-10)){
      hat_phi_gr1 <- roots_fourwise_gr1(q, n_gr1, pi_one, pi_two, pi_three, lower_b = lower_gr1, upper_b = upper_gr1)
    }
    
    # estimate other pi's in group 1 based on constrained MLE for phi(0,0)
    pi_0000 <- hat_phi_gr1
    pi_0001 <- round((pi_three[1,1] - hat_phi_gr1), 10)
    pi_1000 <- round((pi_three[1,2] - hat_phi_gr1), 10)
    pi_1001 <- round((pi_two[1,2] - pi_three[1,1] - pi_three[1,2] + hat_phi_gr1), 10)
    
    
    ### FIND PHI(0,1)=PI_{S,...,S+3}(0,0,1,0) ###
    
    # calculate num observations for ea sequence
    n_0010 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==0) #pi(0,0,1,0)
    n_0011 <- sum(d[,s]==0 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==1) #pi(0,0,1,1)
    n_1010 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==0) #pi(1,0,1,0)
    n_1011 <- sum(d[,s]==1 & d[,s+1]==0 & d[,s+2]==1 & d[,s+3]==1) #pi(1,0,1,1)
    n_gr2 <- c(n_0010, n_0011, n_1010, n_1011)
    
    # calculate bounds on phi(0,1)
    lower_gr2 <- max(0, (pi_three[3,2] + pi_three[2,1] - pi_two[2,2]))
    upper_gr2 <- min(pi_three[2,1], pi_three[3,2])
    
    # estimate phi(0,1)
    if (abs(lower_gr2 - upper_gr2) <= 10^(-10)){
      hat_phi_gr2 <- (lower_gr2 + upper_gr2)/2
    }
    if (abs(lower_gr2 - upper_gr2) > 10^(-10)){
      hat_phi_gr2 <- roots_fourwise_gr2(q, n_gr2, pi_one, pi_two, pi_three, lower_gr2, upper_gr2)
    }
    
    # calculate estimates for group 2 based on constrained MLE of phi(0,1)
    pi_0010 <- hat_phi_gr2
    pi_0011 <- round((pi_three[2,1] - hat_phi_gr2), 10)
    pi_1010 <- round((pi_three[3,2] - hat_phi_gr2), 10)
    pi_1011 <- round((pi_two[2,2] - pi_three[2,1] - pi_three[3,2] + hat_phi_gr2), 10)
    
    ### FIND PHI(1,0)=PI_{S,...,S+3}(0,1,0,0) ###
    
    # calculate num observations for ea sequence
    n_0100 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==0) #pi(0,1,0,0)
    n_0101 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==1) #pi(0,1,0,1)
    n_1100 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==0) #pi(1,1,0,0)
    n_1101 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==0 & d[,s+3]==1) #pi(1,1,0,1)
    n_gr3 <- c(n_0100, n_0101, n_1100, n_1101)
    
    # calculate bounds on phi(1,0)
    lower_gr3 <- max(0, (pi_three[5,2] + pi_three[3,1] - pi_two[3,2]))
    upper_gr3 <- min(pi_three[3,1], pi_three[5,2])
    
    # estimate phi(1,0)
    if (abs(lower_gr3 - upper_gr3) <= 10^(-10)){
      hat_phi_gr3 <- (lower_gr3 + upper_gr3)/2
    }
    if (abs(lower_gr3 - upper_gr3) > 10^(-10)){
      hat_phi_gr3 <- roots_fourwise_gr3(q, n_gr3, pi_one, pi_two, pi_three, lower_gr3, upper_gr3)
    }
    
    # calculate estimates for group 3 based on constrained MLE for phi(1,0)
    pi_0100 <- hat_phi_gr3
    pi_0101 <- round((pi_three[3,1] - hat_phi_gr3), 10)
    pi_1100 <- round((pi_three[5,2] - hat_phi_gr3), 10)
    pi_1101 <- round((pi_two[3,2] - pi_three[3,1] - pi_three[5,2] + hat_phi_gr3), 10)
    
    ### FIND PHI(0,0)=PI_{S,...,S+3}(0,1,1,0) ###
    
    # calculate num observations for ea sequence
    n_0110 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==0) #pi(0,1,1,0)
    n_0111 <- sum(d[,s]==0 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==1) #pi(0,1,1,1)
    n_1110 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==0) #pi(1,1,1,0)
    n_1111 <- sum(d[,s]==1 & d[,s+1]==1 & d[,s+2]==1 & d[,s+3]==1) #pi(1,1,1,1)
    n_gr4 <- c(n_0110, n_0111, n_1110, n_1111)
    
    # calculate bounds on phi(1,1)
    lower_gr4 <- max(0, (pi_three[7,2] + pi_three[4,1] - pi_two[4,2]))
    upper_gr4 <- min(pi_three[4,1], pi_three[7,2])
    
    # estimate phi(1,1)
    if (abs(lower_gr4 - upper_gr4) <= 10^(-10)){
      hat_phi_gr4 <- (lower_gr4 + upper_gr4)/2
    }
    if (abs(lower_gr4 - upper_gr4) > 10^(-10)){
      hat_phi_gr4 <- roots_fourwise_gr4(q, n_gr4, pi_one, pi_two, pi_three, lower_gr4, upper_gr4)
    }
    
    # estimate others in group 4 based on constrained MLE phi(1,1)
    pi_0110 <- hat_phi_gr4
    pi_0111 <- round((pi_three[4,1] - hat_phi_gr4), 10)
    pi_1110 <- round((pi_three[7,2] - hat_phi_gr4), 10)
    pi_1111 <- round((pi_two[4,2] - pi_three[4,1] - pi_three[7,2] + hat_phi_gr4), 10)
    
    # append to fourwise matrix
    fourwise_marginals[,s] <- c(pi_0000, pi_0001, pi_0010, pi_0011, pi_0100, pi_0101, pi_0110, pi_0111, pi_1000, pi_1001, pi_1010, pi_1011, pi_1100, pi_1101, pi_1110, pi_1111)
  }
  return(fourwise_marginals)
}
#######

################
# SUB FUNCTIONS
################

## SubFunctions for m=2 Main Function

a_s <- function(q){
  # order: a_00, a_01, a_10, a_11
  a <- c((1-q)^2, -(1-q)^2, -(1-q)^2 , (1-q)^2)
  return(a)
}

b_s0 <- function(q, pi_one, pi_two){
  b_00 <- q*(1-q)*(pi_one[1,1]*pi_two[1,2] + pi_two[1,1]*pi_one[1,3]) + (q^2)*pi_one[1,1]*pi_one[1,2]*pi_one[1,3]
  b_01 <- ((1-q)^2)*pi_two[1,1] + q*(1-q)*(pi_one[1,1]*pi_two[2,2] + pi_two[1,1]*pi_one[2,3]) + (q^2)*pi_one[1,1]*pi_one[1,2]*pi_one[2,3]
  b_10 <- ((1-q)^2)*pi_two[1,2] + q*(1-q)*(pi_one[2,1]*pi_two[1,2] + pi_two[3,1]*pi_one[1,3]) + (q^2)*pi_one[2,1]*pi_one[1,2]*pi_one[1,3]
  b_11 <- ((1-q)^2)*(pi_one[1,2] - pi_two[1,1] - pi_two[1,2]) + q*(1-q)*(pi_one[2,1]*pi_two[2,2] + pi_two[3,1]*pi_one[2,3]) + (q^2)*pi_one[2,1]*pi_one[1,2]*pi_one[2,3]
  
  return(c(b_00, b_01, b_10, b_11))
}

A_1 <- function(n, a){
  A1 <- sum(n)*prod(a)
  return(A1)
}

A_2 <- function(n, a, b){
  t1 <- prod(a[1], a[2], a[3], b[4])*sum(n[1], n[2], n[3])
  t2 <- prod(a[1], a[2], b[3], a[4])*sum(n[1], n[2], n[4])
  t3 <- prod(a[1], b[2], a[3], a[4])*sum(n[1], n[3], n[4])
  t4 <- prod(b[1], a[2], a[3], a[4])*sum(n[2], n[3], n[4])
  A2 <- sum(c(t1,t2,t3,t4))
  
  return(A2)
}

A_3 <- function(n, a, b){
  t1 <- prod(a[1], a[2], b[3], b[4])* sum(n[1], n[2])
  t2 <- prod(a[1], b[2], a[3], b[4])* sum(n[1], n[3])
  t3 <- prod(a[1], b[2], b[3], a[4])* sum(n[1], n[4])
  t4 <- prod(b[1], a[2], a[3], b[4])* sum(n[2], n[3])
  t5 <- prod(b[1], a[2], b[3], a[4])* sum(n[2], n[4])
  t6 <- prod(b[1], b[2], a[3], a[4])* sum(n[3], n[4])
  A3 <- sum(c(t1, t2, t3, t4, t5, t6))
  
  return(A3)
}

A_4 <- function(n, a, b){
  t1 <- prod(a[1], b[2], b[3], b[4])*n[1]
  t2 <- prod(b[1], a[2], b[3], b[4])*n[2]
  t3 <- prod(b[1], b[2], a[3], b[4])*n[3]
  t4 <- prod(b[1], b[2], b[3], a[4])*n[4]
  A4 <- sum(c(t1, t2, t3, t4))
  
  return(A4)
}


pick_real = function(aaa){
  ## sequence of complex numbers to real numbers
  
  idx  = (abs(Im(aaa))>1e-4);
  return(Re(aaa)[!idx])
}


roots_threewise_0 <- function(q, n, pi_one, pi_two, lower_b, upper_b){
  
  # find a's and b's
  a_0 <- a_s(q)
  b_0 <- b_s0(q, pi_one, pi_two)
  
  # find A's
  A1_0 <- A_1(n, a_0)
  A2_0 <- A_2(n, a_0, b_0)
  A3_0 <- A_3(n, a_0, b_0)
  A4_0 <- A_4(n, a_0, b_0)
  
  # calculate the roots
  all_roots <- cubic(A1_0, A2_0, A3_0, A4_0)[c(1:3)]
  
  # find the real roots 
  real_roots <- pick_real(all_roots)
  
  # if there are 3 real roots, pick the middle one
  if (length(real_roots) == 1){
    candidate <- real_roots
  } else {
    candidate <- sort(real_roots)[2]
  }
  
  
  # make final choice of estimate based on bounds
  estimate <- candidate*(lower_b <= candidate & candidate <= upper_b) + upper_b*(candidate > upper_b) + lower_b*(candidate < lower_b)
  
  return(round(estimate, 10))
}


b_s1 <- function(q, pi_one, pi_two){
  b_00 <- q*(1-q)*(pi_one[1,1]*pi_two[3,2] + pi_two[2,1]*pi_one[1,3]) + (q^2)*pi_one[1,1]*pi_one[2,2]*pi_one[1,3]
  b_01 <- ((1-q)^2)*pi_two[2,1] + q*(1-q)*(pi_one[1,1]*pi_two[4,2] + pi_two[2,1]*pi_one[2,3]) + (q^2)*pi_one[1,1]*pi_one[2,2]*pi_one[2,3]
  b_10 <- ((1-q)^2)*pi_two[3,2] + q*(1-q)*(pi_one[2,1]*pi_two[3,2] + pi_two[4,1]*pi_one[1,3]) + (q^2)*pi_one[2,1]*pi_one[2,2]*pi_one[1,3]
  b_11 <- ((1-q)^2)*(pi_one[2,2] - pi_two[2,1] - pi_two[3,2]) + q*(1-q)*(pi_one[2,1]*pi_two[4,2] + pi_two[4,1]*pi_one[2,3]) + (q^2)*pi_one[2,1]*pi_one[2,2]*pi_one[2,3]
  
  return(c(b_00, b_01, b_10, b_11))
}


roots_threewise_1 <- function(q, n, pi_one, pi_two, lower_b, upper_b){
  
  # find a's and b's
  a_1 <- a_s(q)
  b_1 <- b_s1(q, pi_one, pi_two)
  # print(b_1)
  
  # find A's
  A1_1 <- A_1(n, a_1)
  A2_1 <- A_2(n, a_1, b_1)
  A3_1 <- A_3(n, a_1, b_1)
  A4_1 <- A_4(n, a_1, b_1)
  
  # calculate the roots
  all_roots <- cubic(A1_1, A2_1, A3_1, A4_1)[c(1:3)]
  
  # find the real roots 
  real_roots <- pick_real(all_roots)
  
  # if there are 3 real roots, pick the middle one
  if (length(real_roots) == 1){
    candidate <- real_roots
  } else {
    candidate <- sort(real_roots)[2]
  }
  
  # make final choice of estimate based on bounds
  estimate <- candidate*(lower_b <= candidate & candidate <= upper_b) + upper_b*(candidate > upper_b) + lower_b*(candidate < lower_b)
  
  return(round(estimate, 10))
}


#SubFunctions for m=3 Main Function

# this function is for the coefficients a_{xs, xs+2}^s 
## same function for all four phi's
a_s_four <- function(q){
  a_00 <- (1-q)^3
  a_01 <- -(1-q)^3
  a_10 <- -(1-q)^3
  a_11 <- (1-q)^3
  return(c(a_00, a_01, a_10, a_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(0,0)
b_s_gr1 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[1,2], pi_two[1,1]*pi_two[1,3], pi_three[1,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[1,3], pi_two[1,1]*pi_one[1,3]*pi_one[1,4], pi_one[1,1]*pi_two[1,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[1,3]*pi_one[1,4]))
  b_01 <- ((1-q)^3*pi_three[1,1] + 
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[2,2], pi_two[1,1]*pi_two[2,3], pi_three[1,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[2,3], pi_two[1,1]*pi_one[1,3]*pi_one[2,4], pi_one[1,1]*pi_two[1,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[1,3]*pi_one[2,4]))
  b_10 <- ((1-q)^3*pi_three[1,2] + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[1,2], pi_two[3,1]*pi_two[1,3], pi_three[5,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[1,3], pi_two[3,1]*pi_one[1,3]*pi_one[1,4], pi_one[2,1]*pi_two[1,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[1,3]*pi_one[1,4]))
  b_11 <- ((1-q)^3*(pi_two[1,2] - pi_three[1,1] - pi_three[1,2]) + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[2,2], pi_two[3,1]*pi_two[2,3], pi_three[5,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[2,3], pi_two[3,1]*pi_one[1,3]*pi_one[2,4], pi_one[2,1]*pi_two[1,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[1,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(0,1) 
b_s_gr2 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[3,2], pi_two[1,1]*pi_two[3,3], pi_three[2,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[3,3], pi_two[1,1]*pi_one[2,3]*pi_one[1,4], pi_one[1,1]*pi_two[2,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[2,3]*pi_one[1,4]))
  
  b_01 <- ((1-q)^3*pi_three[2,1] + 
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[4,2], pi_two[1,1]*pi_two[4,3], pi_three[2,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[1,2]*pi_two[4,3], pi_two[1,1]*pi_one[2,3]*pi_one[2,4], pi_one[1,1]*pi_two[2,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[1,2]*pi_one[2,3]*pi_one[2,4]))
  
  
  b_10 <- ((1-q)^3*pi_three[3,2] + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[3,2], pi_two[3,1]*pi_two[3,3], pi_three[6,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[3,3], pi_two[3,1]*pi_one[2,3]*pi_one[1,4], pi_one[2,1]*pi_two[2,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[2,3]*pi_one[1,4]))
  
  
  b_11 <- ((1-q)^3*(pi_two[2,2] - pi_three[2,1] - pi_three[3,2]) + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[4,2], pi_two[3,1]*pi_two[4,3], pi_three[6,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[1,2]*pi_two[4,3], pi_two[3,1]*pi_one[2,3]*pi_one[2,4], pi_one[2,1]*pi_two[2,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[1,2]*pi_one[2,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(1,0)
b_s_gr3 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[5,2], pi_two[2,1]*pi_two[1,3], pi_three[3,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[1,3], pi_two[2,1]*pi_one[1,3]*pi_one[1,4], pi_one[1,1]*pi_two[3,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[1,3]*pi_one[1,4]))
  
  
  b_01 <- ((1-q)^3*pi_three[3,1] + 
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[6,2], pi_two[2,1]*pi_two[2,3], pi_three[3,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[2,3], pi_two[2,1]*pi_one[1,3]*pi_one[2,4], pi_one[1,1]*pi_two[3,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[1,3]*pi_one[2,4]))
  
  b_10 <- ((1-q)^3*pi_three[5,2] + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[5,2], pi_two[4,1]*pi_two[1,3], pi_three[7,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[1,3], pi_two[4,1]*pi_one[1,3]*pi_one[1,4], pi_one[2,1]*pi_two[3,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[1,3]*pi_one[1,4]))
  
  
  b_11 <- ((1-q)^3*(pi_two[3,2] - pi_three[3,1] - pi_three[5,2]) + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[6,2], pi_two[4,1]*pi_two[2,3], pi_three[7,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[2,3], pi_two[4,1]*pi_one[1,3]*pi_one[2,4], pi_one[2,1]*pi_two[3,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[1,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}

# this function is for the constants b_{xs, xs+2}^s
## this function is for phi(1,1)
b_s_gr4 <- function(q, pi_one, pi_two, pi_three){
  b_00 <- (q*((1-q)^2)*sum(pi_one[1,1]*pi_three[7,2], pi_two[2,1]*pi_two[3,3], pi_three[4,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[3,3], pi_two[2,1]*pi_one[2,3]*pi_one[1,4], pi_one[1,1]*pi_two[4,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[2,3]*pi_one[1,4]))
  
  b_01 <- ((1-q)^3*pi_three[4,1] + 
             q*((1-q)^2)*sum(pi_one[1,1]*pi_three[8,2], pi_two[2,1]*pi_two[4,3], pi_three[4,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[1,1]*pi_one[2,2]*pi_two[4,3], pi_two[2,1]*pi_one[2,3]*pi_one[2,4], pi_one[1,1]*pi_two[4,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[1,1], pi_one[2,2]*pi_one[2,3]*pi_one[2,4]))
  
  b_10 <- ((1-q)^3*pi_three[7,2] + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[7,2], pi_two[4,1]*pi_two[3,3], pi_three[8,1]*pi_one[1,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[3,3], pi_two[4,1]*pi_one[2,3]*pi_one[1,4], pi_one[2,1]*pi_two[4,2]*pi_one[1,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[2,3]*pi_one[1,4]))
  
  
  b_11 <- ((1-q)^3*(pi_two[4,2] - pi_three[4,1] - pi_three[7,2]) + 
             q*((1-q)^2)*sum(pi_one[2,1]*pi_three[8,2], pi_two[4,1]*pi_two[4,3], pi_three[8,1]*pi_one[2,4]) +
             (q^2)*(1-q)*sum(pi_one[2,1]*pi_one[2,2]*pi_two[4,3], pi_two[4,1]*pi_one[2,3]*pi_one[2,4], pi_one[2,1]*pi_two[4,2]*pi_one[2,4]) +
             (q^3)*prod(pi_one[2,1], pi_one[2,2]*pi_one[2,3]*pi_one[2,4]))
  return(c(b_00, b_01, b_10, b_11))
}
########

########
# define functions for coefficients in cubic expression

## A1,2,3,4 for phi(0,0); group 1
A1_four <- function(n, a){
  A1 <- sum(n)*prod(a)
  return(A1)
}

A2_four <-function(n, a, b){
  t1 <- sum(n[1], n[2], n[3])*prod(a[1], a[2], a[3], b[4])
  t2 <- sum(n[1], n[2], n[4])*prod(a[1], a[2], b[3], a[4])
  t3 <- sum(n[1], n[3], n[4])*prod(a[1], b[2], a[3], a[4])
  t4 <- sum(n[2], n[3], n[4])*prod(b[1], a[2], a[3], a[4])
  return(sum(t1, t2, t3, t4))
}

A3_four <- function(n, a, b){
  t1 <- sum(n[1], n[2])*prod(a[1], a[2], b[3], b[4])
  t2 <- sum(n[1], n[3])*prod(a[1], b[2], a[3], b[4])
  t3 <- sum(n[1], n[4])*prod(a[1], b[2], b[3], a[4])
  t4 <- sum(n[2], n[3])*prod(b[1], a[2], a[3], b[4])
  t5 <- sum(n[2], n[4])*prod(b[1], a[2], b[3], a[4])
  t6 <- sum(n[3], n[4])*prod(b[1], b[2], a[3], a[4])
  return(sum(t1, t2, t3, t4, t5, t6))
}

A4_four <- function(n, a, b){
  t1 <- prod(n[1], a[1], b[2], b[3], b[4])
  t2 <- prod(n[2], b[1], a[2], b[3], b[4])
  t3 <- prod(n[3], b[1], b[2], a[3], b[4])
  t4 <- prod(n[4], b[1], b[2], b[3], a[4])
  return(sum(t1, t2, t3, t4))
}
########

########
# define functions for finding roots of cubic expression

## FUNCTION FOR PHI(0,0)=PI_{S,...,S+3}
roots_fourwise_gr1 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  
  # find b constants
  b_gr1 <- b_s_gr1(q, pi_one, pi_two, pi_three)
  
  # find cubic coefficients
  A1_gr1 <- A1_four(n_vec, a)
  A2_gr1 <- A2_four(n_vec, a, b_gr1)
  A3_gr1 <- A3_four(n_vec, a, b_gr1)
  A4_gr1 <- A4_four(n_vec, a, b_gr1)
  
  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr1, A2_gr1, A3_gr1, A4_gr1)[c(1:3)]
  
  # select real roots 
  real_roots <- pick_real(all_roots)
  
  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }
  
  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)
  
  # return the final estimate of phi(0,0)
  return(round(estimate, 10))
}


## FUNCTION FOR PHI(0,1)##
roots_fourwise_gr2 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  
  # find b constants
  b_gr2 <- b_s_gr2(q, pi_one, pi_two, pi_three)
  
  # find cubic coefficients
  A1_gr2 <- A1_four(n_vec, a)
  A2_gr2 <- A2_four(n_vec, a, b_gr2)
  A3_gr2 <- A3_four(n_vec, a, b_gr2)
  A4_gr2 <- A4_four(n_vec, a, b_gr2)
  
  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr2, A2_gr2, A3_gr2, A4_gr2)[c(1:3)]
  
  # select real roots 
  real_roots <- pick_real(all_roots)
  
  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }
  
  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)
  
  # return the final estimate of phi(0,1)
  return(round(estimate, 10))
}

## FUNCTION FOR PHI(1,0)##
roots_fourwise_gr3 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  
  # find b constants
  b_gr3 <- b_s_gr3(q, pi_one, pi_two, pi_three)

  # find cubic coefficients
  A1_gr3 <- A1_four(n_vec, a)
  A2_gr3 <- A2_four(n_vec, a, b_gr3)
  A3_gr3 <- A3_four(n_vec, a, b_gr3)
  A4_gr3 <- A4_four(n_vec, a, b_gr3)
  
  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr3, A2_gr3, A3_gr3, A4_gr3)[c(1:3)]

  
  # select real roots 
  real_roots <- pick_real(all_roots)

  
  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }
  
  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)
  
  # return the final estimate of phi(1,0)
  return(round(estimate, 10))
}

## FUNCTION FOR PHI(1,1)
roots_fourwise_gr4 <- function(q, n_vec, pi_one, pi_two, pi_three, lower_b, upper_b){
  # find a coefficients
  a <- a_s_four(q)
  
  # find b constants
  b_gr4 <- b_s_gr4(q, pi_one, pi_two, pi_three)
  
  # find cubic coefficients
  A1_gr4 <- A1_four(n_vec, a)
  A2_gr4 <- A2_four(n_vec, a, b_gr4)
  A3_gr4 <- A3_four(n_vec, a, b_gr4)
  A4_gr4 <- A4_four(n_vec, a, b_gr4)
  
  # find all roots of the cubic expression (real and complex)
  all_roots <- cubic(A1_gr4, A2_gr4, A3_gr4, A4_gr4)[c(1:3)]
  
  # select real roots 
  real_roots <- pick_real(all_roots)

  # pick candidate estimate from real roots
  if (length(real_roots)==1){
    candidate <- real_roots
  }
  else{
    candidate <- sort(real_roots)[2]
  }
  
  # select final estimate based on bounds
  estimate <- candidate*(candidate >= lower_b & candidate <= upper_b) + lower_b*(candidate < lower_b) + upper_b*(candidate > upper_b)
  
  # return the final estimate of phi(1,1)
  return(round(estimate, 10))
}
########
