######################
# SIMULATION FUNCTIONS
######################

#################
# ANCESTOR DIST'N
#################

### MAIN FUNCTION ###

ancestor <- function(L, hapmap_data){
  ###############
  # select a subset of the data for this sequence length
  hapmap_L <- hapmap_data[1:L,]
  # initialize empty row
  haplotypeL <- rep(NA, length = dim(hapmap_L)[2])
  # combine new row with dataset of first three SNPs
  hapmap_L <- rbind(hapmap_L, haplotypeL)
  
  # iterate over COLUMNS
  for (i in 1:dim(hapmap_L)[2]){
    hap <- vector(length = L)
    # iterate over ROWS
    for (j in 1:(dim(hapmap_L)[1]-1)){
      hap[j] <- hapmap_L[j,i]
      # print(hap)
      # print(unr_bin_L3[j,i])
    }
    full_hap <- toString(hap)
    hapmap_L[(L+1),i] <- full_hap}
  
  # vector to store just the strings
  pool_ancestor <- hapmap_L[(L+1),]
  ###############
  
  ###############
  # create a vector of the possible ancestor sequences in bitwise form
  possible_ancestor <- vector(length = (2^L))
  for (j in 0:((2^L)-1)){
    possible_ancestor[j+1] <- j
  }
  # convert pool of ancestor sequences to bitwise form
  pool_ancestor_decimal <- bitwise_vector(pool_ancestor, L)

  # generate simulated ancestor distribution
  
  ancestor_L <- sim_ancestor(possible_ancestor, pool_ancestor_decimal)
  ###############
  return(as.matrix(ancestor_L))
}


### SUB-FUNCTIONS ###

## write a function that calculates the decimal value based on general form (x_s * 2^(L-s))
### parameter x is the binary value at a site on the sequence, parameter e is the exponent (L-s)
bitwise <- function(x, e){
  decimal_val <- x * 2^e
  return(decimal_val)
}

## write a function to create a vector of the decimal forms of the 8 sequences
### parameter b is a a vector of sequences length L such that each item in the vector is of the form "x1, x2, x3, ..."
bitwise_vector <- function(b, l){
  # initialize an empty vector to hold the decimal forms
  for (s in 1:length(b)){
    # for each sequence initialize empty vector to hold the three terms
    bit_wise <- vector(length = l)
    # split sequence string
    sequence_split <- unlist(strsplit(b[s], ", "))
    # iterate over 3 terms in the sequence
    for (i in 1:l){
      x_s <- as.numeric(sequence_split[i])
      bit <- bitwise(x_s, (l-i))
      bit_wise[i] <- bit
    }
    decimal <- sum(bit_wise)
    sequence_bitwise[s] <- decimal 
  }
  return(sequence_bitwise)}

# write a function that calculates the sequence proportions / simulated ancestor distribution
## parameter possible is a vector of the decimal values for the possible ancestor sequences
## parameter pool is a vector of the decimal values for the ancestor pool generated from HapMap
sim_ancestor <- function(possible, pool){
  #initialize vector to hold counts
  sequencecounts <- vector(length = length(possible))
  # initialize all counts at 0
  for (v in 1:length(possible)){
    sequencecounts[v] <- 0
  }
  # iterate over each decimal in the pool
  for (x in 1:length(pool)){
    # for each decimal in the pool, check if it matches each possible and increment counts
    for (s in 1:length(possible)){
      if (pool[x] == possible[s]){
        sequencecounts[s] <- sequencecounts[s] + 1
      }
    }
  }
  # initialize empty vector to hold the proportions
  sequenceprop <- vector(length = length(possible))
  # calculate each proportion as the number of obs in pool with that decimal over the total number of observations
  for (g in 1:length(possible)){
    sequenceprop[g] <- sequencecounts[g] / length(pool)
  }
  # return the vector of proportions
  return(sequenceprop)
}


###################
# DESCENDANT SAMPLE 
###################

descendent_sample <- function (L, q, n, seed, hapmap){
  ######
  # REMEMBER TO SET SEED
  set.seed(seed)
  ######
  
  #######
  # utilize function ancestor defined elsewhere to get the ancestor dist'n for the entered length
  ancestor_pi <- ancestor(L, hapmap)
  colnames(ancestor_pi) <- "proportion"
  #######
  
  #######
  # initialize matrix to hold descendant sequences
  descendant_seq <- matrix(NA, nrow = n, ncol = L)
  #######
  
  #######
  # select whether recombination occurs at each of the possible locations; 0 no recombination, 1 recombination
  ## we need n*L samples because we need one for each descendant and then one for each site on the sequence for each descendant 
  ## probability is probability of recombination
  recombo_indx <- rbinom(n = L*n , size = 1, prob = q )
  head(recombo_indx)
  #######
  
  #######
  # select ancestor sequences
  ## sample from range 1:2^L, indices of rows in ancestor_pi
  ## number we are sampling same logic as above
  ancest_indx <- sample((1:(2^L)), size = L*n , prob = ancestor_pi, replace = TRUE)
  head(ancest_indx)
  #######
  
  #######
  # use a for loop to pick the n descendants
  ## iterate over the rows
  for(i in 1:n){
    # start with first column for each row
    start <- ((i-1)*L )
    # as starting point, set each descendant to initial ancestor sequence
    descendant_seq[i,] <- drop(digitsBase(x = ancest_indx[start + 1] - 1, base = 2, ndigits = L))
    
    # iterate over locations of possible recombinations
    for (k in 2:L){
      descendant_seq[i, k:L] <- (1 - recombo_indx[start + k])*descendant_seq[i, k:L] + recombo_indx[start + k]*drop(digitsBase(ancest_indx[k+start] - 1, base = 2, ndigits = L)[k:L])
    }
  }
  return(descendant_seq)
  #######
}
