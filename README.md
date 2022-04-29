# recombinationMCCL
## this is a preliminary version of the package which can carry out the marginal estimation and joint distribution reconstruction for m=1,2,3


# BUILT IN DATA

This package includes one built-in dataset, both the raw data and a version of the data converted into binary sequences. The raw data comes from the International HapMap Project (https://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2009-02_phaseIII/HapMap3_r2/).

The raw data is available under `yri_trio_1.csv`. This corresponds to data on Chromosome 1 from the YRI population Trios data.

The converted data is available under `yri_trio_1_binary.csv`.

# MARGINAL ESTIMATION FUNCTIONS

## Main Functions

`estimates_m0`: 

This function returns a matrix of the onewise marginal estimates. It takes three parameters: (1) `descendents` is a data frame of the sample of descendant sequences where rows index descendants and columns index SNP sites, (2) `L` is the length of the SNP sequences which must equal the number of columns in `descendents`, and (3) `n` is the number of descendants in the sample which must equal the number of rows in `descendents`. In the output, rows correspond to the parameters in order of bitwise value and the columns correspond to sites.

`estimates_m1`

This function returns a matrix of the pairwise marginal estimates. It takes four parameters: (1) `L` is the length of the SNP sequences in the sample, (2) `q` is the chosen probability of recombination, (3) `n` is the number of descendants in the sample, and (4) `d` is a data frame of the sample of descendant sequences where rows index descendants and columns index SNP sites (number of rows must equal `n` and number of columns must equal `L`). In the output, rows correspond to the parameters in order of bitwise value and the columns correspond to sites.

`estimates_m2`

This function returns a matrix of threewise marginal estimates. It takes four parameters: (1) `L` is the length of the SNP sequences in the sample, (2) `q` is the chosen probability of recombination, (3) `n` is the number of descendants in the sample, and (4) `descend` is a data frame of the sample of descendant sequences where rows index descendants and columns index SNP sites (number of rows must equal `n` and number of columns must equal `L`). In the output, rows correspond to the parameters in order of bitwise value and the columns correspond to sites.

`estimates_m3`

This function returns a matrix of fourwise marginal estimates. It takes four parameters: (1) `L` is the length of the SNP sequences in the sample, (2) `q` is the chosen probability of recombination, (3) `n` is the number of descendants in the sample, and (4) `d` is a data frame of the sample of descendant sequences where rows index descendants and columns index SNP sites (number of rows must equal `n` and number of columns must equal `L`). In the output, rows correspond to the parameters in order of bitwise value and the columns correspond to sites.

## Subfunctions

### m=2 subfunctions

`a_s`

This function calculates the $a_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function. It has one parameter `q`, the recombination probability. It returns a vector of the four coefficients.

`b_s0`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(0)$. It takes two parameters: (1) `q` is the recombination probability, (2) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, and (3) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function. It returns a vector of the four coefficients.

`b_s1`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(1)$. It takes two parameters: (1) `q` is the recombination probability, (2) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, and (3) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function. It returns a vector of the four coefficients.

`A_1`

This function calculates the $$A_1^3$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes two parameters: (1) `n` is the number of descendants in the sample, and (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s`). It returns the numerical value of the coefficient.

`A_2`

This function calculates the $$A_2^3$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes three parameters: (1) `n` is the number of descendants in the sample, (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s`), and (3) `b` is a vector of the $$b_{x_s,x_{s+3}}^s$$ coefficients (the return of `b_s0` or `b_s1`). It returns the numerical value of the coefficient.

`A_3`

This function calculates the $$A_3^3$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes three parameters: (1) `n` is the number of descendants in the sample, (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s`), and (3) `b` is a vector of the $$b_{x_s,x_{s+3}}^s$$ coefficients (the return of `b_s0` or `b_s1`). It returns the numerical value of the coefficient.

`A_4`

This function calculates the $$A_4^3$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes three parameters: (1) `n` is the number of descendants in the sample, (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s`), and (3) `b` is a vector of the $$b_{x_s,x_{s+3}}^s$$ coefficients (the return of `b_s0` or `b_s1`). It returns the numerical value of the coefficient.

`pick_real`

This function finds the real values in a vector of values. This is used to pick the roots of the cubic functions that are real roots. This is done by finding the values which have an imaginary component greater than `1e-4`, then returning the elements in the vector which don't meet that threshold. 

`roots_threewise_0`

This function takes six parameters: (1) `q` is the recombination probability, (2) `n` is a vector of the number of descendants in the sample with the sequences in group 1, (3) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (4) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, (5) `lower_b` is the lower bound of the interval we constrain $$\hat{\phi}_s(0)$$ to, and (6) `upper_b` is the upper bound of the interval we constrain $$\hat{\phi}_s(0)$$ to. This function includes calls to the subfunctions which calculate coefficients. It returns the constrained, final estimate of $$\hat{\phi}_s(0)$$.

 
`roots_threewise_1`

This function takes six parameters: (1) `q` is the recombination probability, (2) `n` is a vector of the number of descendants in the sample with the sequences in group 2, (3) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (4) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, (5) `lower_b` is the lower bound of the interval we constrain $$\hat{\phi}_s(1)$$ to, and (6) `upper_b` is the upper bound of the interval we constrain $$\hat{\phi}_s(1)$$ to. This function includes calls to the subfunctions which calculate coefficients. It returns the constrained, final estimate of $$\hat{\phi}_s(1)$$.

### m=3 subfunctions

`a_s_four`

This function calculates the $a_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the fourwise probability function. It has one parameter `q`, the recombination probability. It returns a vector of the four coefficients.

`b_s_gr1`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(0,0)$. It takes four parameters: (1) `q` is the recombination probability, (2) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (3) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, and (4) `pi_three` is a subset of the matrix of threewise marginal estimates (output of `estimates_m2`) taken iteratively in the main function. It returns a vector of the four coefficients.

`b_s_gr2`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(0,1)$. It takes four parameters: (1) `q` is the recombination probability, (2) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (3) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, and (4) `pi_three` is a subset of the matrix of threewise marginal estimates (output of `estimates_m2`) taken iteratively in the main function. It returns a vector of the four coefficients.

`b_s_gr3`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(1,0)$. It takes four parameters: (1) `q` is the recombination probability, (2) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (3) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, and (4) `pi_three` is a subset of the matrix of threewise marginal estimates (output of `estimates_m2`) taken iteratively in the main function. It returns a vector of the four coefficients.

`b_s_gr4`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(1,1)$. It takes four parameters: (1) `q` is the recombination probability, (2) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (3) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, and (4) `pi_three` is a subset of the matrix of threewise marginal estimates (output of `estimates_m2`) taken iteratively in the main function. It returns a vector of the four coefficients.

`A1_four`

This function calculates the $$A_1^4$$ coefficient for rewriting the score equation of the fourwise log-MCCL. It takes two parameters: (1) `n` is the number of descendants in the sample, and (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s_four`). It returns the numerical value of the coefficient.

`A2_four`

This function calculates the $$A_2^4$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes three parameters: (1) `n` is the number of descendants in the sample, (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s_four`), and (3) `b` is a vector of the $$b_{x_s,x_{s+3}}^s$$ coefficients (the return of `b_s_gr1`, `b_s_gr2`, `b_s_gr3`, or `b_s_gr4`). It returns the numerical value of the coefficient.

`A3_four`

This function calculates the $$A_3^4$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes three parameters: (1) `n` is the number of descendants in the sample, (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s_four`), and (3) `b` is a vector of the $$b_{x_s,x_{s+3}}^s$$ coefficients (the return of `b_s_gr1`, `b_s_gr2`, `b_s_gr3`, or `b_s_gr4`). It returns the numerical value of the coefficient.

`A4_four`

This function calculates the $$A_4^4$$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes three parameters: (1) `n` is the number of descendants in the sample, (2) `a` is a vector of the $$a_{x_s,x_{s+3}}^s$$ coefficients (the return of `a_s_four`), and (3) `b` is a vector of the $$b_{x_s,x_{s+3}}^s$$ coefficients (the return of `b_s_gr1`, `b_s_gr2`, `b_s_gr3`, or `b_s_gr4`). It returns the numerical value of the coefficient.

`roots_fourwise_gr1`

This function takes six parameters: (1) `q` is the recombination probability, (2) `n_vec` is a vector of the number of descendants in the sample with the sequences in group 1 calculated by the main function, (3) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (4) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, (5) `lower_b` is the lower bound of the interval we constrain $$\hat{\phi}_s(0,0)$$ to, and (6) `upper_b` is the upper bound of the interval we constrain $$\hat{\phi}_s(0,0)$$ to. This function includes calls to the subfunctions which calculate coefficients. It returns the constrained, final estimate of $$\hat{\phi}_s(0,0)$$.

`roots_fourwise_gr2`

This function takes six parameters: (1) `q` is the recombination probability, (2) `n_vec` is a vector of the number of descendants in the sample with the sequences in group 2 calculated by the main function, (3) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (4) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, (5) `lower_b` is the lower bound of the interval we constrain $$\hat{\phi}_s(0,1)$$ to, and (6) `upper_b` is the upper bound of the interval we constrain $$\hat{\phi}_s(0,1)$$ to. This function includes calls to the subfunctions which calculate coefficients. It returns the constrained, final estimate of $$\hat{\phi}_s(0,1)$$.

`roots_fourwise_gr3`

This function takes six parameters: (1) `q` is the recombination probability, (2) `n_vec` is a vector of the number of descendants in the sample with the sequences in group 3 calculated by the main function, (3) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (4) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, (5) `lower_b` is the lower bound of the interval we constrain $$\hat{\phi}_s(1,0)$$ to, and (6) `upper_b` is the upper bound of the interval we constrain $$\hat{\phi}_s(1,0)$$ to. This function includes calls to the subfunctions which calculate coefficients. It returns the constrained, final estimate of $$\hat{\phi}_s(1,0)$$.

`roots_fourwise_gr4`

This function takes six parameters: (1) `q` is the recombination probability, (2) `n_vec` is a vector of the number of descendants in the sample with the sequences in group 4 calculated by the main function, (3) `pi_one` is a subset of the matrix of onewise marginal estimates (output of `estimates_m0`) taken iteratively in the main function, (4) `pi_two` is a subset of the matrix of pairwise marginal estimates (output of `estimates_m1`) taken iteratively in the main function, (5) `lower_b` is the lower bound of the interval we constrain $$\hat{\phi}_s(1,1)$$ to, and (6) `upper_b` is the upper bound of the interval we constrain $$\hat{\phi}_s(1,1)$$ to. This function includes calls to the subfunctions which calculate coefficients. It returns the constrained, final estimate of $$\hat{\phi}_s(1,1)$$.

# JOINT ESTIMATION FUNCTIONS

`ancestor_pairs`

`ancestor_three`

`ancestor_four`

# SIMULATING DATA FUNCTIONS

## Main Functions

## Subfunctions
