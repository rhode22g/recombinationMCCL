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

This function calculates the $a_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function. It has one parameter `q`, the recombination probability.

`b_s0`

This function calculates the $b_{x_s,x_{s+3}}^s$ coefficients for the process of rewriting the threewise probability function for $\phi_s(0)$. It takes two parameters: (1) `q` is the recombination probability, (2) `pi_one` is a matrix of onewise marginal estimates (the output of `estimates_m0`), and (3) `pi_two` is a matix of the pairwise marginal estimates (the output of `estimates_m1`). 

`A_1`

This function calculates the $A_1^3$ coefficient for rewriting the score equation of the threewise log-MCCL. It takes two parameters: (1) `n` is the number of descendants in the sample, and (2) `a` is a vector of the $a_{x_s,x_{s+3}}^s$ coefficients (the return of `a_s`).

`A_2`

`A_3`

`A_4`

`pick_real`

`roots_threewise_0`

`b_s1`

`roots_threewise_1`

### m=3 subfunctions

`a_s_four`

`b_s_gr1`

`b_s_gr2`

`b_s_gr3`

`b_s_gr4`

`A1_four`

`A2_four`

`A3_four`

`A4_four`

`roots_fourwise_gr1`

`roots_fourwise_gr2`

`roots_fourwise_gr3`

`roots_fourwise_gr4`

# JOINT ESTIMATION FUNCTIONS

`ancestor_pairs`

`ancestor_three`

`ancestor_four`

# SIMULATING DATA FUNCTIONS

## Main Functions

## Subfunctions
