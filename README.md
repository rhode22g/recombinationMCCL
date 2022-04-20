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

This function returns a matrix of the pairwise marginal estimates. It takes four parameters: (1) `L` is the length of the SNP sequences in the sample, (2) `q` is the chosen probability of recombination, (3) `n` is the number of descendant in the sample, and (4) `d` is a data frame of the sample of descendant sequences where rows index descendants and columns index SNP sites (number of rows must equal `n` and number of columns must equal `L`)

`estimates_m2`

`estimates_m3`

## Subfunctions

# JOINT ESTIMATION FUNCTIONS

`ancestor_pairs`

`ancestor_three`

`ancestor_four`

# SIMULATING DATA FUNCTIONS

## Main Functions

## Subfunctions
