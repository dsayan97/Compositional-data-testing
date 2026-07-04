# EBC: Empirical Bootstrap Test for High-Dimensional Compositional Data

R implementation of the **Empirical Bootstrap for Compositional data (EBC)** test proposed in

> **Testing Microbiome Community Differences in High Dimensions: A Bootstrap Approach for Compositional Data**

The proposed method provides a nonparametric bootstrap test for comparing high-dimensional compositional mean vectors after centered log-ratio (CLR) transformation. The procedure is applicable to two-sample, paired-sample, and multi-sample hypothesis testing and is particularly designed for microbiome studies where the number of taxa greatly exceeds the sample size.

---

## Repository Contents

| File | Description |
|------|-------------|
| `EBC.R` | Implementation of the EBC hypothesis test. Accepts a list of CLR-transformed data matrices and returns the bootstrap p-value. |
| `Run_two.R` | Simulation code for reproducing the **two-sample** simulation results reported in the paper. Includes data generation, zero inflation, zero replacement, CLR transformation, and hypothesis testing. |
| `Run_K.R` | Simulation code for reproducing the **multi-sample** simulation results reported in the paper. Produces results for 3-, 4-, and 5-sample comparisons. |

---

## Requirements

The simulations require R (version 4.0 or later) together with the following packages:

```r
install.packages(c("MASS", "zCompositions"))
```

---

## Running the EBC Test

First source the implementation

```r
source("EBC.R")
```

Suppose the CLR-transformed observations from each group are stored as matrices

```r
X1
X2
X3
```

where each matrix has

- rows = observations
- columns = taxa (features)

Run the test as

```r
pvalue <- EBC(list(X1, X2))
```

or for multiple groups

```r
pvalue <- EBC(list(X1, X2, X3))
```

The function returns the empirical bootstrap p-value.

---

## Input Format

The function

```r
EBC(X, N = 1000)
```

expects

- `X` : list of CLR-transformed matrices, one matrix per population.
- `N` : number of bootstrap resamples (default 1000).

Each matrix should have dimensions

```
n_k × p
```

where

- `n_k` = sample size of group *k*
- `p` = number of taxa/features.

---

## Reproducing the Simulation Results

### Two-sample simulations

Open `Run_two.R` and modify the user-controlled parameters near the beginning of the script.

Typical parameters include

- sample sizes
- dimension (`p`)
- signal strength (`d`)
- covariance structure
- distribution (Normal or t)
- zero inflation level
- number of Monte Carlo iterations.

Then simply run

```r
source("Run_two.R")
```

The script generates one simulation configuration corresponding to a single table entry in the manuscript.

---

### Multi-sample simulations

Similarly,

```r
source("Run_K.R")
```

produces simulation results for the multi-sample experiments.

The script reports empirical rejection probabilities for

- 3 populations
- 4 populations
- 5 populations

under the selected simulation setting.

---

## Simulation Pipeline

Both simulation scripts follow the same workflow:

1. Generate latent log-abundances.
2. Generate sequencing depths from negative binomial distributions.
3. Generate multinomial count data.
4. Introduce additional zero inflation.
5. Perform zero replacement using **zCompositions**.
6. Apply the centered log-ratio (CLR) transformation.
7. Compute the EBC test statistic.
8. Estimate the bootstrap p-value.

---

## Citation

If you use this code, please cite

> Bhattacharjee, M., Chakraborty, N., Das, S., Chakraborty, S., Liu, L., Shi, Y., Wylie, K. M., Wylie, T. N., and Stout, M. J., *Testing Microbiome Community Differences in High Dimensions: A Bootstrap Approach for Compositional Data.*

---

## Contact

For questions regarding the methodology or implementation, please contact the corresponding author listed in the manuscript.
