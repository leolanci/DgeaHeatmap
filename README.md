
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DgeaHeatmap

<!-- badges: start -->
<!-- badges: end -->

The goal of DgeaHeatmap is to enable R users to generate heatmaps more
easily and to help with preprocessing read counts. Furthermore, the
package is aimed to simplify the extraction of raw read counts from .dcc
and .pkc files generated through Nanostring GeoMx DSP.

## Installation

You can install the development version of DgeaHeatmap from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("leolanci/Dgea_Heatmap_Package")
```

## Usage

This is a basic example which shows you how to solve a common problem:

``` r
library(DgeaHeatmap)
x <- 1
matrixCounts <- build_matrix(input_data, x)
```

|               | P0_cortex_Iba1_pos_1 | P0_cortex_Iba1_pos_2 | P0_cortex_Iba1_pos_3 | P0_cortex_Iba1_neg_1 |
|:--------------|---------------------:|---------------------:|---------------------:|---------------------:|
| Casp6         |                 10.5 |                 7.00 |                  2.8 |                 8.10 |
| Atl3          |                  7.0 |                13.99 |                  8.4 |                12.23 |
| C030006K11Rik |                 10.5 |                 7.00 |                  8.4 |                10.16 |
| Cflar         |                 10.5 |                17.49 |                  5.6 |                 7.22 |
| Aftph         |                  7.0 |                10.50 |                  8.4 |                11.78 |
| Tmem41b       |                  3.5 |                 7.00 |                  8.4 |                 6.63 |

``` r
parameter1 = "cortex"
parameter2 = "pos"
factors_for_individual_matrix = list(parameter1, parameter2)
indiMatrix <- individual_matrix(factors_for_individual_matrix, matrixCounts)
```

|               | P0_cortex_Iba1_pos_1 | P0_cortex_Iba1_pos_2 | P0_cortex_Iba1_pos_3 | P0_cortex_layer_6_Iba1_pos_2 | P0_cortex_layer_6_Iba1_pos_3 | P5_cortex_Iba1_pos_1 | P5_cortex_Iba1_pos_2 | P5_cortex_Iba1_pos_3 | P5_cortex_layer_6_Iba1_pos_1 |
|:--------------|---------------------:|---------------------:|---------------------:|-----------------------------:|-----------------------------:|---------------------:|---------------------:|---------------------:|-----------------------------:|
| Casp6         |                 10.5 |                 7.00 |                  2.8 |                         4.66 |                         8.40 |                11.20 |                 4.66 |                16.79 |                        13.99 |
| Atl3          |                  7.0 |                13.99 |                  8.4 |                         7.00 |                         8.40 |                13.99 |                 9.33 |                 5.60 |                         7.00 |
| C030006K11Rik |                 10.5 |                 7.00 |                  8.4 |                         4.66 |                         2.80 |                13.99 |                 4.66 |                11.20 |                         3.50 |
| Cflar         |                 10.5 |                17.49 |                  5.6 |                        16.33 |                        13.99 |                13.99 |                18.66 |                 8.40 |                         3.50 |
| Aftph         |                  7.0 |                10.50 |                  8.4 |                         2.33 |                         8.40 |                 2.80 |                 4.66 |                 5.60 |                         7.00 |
| Tmem41b       |                  3.5 |                 7.00 |                  8.4 |                         4.66 |                         8.40 |                 5.60 |                 4.66 |                 5.60 |                        20.99 |

Filtering the matrix for only an x amount of most variable expressed
genes:

``` r
top_number_of_genes <- 500
varGenesMatrix <- filtering_for_top_exprGenes(indiMatrix, top_number_of_genes)
print(nrow(varGenesMatrix))
#> [1] 500
```

Using Z-score scaling to scale the values of the matrix:

``` r
scaled_counts <- scale_counts(varGenesMatrix)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
