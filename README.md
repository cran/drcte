# drcte

Overview
--------

This package contains several service functions to be used in support to the package 'drc', to analyse datasets obtained from time-to-event and other types of censored data in agriculture

Installation
------------

``` r
# Last stable version (From CRAN)
# install.packages("drcte")

# You can also install the development version of 'drcte' from GitHub
# install.packages("devtools")
# devtools::install_github("onofriAndreaPG/drcte")
```

Known bugs
----------

The use of the 'units' argument in the 'summary' method gives a warning when the pmodels argument is used in the 'drmte' function
