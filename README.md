# CRAbun
Some new methods are developed to make statistical inference for the abundance under capture-recapture models.
Below are key functions in this package.


### CRMiss function
When covariates are missing at random, four estimation methods for the abundance are carried out 
under capture-recapture models.

First method is a new maximum empirical likelihood estimation method (MELE).
Refer to Liu et al. (2019+).

Second method is based on the naive maximum empirical likelihood derived from the complete cases (CC).
See Liu et al. (2017).

The last two are based on the inverse probability weighting (IPW) and multiple imputation (MI) methods.
See Lee et al. (2016).

# References

Lee, S.-M., W.-H. Hwang, and J. de Dieu Tapsoba (2016). 
Estimation in closed capture–recapture models when covariates are missing at random. 
\emph{Biometrics} \strong{72}(4), 1294--1304.

Liu, Y., P. Li, and J. Qin (2017). 
Maximum empirical likelihood estimation for abundance in a closed
population from capture-recapture data. 
\emph{Biometrika} \strong{104}(3), 527--543.

Liu, Y., Y. Liu, P. Li, and L. Zhu (2019+).
Maximum likelihood abundance estimation from capture-recapture data when covariates are missing at random.
Submitted.

# Contact

Yang Liu, liuyangecnu@163.com
