# CRAbun
This package aims to show how to implement empirical likelihood (EL) methods under the Huggin-Alho model in capture-recapture studies.
As an example, we analyze a real data to show how to reproduce the results of Table 3 in Liu et al. (2020+).

## Vignette about CRAbun package
### 1. A data set named "prinia"
Collected in Hongkong during 17 weeks from Junuary to April in 1993,
the yellow-bellied prinia data set consists of 163 observations and 6 columns：
id, number.of.capture, tail.length, fat.index, wing, and wing.index.

Here, the tail.length vaiable has 41 missing values, and other variables are always observed. We expected to illustrate the performance of EL method in the presence of missing data in Liu et al. (2020+).

### 2. Three main functions
- **abun.opt**: used to implement the EL method whenever there is missing data or not. See  Liu et al. (2017) and Liu et al. (2020+).
- **ipw.mar**: used to implement the inverse probability weighting method in the presence of missing data. See Lee et al. (2016).
- **mi2.mar**: used to implement the multiple impution (MI) method in the presence of missing data. See Lee et al. (2016).

### 3. Other functions
## Comment
- **mi2.mar** function is specific to the case where the missing variable is univariate.
- The results in Table 3 in Liu et al. (2020+) can be reproduced by typing the following lines to the R software:


*library(devtools)*

*install_github('ecnuliuyang/CRAbun')*

*library(CRAbun)*

*example("ipw.mar")*

*example("mi2.mar")*

*example("abun.opt")*



## References
Lee, S.-M., W.-H. Hwang, and J. de Dieu Tapsoba (2016). 
Estimation in closed capture–recapture models when covariates are missing at random. 
*Biometrics* **72**(4), 1294--1304.

Liu, Y., P. Li, and J. Qin (2017). 
Maximum empirical likelihood estimation for abundance in a closed
population from capture-recapture data. 
*Biometrika* **104**(3), 527--543.

Liu, Y., Y. Liu, P. Li, and L. Zhu (2020+).
Maximum likelihood abundance estimation from capture-recapture data when covariates are missing at random.
Submitted.

## Contact
Yang Liu, liuyangecnu@163.com
