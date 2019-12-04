# CRAbun
This package aims to show how to implement empirical likelihood (EL) methods under the Huggin-Alho model in capture-recapture studies.
As an example, we analyze a real data to show how to reproduce the results of Table 3 in Liu et al. (2020+).

## Vignette about CRAbun package
### 1. A data set named "prinia"
Collected in Hongkong during 17 weeks from Junuary to April in 1993,
the yellow-bellied prinia data set consists of 163 observations and 6 columns：
id, number.of.capture, tail.length, fat.index, wing, and wing.index.

Here, the tail.length is a continuous variable, so is the attribute wing (short for wing length).
While, the fat.index and wing.index are both categorical variables.
The tail.length vaiable has 41 missing values, and other variables are always observed.

This dataset is used in Liu et al. (2020+) to illustrate the performance of EL method in the presence of missing data.

### 2. Optimization algorithm for EL methods
A function **abun.opt** is provided to obtain the maximum EL estimator and EL ratio confidence interval of the abundance. Std. Error is also given for reference.

### 3. Other functions
To compare the performance of the EL methods and the existing 
inverse probability weighting (IPW) and multiple impution (MI) methods
proposed by Lee et al. (2016),
we give two functions:
- **ipw.mar**: used to implement the IPW method.
- **mi2.mar**: used to implement the MI method.

## Comment
- **abun.opt** function can deal with all cases whenever there is missing data or not.
- **mi2.mar** function is specific to the case where the missing variable *y* is univariate.

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
