# CRAbun
This package aims to show how to implement empirical likelihood (EL) methods under the Huggin-Alho model in capture-recapture studies.
As an example, we analyze a real data to show how to reproduce the results of Table 3 in Liu et al. (2020+).

## Vignette about CRAbun package
### 1. A Data set named "prinia"
Collected in Hongkong during 17 weeks from Junuary to April in 1993,
Yellow-bellied prinia data set consists of 163 observations and 6 columns：
id, number.of.capture, tail.length, fat.index, wing, and wing.index.

Here, the tail.length is a continuous variable, so is the wing column (wing length).
While, the fat.index and wing.index are both categorical variables.

### 2. Optimization algorithm for EL methods
To obtain the maximum EL estimator and EL ratio confidence interval of the abundance
we provide a function **abun.opt** with 8 options.
- d: vector stands for the number of times of being captured.
- K: number stands for the number of capture occasions.
- x: vector, matrix or data.frame stands for individual covariates without missingness.
- y: vector, matrix or data.frame stands for individual covariates with missingness.
- beta.initial: vector has the same length wih the coefficients in Huggins-Alho model. 0 means that the corresponding coefficient is fixed to zero, and non-zero values are seen as initial values when optimizing.
- CI: logistic with TRUE(default) to show the EL ratio confidence interval of *N*
and FALSE to not.
- level: number stands for the confidence level with default of 0.95.
- SE: logistic with TRUE(default) to show the standard error
and FALSE to not.


### 3. Other functions
To compare the performance of the EL methods and the existing 
inverse probability weighting (IPW) and multiple impution (MI) methods
proposed by Lee et al. (2016),
we give two functions:
- **ipw.mar**: used to implement the IPW method.
- **mi2.mar**: used to implement the MI method.

## Comment
**mi2.mar** function is specific to the case where the missing variable *y* is univariate.

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
