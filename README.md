# Identifying the invalid IVs in Mendelian Randomization by using Lasso

`lasso_iv` is designed to detect the invalid IVs in MR using  Lasso. This function requires individual-level data and incorporates the selection process via a bootstrap method. The stopping rule for the Lasso procedure involves comparing the Cochran'Q statistics and the critical value of a Chi-squared distribution with a degree of freedom equal to the number of  valid IVs minus 1. `lasso_iv` is built on `mr_lasso` [1] in `MendelianRandomization` Package. 

---
### Reference:
1. [Rees JM, Wood AM, Dudbridge F, Burgess S. Robust methods in Mendelian randomization via penalization of heterogeneous causal estimates. _PloS one_. 2019 Sep 23;14(9):e0222362.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0222362)
