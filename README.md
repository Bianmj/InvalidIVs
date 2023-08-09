# Identifying the invalid IVs in Mendelian Randomization by using Lasso

`lasso_iv` is designed to detect the invalid IVs in MR using  Lasso. This function requires individual-level data and incorporates the selection process via a bootstrap method. The stopping rule for the Lasso procedure involves comparing the Cochran'Q statistics and the critical value of a Chi-squared distribution with a degree of freedom equal to the number of  valid IVs minus 1. `lasso_iv` is built on `mr_lasso` in `MendelianRandomization` Package. 

---
