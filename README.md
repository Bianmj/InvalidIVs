# Identifying the invalid IVs in Mendelian Randomization with the Lasso procedure

`lasso_iv` is a function for the identification of invalid IVs in MR using the Lasso procedure and it needs the individual-level data. The selection procedure is also taken into account by using the bootstrap procedure. The stopping rule for the Lasso procedure involves comparing the Cochran'Q statistics and the critical value of a Chi-squared distribution with a degree of freedom equal to the number of  valid IVs minus 1. `lasso_iv` is built on `mr_lasso` in `MendelianRandomization` Package. 

---
