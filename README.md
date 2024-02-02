## Differential Transcript Expression (DTE) benchmarking in R 

The main DTE benchmark script is [dte_run_benchmarks.Rmd](https://github.com/gpertea/r-dte/blob/master/dte_run_benchmarks.Rmd)

The prepackaged simulation data RSE is expected in the `sdownstream` folder which can be downloaded from here:
https://drive.google.com/drive/folders/1s6PhdcdUmDP_jiZy5_uNNyAZVAF-roDh?usp=sharing

## Adding custom DTE methods to the benchmark

Example: in order to add `limmaTrend` on the TPR-FDR plot the following actions were needed:
 - add the custom DTE function `limmtrend()` inf [dte_funcs_rse.R](https://github.com/gpertea/r-dte/blob/master/dte_funcs_rse.R#L84), taking the RSE object as parameter
 - add a corresponding color entry in [dte_plots.R](https://github.com/gpertea/r-dte/blob/master/dte_plots.R#L9) `cols`
 - add the code to call the new function and add the results as an element to `padj` list for iCOBRA in [dte_run_benchmarks.Rmd](https://github.com/gpertea/r-dte/blob/master/dte_run_benchmarks.Rmd#L61)

