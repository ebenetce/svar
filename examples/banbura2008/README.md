# Banbura, Giannone and Reichlin (2008)

This folder contains a MATLAB replication driver for `RePEc_eca_wpaper_2008_033.pdf`, "Large Bayesian VARs".

The paper uses the 131-series monthly Stock and Watson (2005a) macroeconomic dataset from January 1959 through December 2003. The data are not committed here. To run the replication, place `banbura2008Data.mat` in this folder with:

- `tt`: a monthly timetable whose variable names include the paper names used in `banbura2008ModelSpecifications`. For the LARGE model, keep the columns in the paper's structural order: slow-moving variables, `FFR`, then fast-moving variables.
- optional `priorMean`: a 1-by-N vector of Minnesota prior means, with `1` for random-walk variables and `0` for stationary variables.

The official JAE replication archive is available at the Queen's Economics `jae` archive for Banbura, Giannone and Reichlin. If `bgrdata.xls` is available locally, run:

```matlab
prepareBanbura2008Data("archive/bgr-xls/bgrdata.xls", "banbura2008Data.mat")
```

`main.m` reads the PDF with Text Analytics Toolbox, estimates the paper's BVAR panels with the existing `minnesotabvarm` prior, evaluates rolling forecast MSFEs, and computes the monetary-policy IRFs and FEVDs.

For a quick smoke run on the real archive data, set `BGR_FAST=1` before running `main.m`. Fast mode uses fewer horizons, fewer posterior draws, and only the first 12 forecast origins, so its MSFE table is a wiring check rather than the paper's reported forecast table.
