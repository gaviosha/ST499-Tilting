# Screening High Dimensional Time Series via Tilting

This repository contains code for building my MSc dissertation, which investigated a new feature screening procedure for high dimensional time series based on [Cho & Fryzlewicz (2012)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2011.01023.x).

## Abstract

Developing effective feature screening procedures when the number of features exceeds the sample size (p > n) is one of the most active research areas in modern statistics. Although much attention has been devoted to the problem of detecting causal and predictive relationships in high dimensional regimes, comparatively little attention has been given to the problem of developing general purpose screening procedures for high dimensional time series data. This is especially surprising given high dimensional time series are becoming increasingly common in many fields including economics, finance, and neuro-science.

The aim of this dissertation is to adapt tilted correlation screening, a variable selection procedure for linear models introduced by [Cho & Fryzlewicz (2012)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2011.01023.x), to settings in which the response, predictor, and error variables are allowed to be time series processes. Overall, the dissertation aims to make the following original contributions:

1. Show that under mild assumptions the main theoretical results in [Cho & Fryzlewicz (2012)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2011.01023.x) hold for time series data.
2. Introduce a generalised least squares variant of the tilted correlation which is more efficient than the original (in the Gauss-Markov / MSE sense).
3. Apply the principle of stability selection introduced by [Meinshausen & Bühlmann (2010)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2010.00740.x) to tilted correlations when the correlation structure of the predictor variables necessitates careful control.

The theoretical results presented are reinforced by simulation studies. Additionally, a real data application is provided in which tilted correlation screening is used to select a model for forecasting CPI inflation, chained GDP, and the Sterling effective exchange rate.

## Usage

### Building the dissertation

To generate the dissertation itself build `main.tex` from within `.\tex`.

### Generating plots and figures

To set up the global environment run the `setup.R` file in `.\R`, then see R scripts for each plot.

### Running parallel code

Please note that to reduce the runtime parallel computing was used to compute stable tilted correlations.
If you are running this code on a machine with less than eight cores go to line **444** in the setup file and run the following to detect the number of cores on your machine:

```
detectCores(all.tests = FALSE, logical = TRUE)
```

On line **445** set the number of cores to one minus the number returned by `detectCores`.

Note also that parallel code was executed via the `doParallel` package on a Linux machine. Linux machines support fork system calls in which the global environment is automatically available to each worker core when a parallel function is called. Fork system calls are not supported on Windows machines.

To run this code on a Windows machine please copy lines **16-647** of the setup file into a new R scrip and make a new package called `functionsIneed`.

Wherever a `foreach` function is used please add the following argument:

```
.packages = c(‘functionsIneed’,…)
```

where `…` stands for the packages loaded in lines **6-14** of the setup file.
