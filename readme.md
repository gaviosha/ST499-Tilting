# READ ME 

Code in this repository reproduces tables and figures in my MSc dissertation titled "Screening High Dimensional Time Series via Tilting".

## Directories

The directory structure for this repository is as follows: 

```
ST499-Tilting-
    |
    |--dissertation 
    |
    |--setup
    |   |-- Setup.R
    |
    |-- data
    |   |-- clean...csv
    |
    |-- section-3
    |   |-- Figure...R
    |   |-- Figure...R
    |
    |-- section-4
    |   |-- Figure...R
    |   |-- Figure...R
```

The function of each sub-directory is as follows: 

* `dissertation`: contains a pdf of the dissertation as submitted for the MSc
* `setup`: a .R file used to set up the R environment needed to produce figures and tables in the dissertation
* `data`: data produced according to section 4.2.1
* `section-3`: scripts used to produce results in section 3 
* `section-4`: scripts used to produce results in section 4

## Depends 

R >= 3.0.0

## Running parallel code 

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
