###########
# READ ME #
###########

Please note that to reduce the runtime parallel computing was used to compute stable tilted correlations. 
If you are running this code on a machine with less than eight cores go to line 444 in the setup file and run the following to detect the number of cores on your machine: 

detectCores(all.tests = FALSE, logical = TRUE)

On line 445 set the number of cores to one minus the number returned by detectCores. 
Note also that parallel code was executed via the ‘doParallel’ package on a Linux machine. 
Linux machines support fork system calls, hence the global environment is automatically available to each worker core when a parallel function is called. 
Fork system calls are not supported on Windows machines. 
To run this code on a Windows machine please copy lines 16-647 of the setup file into a new R scrip and make a new package called ‘functionsIneed’. 
Wherever a ‘foreach’ function is used please add the following argument: 

.packages = c(‘functionsIneed’,…)

where '…' stands for the packages loaded in lines 6-14 of the setup file. 


