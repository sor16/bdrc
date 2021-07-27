## Test environments
* local OS X install, R 4.1.0
* macOS-latest (4.1) on Github Actions
* windows-latest (4.1) on Github Actions
* ubuntu 18.04 (devel,4.1,4.0,3.6,3.5) on Github Actions

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## CRAN Policy issues
* One helper function, save_report() writes a pdf file to the user’s file system, but it requires permission from the user in an interactive R session
* The examples are slow, because in order to demonstrate the usage of most of the functions, a model has to be run. The total running time for R CMD CHECK where the examples are run takes around 3-4 hours. Therefore, after testing the examples extensively, I decided to wrap them in \dontrun{} to cut the checking time drastically.


