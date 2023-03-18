## Package update
New version number: 1.1.0

## Test environments
* local OS X install, R 4.2.0
* macOS-latest (release) on Github Actions
* windows-latest (release) on Github Actions
* ubuntu-latest (devel,release,old-rel) on Github Actions

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## CRAN Policy issues
* This update is a slight update of the package but also a reaction to a warning in CMD CHECK appearing in the development version of R regarding consistency of generics in the package. This issue has been fixed.
* One helper function, save_report() writes a pdf file to the user’s file system, but it requires permission from the user in an interactive R session
* The examples are slow, because in order to demonstrate the usage of most of the functions, a model has to be run. Therefore, I wrapped them in \donttest{}. Checking time without the examples is well under 10 minutes, around 5 minutes on my machine.


