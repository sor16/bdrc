# RCmodels

An R package called RCmodels was made for this project. The package includes the following functions written to run the models:

Both models
```clean()``` - Cleans the file input data. It can take in the two file formats shown in the test files. The txt format is the exported txt from WISKI database that the Icelandic Met Office uses. The xslx format is for other users that have data on other formats.
```priors()``` - Defines the prior parameters for a given country. If you want to add your country to the list that have specified prior parameters, please contact us.
Model 1

```plotmodel1()``` - makes use of clean() and model1BH() to plot a rating curve calculated with the methods in model 1.
```model1BH()``` - Makes use of Densevalm11() to determine the fit and confidence intervals of a rating curve. + ```Densevalm11()``` - Evaluates the log posterior density for a given parameter vector theta given the data.
Model 2

```plotmodel2()``` - makes use of clean() and model2BH() to plot a rating curve calculated with the methods in model 2.
```model2BH()``` - Makes use of Densevalm22() along with the other listed functions to determine the fit and confidence intervals of a rating curve.
```Densevalm22()``` - Evaluates the log posterior density for a given parameter vector theta given the data.
```Adist()``` - Extracts unique elements of water level measurements and creates a distance matrix from them.
```B_splines()``` - Tests the B-splines in the data.
```W_unobserved()``` - Returns the stages needed to make an equally spaced grid of stages.
```predict_u()``` - Calculates predictive values for unobserved stages.


Generalized Bayesian Power Law Model

```gbplm() ```- Infers a rating curve for paired measurements of stage and discharge using a generalized power law model described in Hrafnkelsson et al.


If you have any questions or bug fixes, please contact the developers, at solviro@gmail.com or axelorn94@gmail.com or send a pull request on the RCmodels repository on github. To download the package, type the following in your R console:
```
#Check if devtools is already installed
if(!require("devtools")){
install.packages('devtools')
}
devtools::install_github('sor16/RCmodels')
library(RCmodels)
```
Now you are ready to go and use the functions in the package. Check out the documentation for each function by typing ? in front of the function name.
