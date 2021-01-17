# bdrc
Bayesian Discharge Rating Curves

The package implements the following Bayesian hierarchical discharge rating curve models for paired measurements of stage and discharge in rivers described in Hrafnkelsson et al.:

```bplm0() ``` - Power-law model with constant variance 

```bplm() ``` - Power-law model with variance that may vary with stage

```bgplm0() ```- Generalized power-law model with constant variance 

```bgplm() ```- Generalized power-law model with variance that may vary with stage

If you have any questions or bug fixes, please contact the maintainer, Birgir Hrafnkelsson at birgirhr@hi.is, or send a pull request on the RCmodels repository on github. To download the developmental version of the package, type the following in your R console:
```
#Check if devtools is already installed
if(!require("devtools")){
install.packages('devtools')
}
devtools::install_github('sor16/bdrc')
```
