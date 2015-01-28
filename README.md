<!-- README.md is generated from README.Rmd. Please edit that file -->
SteeleVR003MSAPigments
======================

<!-- DOI here -->
### Author

Ben Marwick (<benmarwick@gmail.com>)

### Contents

This repository contains the research compendium of my work on mineral pigments from Middle Stone Age archaeological sites in South Africa. The compendium contains all data, code, and text associated with this section of the publication and has been permanently archived at the DOI indicated by the above badge.

### How to use

#### Read the manuscript

See the [`vignettes`](https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments/tree/master/vignettes) directory here on GitHub for source code and data of the manuscript, and a HTML file that includes the text and figures which you can view in the browser here.

#### Install the R package

[![Build Status](https://travis-ci.org/benmarwick/Steele_et_al_VR003_MSA_Pigments.svg?branch=master)](https://travis-ci.org/benmarwick/Steele_et_al_VR003_MSA_Pigments.svg)

This repository is organized as an R package, providing functions and raw data to reproduce and extend the analysis reported in the publication. Note that this package has been written explicitly for this project and may not be suitable for more general use. To download the package source as you see it here on GitHub, use this line at the shell prompt:

``` r
git clone https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments.git
```

To install and build the package within R, use this line at the R prompt:

``` r
devtools:install_github("benmarwick/Steele_et_al_VR003_MSA_Pigments")
```

The package has several depedencies that are listed below, some of which need to be installed manually if using this package from your local R installation.

#### Run the Docker container

This compendium is also available as a [Docker](https://docs.docker.com/installation) container. The advantage of this format is that it includes this package and all its dependencies, so you don't have to worry about those. OSX & Windows users should launch [`boot2docker`](http://boot2docker.io/) to access the Docker terminal, Linux users can just open any terminal). You can either generate the Docker container yourself using the [Dockerfile](https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments/blob/master/vignettes/Dockerfile) included here, or for a quicker start, pull the image and run the container using this line at the Docker prompt:

``` r
docker run -dp 8787:8787 benmarwick/steeleetalvr003msapigments
```

Then you can interact with RStudio via your browser at localhost:8787 (on Linux) or <http://192.168.59.103:8787/> (on Windows/OSX, or whatever address boot2docker ip returns). Log in to RStudio with user: `rstudio` and password: `rstudio`. See the [rocker-org Wiki](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image) for more details. In RStudio you'll see the `Rmd` file for the manuscript and a directory for the raw data. You can knit the `Rmd` file to produce the HTML file that reproduces the text, plots and other results of the analysis. You can also edit and run the `Rmd` interactively in RStudio to explore the analysis further.

### Licenses:

Manuscript: CC-BY-4.0 <http://creativecommons.org/licenses/by/4.0/>

Code: MIT <http://opensource.org/licenses/MIT> year: 2014, copyright holder: Ben Marwick)

Data: CC0 <http://creativecommons.org/publicdomain/zero/1.0/>

### Notes and resources

-   The [issues tracker](https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments) is the place to report problems or ask questions

-   See the repository [history](https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments/commits/master) for a fine-grained view of progress and changes.

### R Session Information

``` r
sessionInfo()
#> R version 3.1.2 (2014-10-31)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.1252 
#> [2] LC_CTYPE=English_United States.1252   
#> [3] LC_MONETARY=English_United States.1252
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.1252    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#> [1] digest_0.6.8      evaluate_0.5.5    formatR_1.0       htmltools_0.2.6  
#> [5] knitr_1.8         rmarkdown_0.5.0.1 stringr_0.6.2     tools_3.1.2      
#> [9] yaml_2.1.13
```
