FROM rocker/hadleyverse
MAINTAINER Ben Marwick <benmarwick@gmail.com>
RUN apt-get update \
  && apt-get install -y r-cran-Hmisc  r-cran-cluster r-cran-Bolstad2 r-cran-rjags \
  && installGithub.r --deps TRUE \
    mkoohafkan/kfigr \
    vqv/ggbiplot \
    rasmusab/bayesian_first_aid \
    benmarwick/Steele_et_al_VR003_MSA_Pigments \
 && mkdir /home/rstudio/Steele_et_al_VR003_MSA_Pigments \
 && git clone https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments.git  /home/rstudio/Steele_et_al_VR003_MSA_Pigments \



WORKDIR /home/rstudio/Steele_et_al_VR003_MSA_Pigments/vignettes

# RUN Rscript -e 'rmarkdown::render("/home/rstudio/Steele_et_al_VR003_MSA_Pigments/Marwick-WII13-ochre.Rmd")'
