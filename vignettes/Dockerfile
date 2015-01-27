FROM rocker/hadleyverse
MAINTAINER Ben Marwick <benmarwick@gmail.com>
RUN apt-get update \
  && apt-get install -y r-cran-Hmisc r-cran-fpc r-cran-cluster r-cran-Bolstad2 \
  && installGithub.r --deps TRUE \
    mkoohafkan/kfigr \
    vqv/ggbiplot \
    rasmusab/bayesian_first_aid \
    benmarwick/Steele_et_al_VR003_MSA_Pigments

COPY . /home/rstudio/Steele_et_al_VR003_MSA_Pigments
WORKDIR /home/rstudio/Steele_et_al_VR003_MSA_Pigments

# RUN Rscript -e 'rmarkdown::render("/home/rstudio/Steele_et_al_VR003_MSA_Pigments/Marwick-WII13-ochre.Rmd")'
