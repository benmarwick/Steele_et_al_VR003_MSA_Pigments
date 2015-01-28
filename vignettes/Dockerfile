FROM rocker/hadleyverse
MAINTAINER Ben Marwick <benmarwick@gmail.com>
RUN apt-get update \
  && apt-get install -y r-cran-Hmisc  r-cran-cluster r-cran-Bolstad2 r-cran-rjags \
  && installGithub.r --deps TRUE \
    mkoohafkan/kfigr \
    vqv/ggbiplot \
    rasmusab/bayesian_first_aid \
    benmarwick/Steele_et_al_VR003_MSA_Pigments

# make a directory to hold the github repo
RUN mkdir /home/rstudio/Steele_et_al_VR003_MSA_Pigments

RUN git clone https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments.git   /home/rstudio/Steele_et_al_VR003_MSA_Pigments \
  && cd /home/rstudio/Steele_et_al_VR003_MSA_Pigments/vignettes

WORKDIR /home/rstudio/Steele_et_al_VR003_MSA_Pigments/vignettes

# RUN Rscript -e 'rmarkdown::render("/home/rstudio/Steele_et_al_VR003_MSA_Pigments/Marwick-WII13-ochre.Rmd")'

# to build this image:
# docker build -t benmarwick/steeleetalvr003msapigments https://raw.githubusercontent.com/benmarwick/Steele_et_al_VR003_MSA_Pigments/master/vignettes/Dockerfile

# to run this container:
# docker -dp 8787:8787 benmarwick/steeleetalvr003msapigments
# then open broswer at localhost:8787 or http://192.168.59.103:8787/

