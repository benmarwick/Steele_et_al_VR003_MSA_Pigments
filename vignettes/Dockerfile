# get the base image, this one has R, RStudio, pandoc, and a bunch of R packages that I use often
FROM rocker/hadleyverse

# required
MAINTAINER Ben Marwick <benmarwick@gmail.com>

# install some packages that not in the base image
RUN apt-get update \
  && apt-get install -y r-cran-Hmisc  r-cran-cluster r-cran-Bolstad2 r-cran-rjags \
  # install a few packages from GitHub (not on CRAN)
  && installGithub.r --deps TRUE \
    mkoohafkan/kfigr \
    vqv/ggbiplot \
    rasmusab/bayesian_first_aid \
    # install my package that is the focus of this image
    benmarwick/Steele_et_al_VR003_MSA_Pigments

# make a directory to hold the github repo
RUN mkdir /home/rstudio/Steele_et_al_VR003_MSA_Pigments

# Get repository from github and onto Docker
RUN git clone https://github.com/benmarwick/Steele_et_al_VR003_MSA_Pigments.git   /home/rstudio/Steele_et_al_VR003_MSA_Pigments

# Keep only the vignette folder that has the Rmd file and data files
RUN cd /home/rstudio/Steele_et_al_VR003_MSA_Pigments \
  && cp vignettes ../rstudio \
  && cd .. \
  && rm -r Steele_et_al_VR003_MSA_Pigments

# Set the working directory
WORKDIR /home/rstudio/vingettes/

# Open the package in RStudio
# RUN rstudio ~/Steele_et_al_VR003_MSA_Pigments.Rproj

# RUN Rscript -e 'rmarkdown::render("/home/rstudio/Steele_et_al_VR003_MSA_Pigments/Marwick-WII13-ochre.Rmd")'

# to build this image:
# docker build -t benmarwick/steeleetalvr003msapigments https://raw.githubusercontent.com/benmarwick/Steele_et_al_VR003_MSA_Pigments/master/vignettes/Dockerfile

# to run this container:
# docker -dp 8787:8787 benmarwick/steeleetalvr003msapigments
# then open broswer at localhost:8787 or http://192.168.59.103:8787/

