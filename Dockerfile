
FROM rocker/r-ver:4.3.2

## Install the os packages
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    lbzip2 \
    libssl-dev \
    libudunits2-dev \
    imagemagick \
    pandoc \
    pandoc-citeproc \
    make \
    ghostscript \
    poppler-utils \
    zip \
    wget \
    fonts-dejavu-extra \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    curl \
    tk \
    libglpk-dev \ 
  && rm -rf /var/lib/apt/lists/*


RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('dplyr'); \
  install.packages('lubridate'); \
  install.packages('ggplot2'); \
  install.packages('readr'); \
  install.packages('stringr'); \
  install.packages('tidyr'); \
  install.packages('tidyverse'); \
"  


## Install extra packages required for quarto 
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    curl \
    gdebi-core \
  && rm -rf /var/lib/apt/lists/*

ARG QUARTO_VERSION="1.4.550"
RUN curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb

RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('quarto'); \
"

RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('remotes'); \
"

RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('rstan');  \ 
  install.packages('brms');   \
  install.packages('tidybayes'); 	\
  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos'))); \
  remotes::install_github('stan-dev/cmdstanr'); \
  library(cmdstanr); \
  check_cmdstan_toolchain(); \
  install_cmdstan(cores = 2); \
"  

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  cmake \
  && rm -rf /var/lib/apt/lists/*

## Other packages that are dependencies of INLA
RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('sp');   \
  install.packages('fmesher');   \
"  

## Install INLA
RUN  wget https://inla.r-inla-download.org/R/stable/src/contrib/INLA_24.02.09.tar.gz \ 
  && R CMD INSTALL --clean --no-multiarch --without-keep.source --byte-compile --resave-data --compact-docs --no-demo INLA_24.02.09.tar.gz \ 
  && rm INLA_24.02.09.tar.gz 


RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('furrr');  \
  install.packages('progressr');  \
"

RUN R -e "options(repos = \
  list(CRAN = 'https://packagemanager.posit.co/cran/2024-02-20/')); \
  install.packages('glmmTMB');   \
  install.packages('emmeans');   \
  install.packages('DHARMa');   \
  install.packages('patchwork');   \
"  
## Create project directory in docker image 
RUN mkdir ~/Project

## Copy scripts and parameters (folders and contents) into docker image project directory
COPY R/ ~/Project/R/ 
## COPY docs/ ~/Project/docs/ 
WORKDIR ~/Project/ 

## ENTRYPOINT ["make", "-i","all"]