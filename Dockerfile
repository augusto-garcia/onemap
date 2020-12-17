FROM rocker/rstudio:latest

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libgfortran5 libxml2-dev libjpeg-dev libbz2-dev liblzma-dev

RUN R -e "install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/', dependencies=T);\
          devtools::install_github('Cristianetaniguti/onemap')"