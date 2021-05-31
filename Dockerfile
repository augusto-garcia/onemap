FROM rocker/rstudio:4.0.1

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libgfortran5 libxml2-dev libjpeg-dev libbz2-dev liblzma-dev libgit2-dev libfontconfig1-dev 

RUN R -e "install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/', dependencies=T)"
          
RUN R -e "devtools::install_github('Cristianetaniguti/onemap')"
