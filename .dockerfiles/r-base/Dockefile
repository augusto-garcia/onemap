FROM r-base:3.6.0

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libxml2-dev

RUN R -e "install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/');\
          devtools::install_github('Cristianetaniguti/onemap')"