FROM rocker/r-devel

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libxml2-dev 

RUN R -e "install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/')"

RUN R -e "devtools::install_github('Cristianetaniguti/onemap')"

RUN R -e "install.packages('rmarkdown', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('qtl', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/')"
