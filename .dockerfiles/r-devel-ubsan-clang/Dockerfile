FROM cristaniguti/r-devel-ubsan-clang

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libxml2-dev pandoc qpdf

RUN R -e "install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/')"

RUN R -e "install.packages('rmarkdown', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('qtl', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('survival', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('Matrix', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('lattice', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('cluster', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('foreign', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('nnet', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('rpart', repos='http://brieger.esalq.usp.br/CRAN/')"
RUN R -e "install.packages('onemap', repos='http://brieger.esalq.usp.br/CRAN/')"
