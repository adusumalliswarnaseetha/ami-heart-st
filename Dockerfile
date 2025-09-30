FROM --platform=linux/amd64 rocker/shiny:latest

# install system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev libhdf5-dev libglpk-dev libmagick++-dev libproj-dev libgdal-dev \
    fftw3-dev pkg-config

# install required R packages
RUN R -e "install.packages(c('dplyr','purr', 'zip','DT', 'Seurat', \
            'devtools', 'remotes', 'bslib','shinyjs', 'fftwtools', \
            'shinybusy', 'shinyvalidate', 'shinythemes', 'hdf5r', 'BiocManager','sfdep'))"
RUN R -e "remotes::install_github('jbergenstrahle/STUtility')"
RUN R -e "BiocManager::install('biomaRt')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('org.Ss.eg.db')"
RUN R -e "BiocManager::install('clusterProfiler')"

# make directory for the app
RUN mkdir /home/mi-heart-st

# copy the app directory
COPY ./ /home/mi-heart-st/

# set workdir
WORKDIR /home/mi-heart-st

# Expose the application port
EXPOSE 8180
  
# run app
CMD R -e "shiny::runApp('/home/mi-heart-st')"