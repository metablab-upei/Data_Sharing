FROM rocker/shiny:4.0.5

# system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2
  

# install R packages required 
RUN R -e 'install.packages(c(\
              "shiny", \
              "shinydashboard", \
              "ggplot2", \
              "rsconnect", \
              "pacman", \
              "stringr", \
              "reshape2", \
              "viridis", \
              "tidyverse", \
              "ggthemes", \
              "clipr", \
              "matrixStats", \
              "ggrepel", \
              "cowplot", \
              "readr", \
              "magrittr", \
              "dplyr", \
              "plotly", \
              "data.table", \
              "tidyr", \
              "shinyWidgets", \
              "shinyjs" \
            ), \
            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23", \
            dependencies=TRUE \
          )'


# copy the app directory into the image
COPY ./ASAppRead-only/*.R /srv/shiny-server/asapp/
COPY ./ASAppRead-only/*.csv /srv/shiny-server/asapp/

# run app
CMD ["/usr/bin/shiny-server"]
