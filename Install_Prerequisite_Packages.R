#!/usr/bin/env Rscript

install.packages(c('shiny', 'shinydashboard', 'data.table', 'ggplot2',
                   'vegan'), repos='http://cran.us.r-project.org')
source("https://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
