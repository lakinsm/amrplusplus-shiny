#!/usr/bin/env Rscript

install.packages(c('shiny', 'shinydashboard', 'data.table', 'ggplot2',
                   'vegan'))
source("https://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")