#!/usr/bin/env Rscript

require(shiny)
require(utils)
sourceDir <- utils::getSrcDirectory(function(dummy) {dummy})
print(sourceDir)
setwd(sourceDir)
shiny::runApp('amrplusplus-shiny', launch.browser = TRUE)
