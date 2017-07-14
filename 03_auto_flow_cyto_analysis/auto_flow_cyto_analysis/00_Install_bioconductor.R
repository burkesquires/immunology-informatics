# Bioconductor installs

source("https://bioconductor.org/biocLite.R")

biocLite ()

biocLite("flowCore")
biocLite("flowType")
biocLite("RchyOptimyx")
biocLite("flowQ")

install.packages("ROCR")
install.packages("sfsmisc")
install.packages("GEOmap")

install.packages("car"); install.packages("RFOC")

library(flowCore)
library(flowDensity)
library(flowType)
library(RchyOptimyx);
library(ROCR)
library(sfsmisc)
library(GEOmap)


# FOR GATING

biocLite('rrcov')
biocLite('codetools')
biocLite('foreach')
biocLite('flowMerge')
biocLite('flowType')
biocLite('RchyOptimyx')
