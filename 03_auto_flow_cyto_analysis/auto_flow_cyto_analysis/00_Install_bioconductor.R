# Bioconductor installs
# R. Burke Squires
# Source: https://github.com/burkesquires/immunology-informatics/tree/master/03_auto_flow_cyto_analysis/auto_flow_cyto_analysis

################################################################################
# Normally to install package using R and Bioconductor you will type:
################################################################################

# You tell R where to find the source code for the Bioconductor packages

source("https://bioconductor.org/biocLite.R")

# Now that R now where to find the package information,
# you type _biocLite()_ to retrieve teh most recent list of packages

biocLite ()

# Now that R nows what packages are available, we can install a package.
# Lets install the _flowCore_ package

biocLite("flowCore")

# You only need to install a package once. You can update it later on by installing again,
# which will install the latest version

# Before you can actually USE a package you have to load it using the library command like this:

library(flowCore)

# We can install all of our packages

biocLite('flowMerge')
biocLite("flowQ")
biocLite('flowType')

biocLite('RchyOptimyx')

biocLite('rrcov')
biocLite('codetools')
biocLite('foreach')


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

# Before leaving for lunch

source("http://bioconductor.org/biocLite.R")
biocLite('rrcov')
biocLite('codetools')
biocLite('foreach')
biocLite('flowMerge')
biocLite('flowType')
biocLite('RchyOptimyx')
biocLite('Rcpp')

library(flowType)
library(RchyOptimyx)



# This is a function that just makes sure you have a package, or installs it for you without prompting
# Gresham Lab Flow Core Guide
# https://github.com/GreshamLab/flow

requireInstall <- function(packageName,isBioconductor=F) {
  if ( !try(require(packageName,character.only=T)) ) {
    print(paste0("You don't have ",packageName," accessible, ",
                 "I'm gonna install it"))
    if (isBioconductor) {
      source("http://bioconductor.org/biocLite.R")
      biocLite(packageName)
    } else {
      install.packages("packageName", repos = "http://cran.us.r-project.org")
    }
  }
  return(1)
}

#Load libraries
requireInstall("flowCore",isBioconductor=T)
