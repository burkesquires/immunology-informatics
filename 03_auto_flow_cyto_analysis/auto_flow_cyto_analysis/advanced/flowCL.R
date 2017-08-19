### R code from vignette source 'flowCL.Rnw'

###################################################
### code chunk number 1: loadlibs
###################################################
library("flowCL")


###################################################
### code chunk number 2: Load_archive
###################################################
flowCL("Archive")


###################################################
### code chunk number 3: Date
###################################################
flowCL("Date")


###################################################
### code chunk number 4: CCR7+CD45RA+_with_visual
###################################################
Res <- flowCL("CCR7+CD45RA+")
Res$Table
tmp <- Res$'CCR7+CD45RA+'
plot(tmp[[1]], nodeAttrs=tmp[[2]], edgeAttrs=tmp[[3]], attrs=tmp[[4]])


###################################################
### code chunk number 5: CD3-CD19-CD20-CD14+
###################################################
x <-"CD3-CD19-CD20-CD14+"
Res <- flowCL(x)
Res$Cell_Label[[x]][[1]]


