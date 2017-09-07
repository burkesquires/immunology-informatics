### R code from vignette source 'FlowSOM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: FlowSOM.Rnw:59-76
###################################################
library(FlowSOM)

fileName <- system.file("extdata","lymphocytes.fcs",
                        package="FlowSOM")
fSOM <- FlowSOM(fileName,
                # Input options:
                compensate = TRUE,transform = TRUE,toTransform=c(8:18),
                scale = TRUE,
                # SOM options:
                colsToUse = c(9,12,14:18), xdim = 7, ydim = 7,
                # Metaclustering options:
                nClus = 10,
                # Seed for reproducible results:
                seed = 42)

PlotStars(fSOM$FlowSOM,
            backgroundValues = as.factor(fSOM$metaclustering))


###################################################
### code chunk number 3: FlowSOM.Rnw:109-120
###################################################
set.seed(42)
library(flowCore)
library(FlowSOM)

fileName <- system.file("extdata","lymphocytes.fcs",
                        package="FlowSOM")
fSOM <- ReadInput(fileName,compensate = TRUE,transform = TRUE, 
                    toTransform=c(8:18),scale = TRUE)

ff <- suppressWarnings(flowCore::read.FCS(fileName))
fSOM <- ReadInput(ff,compensate = TRUE,transform = TRUE, scale = TRUE)


###################################################
### code chunk number 4: FlowSOM.Rnw:131-132
###################################################
str(fSOM,max.level = 2)


###################################################
### code chunk number 5: FlowSOM.Rnw:150-152
###################################################
fSOM <- BuildSOM(fSOM,colsToUse = c(9,12,14:18))
str(fSOM$map,max.level = 2)


###################################################
### code chunk number 6: FlowSOM.Rnw:162-164
###################################################
fSOM <- BuildMST(fSOM,tSNE=TRUE)
str(fSOM$MST)


###################################################
### code chunk number 7: FlowSOM.Rnw:172-173
###################################################
PlotStars(fSOM)


###################################################
### code chunk number 8: FlowSOM.Rnw:175-176
###################################################
PlotStars(fSOM,view="grid")


###################################################
### code chunk number 9: FlowSOM.Rnw:178-179
###################################################
PlotStars(fSOM,view="tSNE")


###################################################
### code chunk number 10: FlowSOM.Rnw:185-189
###################################################
fSOM <- UpdateNodeSize(fSOM, reset=TRUE)
fSOM$MST$size <- fSOM$MST$size/2
PlotStars(fSOM)
fSOM <- UpdateNodeSize(fSOM)


###################################################
### code chunk number 11: FlowSOM.Rnw:196-212
###################################################
#<<>>= 
library(flowUtils)
flowEnv <- new.env()
ff_c <- compensate(ff,description(ff)$SPILL)
colnames(ff_c)[8:18] <- paste("Comp-",colnames(ff_c)[8:18],sep="")
gatingFile <- system.file("extdata","manualGating.xml", 
                        package="FlowSOM")
gateIDs <- c( "B cells"=8,
                "ab T cells"=10,
                "yd T cells"=15,
                "NK cells"=5,
                "NKT cells"=6)
cellTypes <- names(gateIDs)
gatingResult <- ProcessGatingML(ff_c, gatingFile, gateIDs, cellTypes)

PlotPies(fSOM,cellTypes=gatingResult$manual)


###################################################
### code chunk number 12: FlowSOM.Rnw:218-220
###################################################
print(colnames(fSOM$map$medianValues))
PlotMarker(fSOM,"Pacific Blue-A")


###################################################
### code chunk number 13: FlowSOM.Rnw:224-225
###################################################
PlotNumbers(UpdateNodeSize(fSOM,reset=TRUE))


###################################################
### code chunk number 14: FlowSOM.Rnw:229-230
###################################################
PlotClusters2D(fSOM,"PE-Texas Red-A","Pacific Blue-A",c(81,82,91,92,93))


###################################################
### code chunk number 15: FlowSOM.Rnw:243-247
###################################################
#<<>>=
metaClustering <- metaClustering_consensus(fSOM$map$codes,k=7)
PlotPies(fSOM,cellTypes=gatingResult$manual,
        backgroundValues = as.factor(metaClustering))


###################################################
### code chunk number 16: FlowSOM.Rnw:251-252
###################################################
metaClustering_perCell <- metaClustering[fSOM$map$mapping[,1]]


###################################################
### code chunk number 17: FlowSOM.Rnw:262-274
###################################################
# Look for CD8+ ab T cells
query <- c("PE-Cy7-A" = "high", #CD3
            "APC-Cy7-A" = "high", #TCRb
            "Pacific Blue-A" = "high") #CD8
query_res <- QueryStarPlot(UpdateNodeSize(fSOM,reset=TRUE), query, 
                            plot = FALSE)

cellTypes <- factor(rep("Unknown",49),levels=c("Unknown","CD8 T cells"))
cellTypes[query_res$selected] <- "CD8 T cells"
PlotStars(fSOM,
            backgroundValues=cellTypes,
            backgroundColor=c("#FFFFFF00","#0000FF22"))


###################################################
### code chunk number 18: FlowSOM.Rnw:280-313
###################################################
library(FlowSOM)

# Build the FlowSOM tree on the example file
fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
flowSOM.res <- FlowSOM(fileName, compensate=TRUE,transform=TRUE,
scale=TRUE,colsToUse=c(9,12,14:18),nClus = 10, seed=1)

# Have a look at the resulting tree
# PlotStars(flowSOM.res[[1]],backgroundValues = as.factor(flowSOM.res[[2]]))

# Select all cells except the branch that corresponds with automated
# cluster 7 (CD3+ TCRyd +) and write te another file for the example
# In practice you would not generate any new file but 
# use your different files from your different groups
ff <- flowCore::read.FCS(fileName)
ff_tmp <- ff[flowSOM.res[[1]]$map$mapping[,1] %in% which(flowSOM.res[[2]] != 7),]
flowCore::write.FCS(ff_tmp,file="ff_tmp.fcs")
# Make an additional file without cluster 7 and double amount of cluster 10
ff_tmp <- ff[c(which(flowSOM.res[[1]]$map$mapping[,1] %in% which(flowSOM.res[[2]] != 7)),
which(flowSOM.res[[1]]$map$mapping[,1] %in% which(flowSOM.res[[2]] == 5))),]
flowCore::write.FCS(ff_tmp,file="ff_tmp2.fcs")

# Compare the original file with the two new files we made
groupRes <- CountGroups(flowSOM.res[[1]], 
                        groups=list("AllCells"=c(fileName),
                                    "Without_ydTcells"=c("ff_tmp.fcs","ff_tmp2.fcs")))
# PlotGroups(flowSOM.res[[1]], groupRes)

# Compare only the file with the double amount of cluster 10
groupRes <- CountGroups(flowSOM.res[[1]], 
                        groups=list("AllCells"=c(fileName),
                                    "Without_ydTcells"=c("ff_tmp2.fcs")))
PlotGroups(flowSOM.res[[1]], groupRes)


###################################################
### code chunk number 19: FlowSOM.Rnw:325-326
###################################################
sessionInfo()


