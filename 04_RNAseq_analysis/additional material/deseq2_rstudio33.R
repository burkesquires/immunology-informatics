setwd("C:/Users/leerkesm/Desktop/FAES/WORKSHOP_RNAseq/RNAseq_part_1")
#source("C:/Users/leerkesm/Desktop/FAES/WORKSHOP_RNAseq/RNAseq_part_1/deseq2_rstudio33.R")

source("http://bioconductor.org/biocLite.R")

biocLite("genefilter", suppressUpdates=TRUE, ask=FALSE)
library("genefilter")

### aangepast voor genefilter:
### nog checken:
### Sys.getenv("R_LIBS_USER") # value of the environment variable R_LIBS_USER
### [1] "C:\\Users\\User\\Documents/R/win-library/3.1"  
### .libPaths() # the library trees within which packages are looked for
### [1] "C:/Revolution/R-Enterprise-7.3/R-3.1.1/library" 
### http://stackoverflow.com/questions/15217758/remove-a-library-from-libpaths-permanently-without-rprofile-site

#source("http://bioconductor.org/workflows.R")
#setwd("/Users/leerkesm/Desktop/RNAseq/BATCH/class_data_ruv2_deseq2/")
#http://www.bioconductor.org/help/workflows/rnaseqGene/
#jpeg('plotTracks.jpg')
####dev.off()
#readline(prompt="Press [enter] to continue 0")

biocLite("airway", suppressUpdates=TRUE, ask=FALSE)
library("airway")
data("airway")
se <- airway

dir <- system.file("extdata", package="airway", mustWork=TRUE)

#vignette("airway")

csvfile <- file.path(dir,"sample_table.csv")
list.files(dir)
csvfile <- file.path(dir,"sample_table.csv")
sampleTable <- read.csv(csvfile,row.names=1)

filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))

biocLite("Rsamtools", suppressUpdates=TRUE, ask=FALSE)
library("Rsamtools")

#this is really important otherwise the AnnotationDbi
#will not get installed down the line:
#NB notabene: AnnotationDbi has to be installed before GenomicRanges
biocLite("AnnotationDbi", suppressUpdates=TRUE, ask=FALSE)
biocLite("org.Hs.eg.db", suppressUpdates=TRUE, ask=FALSE)

#this is really important otherwise the function summarizeOverlaps
#will not get installed down the line:
#NB notabene: GenomicAlignments has to be installed before GenomicFeatures
biocLite("GenomicAlignments", suppressUpdates=TRUE, ask=FALSE)
library("GenomicAlignments")

biocLite("GenomicFeatures", suppressUpdates=TRUE, ask=FALSE)
library("GenomicFeatures")

gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")

file.exists(filenames)

biocLite("GenomicRanges", suppressUpdates=TRUE, ask=FALSE)
library("GenomicRanges")

bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])

biocLite("GenomicAlignments", suppressUpdates=TRUE, ask=FALSE)
library("GenomicAlignments")

txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
ebg <- exonsBy(txdb, by="gene")

#library("DESeq2")
#library("GenomicFeatures")
#library("Rsamtools")
#library("GenomicAlignments")
#library("GenomicRanges”)

.libPaths()
#> .libPaths()
#[1] "C:/Users/leerkesm/Documents/R/R-3.3.2/library"
#.libPaths("C:/Users/leerkesm/Documents/R/R-3.3.2/library")

#solution is here:
#http://stackoverflow.com/questions/19407092/r-not-finding-package-even-after-package-installation

#http://stackoverflow.com/questions/26570912/error-in-installation-a-r-package
#http://seqanswers.com/forums/showthread.php?p=142291
#https://support.bioconductor.org/p/63875/

library(BiocParallel)
register(SerialParam())

#se <- summarizeOverlaps(features=ebg, reads=bamfiles,
#                        ignore.strand=TRUE, BPPARAM=SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

#se <- summarizeOverlaps(features=ebg, reads=bamfiles,
#                        mode="Union",
#                        singleEnd=FALSE,
#                        ignore.strand=TRUE,
#                        fragments=TRUE )

assayNames(se)
head(assay(se), 3)
colSums(assay(se))

colSums(assay(se))
rowRanges(se)

data("airway")
se <- airway

str(metadata(rowRanges(se)))
colData(se)
colData(se) <- DataFrame(sampleTable)
se$cell
se$dex

str(se)

se$dex <- relevel(se$dex, "untrt")
se$dex

round( colSums(assay(se)) / 1e6, 1 )

colData(se)

biocLite("DESeq2", suppressUpdates=TRUE, ask=FALSE)
library("DESeq2")

dds <- DESeqDataSet(se, design = ~ cell + dex)
countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)

plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.3)
readline(prompt="Press [enter] to continue 1")

plot(assay(rld)[,1:2], pch=16, cex=0.3)
readline(prompt="Press [enter] to continue 2")

sampleDists <- dist( t( assay(rld) ) )

biocLite("pheatmap", suppressUpdates=TRUE, ask=FALSE)
library("pheatmap")

biocLite("RColorBrewer", suppressUpdates=TRUE, ask=FALSE)
library("RColorBrewer")

#Sample distances

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
readline(prompt="Press [enter] to continue 3")

#check_outliers

plotPCA(rld, intgroup = c("dex", "cell"))
readline(prompt="Press [enter] to continue 4")
#dev.off()

data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE)

percentVar <- round(100 * attr(data, "percentVar"))

biocLite("ggplot2", suppressUpdates=TRUE, ask=FALSE)
library("ggplot2")

#PCA plot using the rlog-transformed values. 

ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
readline(prompt="Press [enter] to continue 5")

#MDS plot

#check_outliers

mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)
readline(prompt="Press [enter] to continue 6")
#dev.off()

#Running the differential expression pipeline

dds <- DESeq(dds)

#Building the results table

res <- results(dds)

mcols(res, use.names=TRUE)

summary(res)

res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.001)

#Other comparisons

results(dds, contrast=c("cell", "N061011", "N61311"))

sum(res$pvalue < 0.05, na.rm=TRUE)

sum(!is.na(res$pvalue))

sum(res$padj < 0.1, na.rm=TRUE)

resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

#Plotting results

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))
readline(prompt="Press [enter] to continue 7")
#dev.off()

#Normalized counts for a single gene over treatment group.

data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE)
ggplot(data, aes(x=dex, y=count, color=cell)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3)
readline(prompt="Press [enter] to continue 8")
#dev.off()

#Normalized counts indicating cell line with color.

ggplot(data, aes(x=dex, y=count, fill=dex)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")
readline(prompt="Press [enter] to continue 9")
#dev.off()

#Normalized counts using a more structural arrangement. Here the color indicates treatment.

ggplot(data, aes(x=dex, y=count, color=cell, group=cell)) +
  scale_y_log10() + geom_point(size=3) + geom_line()
readline(prompt="Press [enter] to continue 10")
#just so that the before figure will not overlap with the next figure, do dev.off
dev.off()

#An MA-plot

plotMA(res, ylim=c(-5,5))
readline(prompt="Press [enter] to continue 11")
#dev.off()

plotMA(resLFC1, ylim=c(-5,5))
readline(prompt="Press [enter] to continue 12")
#dev.off()

topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Histogram of p values for genes with mean normalized count larger than 1.

hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")
readline(prompt="Press [enter] to continue 13")
#dev.off()

#Gene clustering

Sys.getenv("R_LIBS_USER")
.libPaths()

#> .libPaths()
#[1] "C:/Users/leerkesm/Documents/R/R-3.3.2/library"
#.libPaths("C:/Users/leerkesm/Documents/R/R-3.3.2/library")
#https://support.bioconductor.org/p/59290/

#There could be a few things happening here. Start by first figuring out your library location:
#Sys.getenv("R_LIBS_USER")
#http://stackoverflow.com/questions/26570912/error-in-installation-a-r-package
#library(caret, lib.loc = "C:/Users/jkuusik/Documents/R/win-library/2.14")

#http://r.789695.n4.nabble.com/xts-zoo-quot-cannot-remove-prior-installation-of-package-quot-td4671096.html
#If you look at the help page for remove.packages(), it states that the second argument is taken to be the first element in .libPaths(). You have not told us what .libPaths() is returning and you should do so (even though it wasn't mentioned in the Posting Guide.)  It seems possible that you set up an additonal library location which your new installation of R is "finding" but which is not the first item in the list of directories with library entries. 
#It seems like there was an issue with the libPath. The path is setup to reference a cached copy of a shared drive and somehow the cache and the shared drive were not synching. I finally located the real place (i.e. the cached area) where R was looking and manually removed everything. 

#?remove.packages
#http://127.0.0.1:24146/library/utils/html/remove.packages.html
#remove.packages(pkgs, lib)

#remove.packages(genefilter, "C:/Users/leerkesm/Documents/R/R-3.3.2/library")

#Installing R packages without admin rights on MS Windows 
#http://www.magesblog.com/2012/04/installing-r-packages-without-admin.html

#http://stackoverflow.com/questions/2615128/where-does-r-store-packages

#installed.packages(lib.loc = NULL, priority = NULL,
#                   noCache = FALSE, fields = NULL,
#                   subarch = .Platform$r_arch)

#http://stackoverflow.com/questions/2615128/where-does-r-store-packages
#.GlobalEnv
#<environment: R_GlobalEnv>
#http://stackoverflow.com/questions/21209230/assigning-a-locked-variable-in-an-r-package
#https://gist.github.com/wch/3280369
#https://revolutionanalytics.zendesk.com/hc/en-us/community/posts/203901043-Can-not-install-caret-package

#install as type source, then change library path
#biocLite("IRanges", type="source")
#The downloaded source packages are in
#        ‘C:\Users\leerkesm\AppData\Local\Temp\1\RtmpIRJ6y2\downloaded_packages’
#Warning messages:
#1: running command '"C:/Users/leerkesm/Documents/R/R-3.3.2/bin/x64/R" CMD INSTALL -l "C:\Users\leerkesm\Documents\R\R-3.3.2\library" C:\Users\leerkesm\AppData\Local\Temp\1\RtmpIRJ6y2/downloaded_packages/genefilter_1.56.0.tar.gz' had status 1 
#2: In install.packages(pkgs = doing, lib = lib, ...) :
#  installation of package ‘genefilter’ had non-zero exit status

#To install a package, I have to specify a library location:
#install.packages("zoo", lib="C:/software/Rpackages")

#To load a package, I also have to specify the library location:
#library("genefilter", lib.loc="C:/software/Rpackages")

#library("zoo", lib.loc="C:/Users/leerkesm/Documents/R/R-3.3.2/library")

#try to find workaround for system to find libraries in path
#http://www.magesblog.com/2012/04/installing-r-packages-without-admin.html
#.libPaths(c("C:\\Users\\MyNAME\\R", .libPaths()))

#.libPaths(c("C:/Users/leerkesm/Documents/R/R-3.3.2/library/", .libPaths()))
########################################biocLite("genefilter", suppressUpdates=TRUE, ask=FALSE)
#biocLite("genefilter", type="source", suppressUpdates=TRUE, ask=FALSE)
#library(genefilter, lib.loc = "C:/Users/leerkesm/Documents/R/R-3.3.2/library/genefilter")
#library(genefilter, lib.loc = "C:/Users/leerkesm/Documents/R/R-3.3.2/library/genefilter")
#biocLite("genefilter", type="source", suppressUpdates=TRUE, ask=FALSE)

####execute in beginning, delete first from libraries
####biocLite("genefilter", suppressUpdates=TRUE, ask=FALSE)
####library("genefilter")

#library("genefilter", lib.loc="C:/software/Rpackages")

topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)

mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])

#jpeg('rplotpheatmap.jpg')
pheatmap(mat, annotation_col=df)
readline(prompt="Press [enter] to continue pheatmap14")
#dev.off()

#Independent filtering

qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
readline(prompt="Press [enter] to continue 15")
#dev.off()

#Annotating and exporting results

#biocLite("AnnotationDbi", suppressUpdates=TRUE, ask=FALSE)
#biocLite("org.Hs.eg.db", suppressUpdates=TRUE, ask=FALSE)

library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$padj),]

#Exporting results

resOrderedDF <- as.data.frame(resOrdered)[1:100,]
#write.csv(resOrderedDF, file="myresults.csv")
#write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="logFC", number=nrow(d_vsn)), file="blm-cd11.xls", row.names=F, sep="\t")
write.table(resOrderedDF, file="myresults33_nhlbi.xls", row.names=F, sep="\t")
getwd()

#Plotting fold changes in genomic space
resGR <- results(dds, lfcThreshold=1, format="GRanges")
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

biocLite("Gviz", suppressUpdates=TRUE, ask=FALSE)
library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))

options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
               type="h", name="log2 fold change", strand="+")

#log2 fold changes in genomic region surrounding the gene with smallest adjusted p value. Genes highlighted in pink have adjusted p value less than 0.1.


#jpeg('plotTracks.jpg')
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
readline(prompt="Press [enter] to continue pheatmap16a")
#dev.off()

#readline(prompt="Press [enter] to continue 16 svaseq now")

#Below we obtain a matrix of normalized counts for which the average count across 
#samples is larger than 1. As we described above, we are trying to recover any hidden 
#batch effects, supposing that we do not know the cell line information. So we use a 
#full model matrix with the dex variable, and a reduced, or null, model matrix with only 
#an intercept term. Finally we specify that we want to estimate 2 surrogate variables. 
#For more information read the manual page for the svaseq function by typing ?svaseq.

biocLite("sva", suppressUpdates=TRUE, ask=FALSE)
library("sva")
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)

svseq$sv

#Because we actually do know the cell lines, we can see how well the SVA method did 
#at recovering these variables (Figure below).

#jpeg('plotTracks.jpg')
#plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
##dev.off()

#jpeg('stripchart1b.jpg')
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
abline(h=0)
readline(prompt="Press [enter] to continue 17 svaseq now")
#dev.off()

#jpeg('stripchart2b.jpg')
stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
abline(h=0)
readline(prompt="Press [enter] to continue 18 svaseq now")
#dev.off()

#Surrogate variables 1 and 2 plotted over cell line. Here, we know the hidden source of variation (cell line), and therefore can see how the SVA procedure is able to identify a source of variation which is correlated with cell line.
#Finally, in order to use SVA to remove any effect on the counts from our surrogate variables, we simply add these two surrogate variables as columns to the DESeqDataSet and then add them to the design:

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex

#We could then produce results controlling for surrogate variables by running DESeq with the new design:
#ddssva <- DESeq(ddssva)
#write.table(ddssva, file="ddssva.xls")
#We could then produce results controlling for surrogate variables by running DESeq with the new design:
ddssva <- DESeq(ddssva)
#Building the results table
res_sva <- results(ddssva)
mcols(res_sva, use.names=TRUE)
summary(res_sva)

res_sva.05 <- results(ddssva, alpha=.05)
table(res_sva.05$padj < .05)

res_svaLFC1 <- results(ddssva, lfcThreshold=1)
table(res_svaLFC1$padj < 0.001)

res_sva$symbol <- mapIds(org.Hs.eg.db,
                         keys=row.names(res_sva),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
res_sva$entrez <- mapIds(org.Hs.eg.db,
                         keys=row.names(res_sva),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")

res_svaOrdered <- res_sva[order(res_sva$padj),]

#Exporting results

res_svaOrderedDF <- as.data.frame(res_svaOrdered)[1:100,]
write.table(res_svaOrderedDF, file="my_ddssva_results33_nhlbi.xls", row.names=F, sep="\t")
