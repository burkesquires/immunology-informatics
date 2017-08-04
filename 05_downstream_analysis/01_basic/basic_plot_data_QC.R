rm(list=ls());

## setting working directory
# setwd("/Users/yourname/R");
setwd("C:/Users/leerkesm/Desktop/OmicCircos/all_files")
getwd();

## loading cnv data
cnv <- read.table("TCGA.BC.cnv.2k.60.txt", header=T); 
## row/column number
dim(cnv)
## data structure
str(cnv)
## list first six rows and columns
head(cnv[,c(1:6)])

## data distributions
cnv.s  <- cnv[,c(4:ncol(cnv))];
summary(cnv.s[,c(1:6)]);

## output pdf file
pdf("basic_plot_data_QC.pdf", 6, 6);
## boxplot
boxplot(cnv.s, col=rainbow(ncol(cnv.s)));
## hist
hist(cnv.s[,1], col="lightblue");

## count number per chr
table(cnv[,1]);
barplot(table(cnv[,1]), col=rainbow(23), xlab="chromosome", ylab="Number")

## loading expression data
expr <- read.table("TCGA.BC.gene.exp.2k.60.txt", header=T);
expr.s <- expr[,4:ncol(expr)];

## loading sample information
tumor.type <- read.table("TCGA.BC.sample60.txt", header=T);
colnames(cnv.s) == tumor.type[,1];
colnames(expr.s) == tumor.type[,1];

## pca 
pca    <- prcomp(t(expr.s), scale=T);
cols <- rainbow(10, alpha=0.8)[c(1,2,7,9)];
plot(pca$x[,1], pca$x[,2], col=cols[as.numeric(as.factor(tumor.type[,2]))], pch=19, xlab="PC1", ylab="PC2");
legend("topright", as.character(unique(tumor.type[,2])), col=cols, pch=19)

## cluster
expr.d <- dist(t(expr.s));
expr.h <- hclust(expr.d);
plot(expr.h, cex=0.6);

dev.off()

