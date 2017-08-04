rm(list=ls());
dev.off()

## setting working directory
# setwd("/Users/yourname/R");

library(OmicCircos);
options(stringsAsFactors = FALSE);

set.seed(1234);

data("UCSC.hg19.chr");

expr <- read.table("TCGA.BC.gene.exp.2k.60.txt", header=T);
i    <- sample(1:nrow(expr), 400);
expr <- expr[i,]
tran <- read.table("Data_eQTL_trans.txt", header=T);
cnv  <- read.table("TCGA.BC.cnv.2k.60.txt", header=T);
col1 <- rainbow(23);
col2 <- rainbow(10, alpha=0.5);

pdf("advance_zoom1.pdf", 8, 8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="advance_zoom1.R");

angle.start <- 0;
angle.end   <- 360;

chr.i     <- c(2,4,9,8,17,22)
chr.s     <- paste0("chr", chr.i);

## select chromsomes in zoom regions
chr.i1    <- which(tran[,1] %in% chr.i);
chr.i2    <- which(tran[,4] %in% chr.i);
chr.j     <- intersect(chr.i1, chr.i2);
chr.db <- segAnglePo(seg.dat=UCSC.hg19.chr, seg=chr.s, 
                     angle.start=angle.start, angle.end=angle.end);
circos(R=390, type="chr", cir=chr.db, W=10, print.chr.lab=T, col=col1[chr.i]);
circos(R=290, cir=chr.db, W=100, mapping=expr,  col.v=50,  type="heatmap2", 
       cluster=F, col.bar=F, lwd=0.1, col="blue");
circos(R=230, cir=chr.db, W=60,  mapping=cnv,  col.v=50,  type="ml3", B=F, lwd=1, cutoff=0);
circos(R=170, cir=chr.db, W=60,  mapping=cnv,  col.v=50,  type="b2",  cutoff=0, B=T, lwd=1, col=col2[c(7,9)]);
circos(R=160, cir=chr.db, W=40, mapping=tran[chr.j,], type="link", lwd=2, col=col2);

dev.off()