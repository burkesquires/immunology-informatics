rm(list=ls());

## setting working directory
# setwd("/Users/yourname/R");

library(OmicCircos);
options(stringsAsFactors = FALSE);

expr <- read.table("TCGA.BC.gene.exp.2k.60.txt", header=T);
expr <- expr[-nrow(expr),];
i    <- sample(1:nrow(expr), 400);
expr <- expr[i,]
tran <- read.table("Data_eQTL_trans.txt", header=T);
cnv  <- read.table("TCGA.BC.cnv.2k.60.txt", header=T);
cols <- rainbow(10, alpha=0.8);

pdffile  <- "advance_heatmap.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="advance_heatmap.R");

circos(R=390, cir="hg19", type="chr", print.chr.lab=T, W=10);
circos(R=290, cir="hg19", W=100, mapping=expr,  col.v=20,  type="heatmap2", 
       cluster=T, col.bar=T, col.bar.po="bottomright", lwd=0.1, col="blue");
circos(R=230, cir="hg19", W=60,  mapping=cnv,  col.v=60,  type="ml3", B=F, lwd=1, cutoff=0);
circos(R=170, cir="hg19", W=60,  mapping=cnv,  col.v=61,  type="b2",  cutoff=0, B=T, lwd=1, col=cols[c(7,9)]);

circos(R=160, cir="hg19", W=40, mapping=tran[c(1:20),], type="link", lwd=2, col=cols);

dev.off()
