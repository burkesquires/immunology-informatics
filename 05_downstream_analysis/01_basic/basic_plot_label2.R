rm(list=ls());

## setting working directory
# setwd("/Users/yourname/R");

library(OmicCircos);
options(stringsAsFactors = FALSE);

gene <- read.table("TCGA.PAM50_subset.txt", header=T);
tran <- read.table("Data_eQTL_trans.txt", header=T);
cnv  <- read.table("TCGA.BC.cnv.2k.60.txt", header=T);
cols <- rainbow(10, alpha=0.8);

pdffile  <- "basic_plot_label2.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="basic_plot_label2.R");

circos(R=310, cir="hg19", W=80, mapping=cnv, col.v=60, type="b2", B=T, cutoff=0, col=cols[c(1,7)], lwd=2);
circos(R=300, cir="hg19", type="chr", print.chr.lab=T, W=6);
circos(R=290, cir="hg19", W=20, mapping=gene, type="label", side="in", col=c("black", "blue"), cex=c(0.4,0.6));
circos(R=130, cir="hg19", W=80, mapping=cnv, col.v=60, type="ml3", B=T, cutoff=0, col=cols[c(1,7)], lwd=2);
circos(R=110, cir="hg19", W=40, mapping=tran[c(1:10),], type="link", lwd=2, col=cols);

dev.off()
