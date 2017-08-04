rm(list=ls());
dev.off()

## setting working directory
# setwd("/Users/yourname/R");

library(OmicCircos);
options(stringsAsFactors = FALSE);

gene <- read.table("TCGA.PAM50_subset.txt", header=T);
tran <- read.table("Data_eQTL_trans.txt", header=T);
cnv  <- read.table("TCGA.BC.cnv.2k.60.txt", header=T);
cols <- rainbow(10, alpha=0.8);

pdffile  <- "basic_plot_lable1.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="basic_plot_label1.R");

circos(R=300, cir="hg19", type="chr", print.chr.lab=F, W=10);
circos(R=310, cir="hg19", W=20, mapping=gene, type="label", side="out", col=c("black", "blue"), cex=c(0.4,0.6));
circos(R=230, cir="hg19", W=60, mapping=cnv, col.v=60, type="b2", B=T, cutoff=0, col=cols[c(1,7)], lwd=2);
circos(R=150, cir="hg19", W=60, mapping=cnv, col.v=60, type="ml3", B=F, cutoff=0, col=cols[c(1,7)], lwd=2);
circos(R=140, cir="hg19", W=40, mapping=tran[c(1:10),], type="link", lwd=2, col=cols);

dev.off()
