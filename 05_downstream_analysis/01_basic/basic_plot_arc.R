rm(list=ls());

## setting working directory
# setwd("/Users/yourname/R");

library(OmicCircos);
options(stringsAsFactors = FALSE);

vaf <- read.table("Data_vaf.txt", header=T);
i   <- sample(1:nrow(vaf), 200);
vaf <- vaf[i,];

tran <- read.table("Data_eQTL_trans.txt", header=T);
cis  <- read.table("Data_eQTL_cis.txt", header=T);
cnv  <- read.table("TCGA.BC.cnv.2k.60.txt", header=T);
cols <- rainbow(10, alpha=0.8);

pdffile  <- "basic_plot_arc.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="basic_plot_arc.R");

circos(R=380, cir="hg19", type="chr", print.chr.lab=T, W=6);
circos(R=320, cir="hg19", W=50, mapping=vaf, col.v=4, type="b3",   B=T, col=cols[9], lwd=2);
circos(R=270, cir="hg19", W=50, mapping=vaf, col.v=4, type="s2",  B=F, col=cols, lwd=2);
circos(R=220, cir="hg19", W=50, mapping=cis, col.v=4, type="arc2",  B=T, col=cols, lwd=2);
circos(R=170, cir="hg19", W=50, mapping=cnv, col.v=60, type="b2", B=F, cutoff=0, col=cols[c(1,7)], lwd=2);
circos(R=120, cir="hg19", W=50, mapping=cis, col.v=4, type="arc", col=cols, B=F, cutoff=5, lwd=2);
circos(R=110, cir="hg19", W=40, mapping=tran[c(1:10),], type="link", lwd=2, col=cols);

dev.off()
