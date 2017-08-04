rm(list=ls());
dev.off()

## setting working directory
#setwd("/Users/yhu/Documents/Projects/CBIIT/OmicCircos_NIH/20160405/R");

library(OmicCircos);
options(stringsAsFactors = FALSE);

expr <- read.table("TCGA.BC.gene.exp.2k.60.txt", header=T);
tran <- read.table("Data_eQTL_trans.txt", header=T);
cols <- rainbow(10, alpha=0.8);

pdffile  <- "basic_plot_varitation.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="basic_plot_variation.R");

circos(R=380, cir="hg19", type="chr", print.chr.lab=T, W=6);
circos(R=330, cir="hg19", W=40, mapping=expr, col.v=4, type="quant90",   B=T, col=cols[1], lwd=2);
circos(R=290, cir="hg19", W=40, mapping=expr, col.v=4, type="sv",  B=F, col=cols[9], lwd=2);
circos(R=250, cir="hg19", W=40, mapping=expr, col.v=4, type="ss",  B=T, col=cols[7], lwd=2);
circos(R=210, cir="hg19", W=40, mapping=expr, col.v=60, type="heatmap", B=F, cutoff=5, lwd=2);
circos(R=170, cir="hg19", W=40, mapping=expr, col.v=60, type="s.sd", B=F, cutoff=5, lwd=2);
circos(R=130, cir="hg19", W=40, mapping=expr, col.v=60, type="ci95", B=F, cutoff=5, lwd=2);
circos(R=120, cir="hg19", W=40, mapping=tran[c(1:10),], type="link", lwd=2, col=cols);

dev.off()
