rm(list=ls());

library(OmicCircos);
options(stringsAsFactors = FALSE);
data("UCSC.hg19");

## read G banding data
hg19.b <- read.table("cytoBand.txt", header=T);
## simulation value
hg19.b <- data.frame(hg19.b, value=rnorm(nrow(hg19.b)));

## for G banding colors
col.c     <- NULL;
c1        <- which(hg19.b[,5]=="gneg");
col.c[c1] <- colors()[351];
c2        <- which(hg19.b[,5]=="acen" | hg19.b[,5]=="gvar" | hg19.b[,5]=="stalk");
col.c[c2] <- colors()[26];
c3        <- which(hg19.b[,5]=="gpos50" | hg19.b[,5]=="gpos75" | hg19.b[,5]=="gpos100");
col.c[c3] <- colors()[300];
c4        <- which(hg19.b[,5]=="gpos25");
col.c[c4] <- colors()[351];

col2   <- rainbow(10, alpha=0.5);

## chromosome labels
chr.l   <- data.frame(chr=UCSC.hg19[,1], po=as.numeric(UCSC.hg19[,7]/2), name=gsub("chr","", UCSC.hg19[,1]));

pdffile  <- "advance_cytoBand.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");

circos(R=370, cir="hg19", type="label2", W=30, mapping=chr.l, side="out", cex=1.2, col=c("blue"));
circos(R=350, cir="hg19", W=20, mapping=hg19.b, type="arc2", col=col.c, lwd=8);
circos(R=270, cir="hg19", W=80, mapping=hg19.b, type="b3", B=T, col=col2[7]);
circos(R=190, cir="hg19", W=80, mapping=hg19.b, type="s",  B=T, col=col2, col.v=6, cex=abs(hg19.b[,6]));
circos(R=110, cir="hg19", W=80, mapping=hg19.b, type="b",  B=T, col=col2, col.v=6, cex=abs(hg19.b[,6]));

dev.off()
