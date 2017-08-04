rm(list=ls());
dev.off()

## setting working directory
# setwd("/Users/yourname/R");
library(OmicCircos);
options(stringsAsFactors = FALSE);

set.seed(1234);

## EColi_K-12_GM4792
inf <- read.table("EColi_db.txt", header=T);
db  <- segAnglePo(inf, seg=inf$seg.name);

## genes on EColi_K-12_GM4792
dat <- read.table("EColi_K-12_GM4792.txt", header=T);

## simulation data
s1  <- sample(1:nrow(dat), 500);
s2  <- sample(1:nrow(dat), 400);
s3  <- sample(1:nrow(dat), 300);
s4  <- sample(1:nrow(dat), 200);
d1  <- data.frame(chr=dat[s4,1], po=dat[s4,2], val=rnorm(200));
d2  <- data.frame(chr=dat[s3,1], po=dat[s3,2], val=rnorm(300));
lab <- sample(1:nrow(dat), 3);

## colors
cols <- rainbow(10, alpha=0.8);

pdf("advance_user_genome.pdf", 8, 8)
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="advance_user_genome.R");

## chromosome
circos(R=300, type="chr", cir=db, W=30, scale=T, print.chr.lab=F);

## simulated gene positions
circos(R=270, cir=db, W=30, mapping=dat[s1,], type="b3", B=T, col=cols, lwd=1);
circos(R=240, cir=db, W=30, mapping=dat[s2,], type="b3", B=F, col=cols, lwd=1);
circos(R=210, cir=db, W=30, mapping=dat[s3,], type="b3", B=T, col=cols, lwd=1);

## simulated gene expressions 
circos(R=160, cir=db, W=50, mapping=d1, col.v=3, type="b2", B=F, 
       col=cols[c(7,9)], cutoff=0, lwd=2, scale=T);
circos(R=110, cir=db, W=50, mapping=d2, col.v=3, type="b2", B=T, 
       col=cols[c(6,1)], cutoff=0, lwd=2, scale=T);

## gene labels
circos(R=100, type="chr", cir=db, W=30, scale=F, print.chr.lab=F);
circos(R=50, cir=db, W=20, mapping=dat[lab,],  type="label2", col.v=5, col="blue", cex=0.8);

legend(160, 50, paste0("n", 0:4), col=rainbow(10)[c(1:5)], lwd=4, bty="n", horiz=TRUE, cex=1);
legend(160, 20, paste0("n", 5:9), col=rainbow(10)[c(6:10)], lwd=4, bty="n", horiz=TRUE, cex=1);

dev.off()
