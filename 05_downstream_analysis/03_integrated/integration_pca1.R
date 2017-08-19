rm(list=ls());

library(OmicCircos);
options(stringsAsFactors = FALSE);

data("TCGA.BC.fus");
data("TCGA.BC.cnv.2k.60");
data("TCGA.BC.gene.exp.2k.60");
data("TCGA.BC.sample60");

## gene expression data for PCA
exp.m <- TCGA.BC.gene.exp.2k.60[,c(4:ncol(TCGA.BC.gene.exp.2k.60))];

cnv   <- TCGA.BC.cnv.2k.60;

## PCA
type.n  <- unique(TCGA.BC.sample60[,2]);
colors  <- rainbow(length(type.n), alpha=0.5);

pca.col <- rep(NA, nrow(TCGA.BC.sample60));
for (i in 1:length(type.n)){
  n   <- type.n[i];
  n.i <- which(TCGA.BC.sample60[,2] == n);
  n.n <- TCGA.BC.sample60[n.i,1];
  g.i <- which(colnames(exp.m) %in% n.n);
  pca.col[g.i] <- colors[i];
}

exp.m   <- na.omit(exp.m);
pca.out <- prcomp(t(exp.m), scale = TRUE);

## subtype cnv
cnv.i <- c();
for (i in 1:length(type.n)){
  n     <- type.n[i];
  n.i   <- which(TCGA.BC.sample60[,2] == n);
  n.n   <- TCGA.BC.sample60[n.i,1];
  cnv.i <- which(colnames(cnv) %in% n.n);
}

## main
pdf("integration_pca1.pdf", 8,8);
par(mar=c(5, 5, 5, 5));

plot(pca.out$x[,1]*10, pca.out$x[,2]*10, pch=19, col=pca.col, main="integration_pca1.R",  
     cex=2, xlab="PC1", ylab="PC2", ylim=c(-200, 460), xlim=c(-200,460));
legend(200,0, c("Basal","Her2","LumA","LumB"), pch=19, col=colors[c(2,4,1,3)], cex=1, 
       title ="Gene Expression (PCA)", box.col="white");

circos(xc=260, yc=260, R=200, cir="hg18", W=4, type="chr", print.chr.lab=T);
R.v <- 160;
for (i in 1:length(type.n)){
  n     <- type.n[i];
  n.i   <- which(TCGA.BC.sample60[,2] == n);
  n.n   <- TCGA.BC.sample60[n.i,1];
  cnv.i <- which(colnames(cnv) %in% n.n);
  cnv.v <- cnv[,cnv.i];
  cnv.v[cnv.v > 2]  <- 2;
  cnv.v[cnv.v < -2] <- -2;
  cnv.m <- cbind(cnv[,c(1:3)], cnv.v);
  if (i %% 2 == 1){
    circos(xc=260, yc=260, R=R.v, cir="hg18", W=30, mapping=cnv.m, col.v=4,  
           type="ml3", B=T, lwd=0.5, cutoff=0);
  } else {
    circos(xc=260, yc=260, R=R.v, cir="hg18", W=30, mapping=cnv.m, col.v=4,  
           type="ml3", B=F, lwd=0.5, cutoff=0);
  }
  
  R.v <- R.v - 30;
}

legend(-140, 460, c("1 Basal", "2 Her2", "3 LumA", "4 LumB", "(center)"), cex=1, 
       title ="CNV (OmicCircos)", box.col="white");

dev.off()



