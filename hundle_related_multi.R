setwd("~/git/IBDplot/")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt")

library(gdsfmt)
library(SNPRelate)
# UNR
#######################
locus_number = 10000
error_rate = 0.03
REP = 1000
K0_tmp = K1_tmp = rep(0,REP)
for (i in 1: REP){
  cat("rep =",i,"\n")
  INLINE1 = "gcc -DHAVE_INLINE -lgsl -lm -lgslcblas  -Wall -o kinfer_unr kinfer_unr.c"
  system(INLINE1)
  INLINE2 = paste0("./kinfer_unr ",locus_number," ",error_rate)
  system(INLINE2)
  tmp = read.delim("unr.txt",header = F)
  
  genotype = cbind(tmp$V1,tmp$V2)
  snp.id = seq(1, locus_number)
  sample.id = c("1", "2")
  snp.position = rep(0, locus_number)
  snp.chromosome = rep(1, locus_number)
  
  #file.remove("hsp.gds")
  snpgdsCreateGeno(gds.fn = "unr.gds", genmat = genotype, sample.id=sample.id, snp.id=snp.id, snp.rs.id=NULL,
                   snp.chromosome=snp.chromosome, snp.position=snp.position, snp.allele=NULL, snpfirstdim=TRUE,
                   compress.annotation="ZIP_RA.max", compress.geno="", other.vars=NULL)
  
  
  genofile <- snpgdsOpen("~/git/IBDplot/unr.gds")
  tmp = read.delim("AF.txt",header = F)
  AF = tmp$V1
  ibd <- snpgdsIBDMLE(genofile,sample.id=sample.id, snp.id=snp.id, missing.rate=0.0, allele.freq = AF, num.thread = 12)
  
  ibd.coeff <- snpgdsIBDSelection(ibd)
  K0_tmp[i] = ibd.coeff$k0
  K1_tmp[i] = ibd.coeff$k1
  
  snpgdsClose(genofile)
  
}
save(list=ls(),file = paste0("res_unr.Rdata"))
par(pty="s")
plot(NULL, xlim = c(0,1), ylim = c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col=1, lty=2)
points(K0_tmp, K1_tmp, xlim=c(0,1), ylim=c(0,1), pch = 20, col="black")
title(main = paste0("locus_number = ",locus_number,"\n error_rate = ",error_rate))
#hist(K1_tmp)
#######################
# PO
#######################
REP = 100
K0_tmp = K1_tmp = rep(0,REP)
for (i in 1: REP){
  cat("rep =",i,"\n")
  INLINE1 = "gcc -DHAVE_INLINE -lgsl -lm -lgslcblas  -Wall -o kinfer_po kinfer_po.c"
  system(INLINE1)
  INLINE2 = paste0("./kinfer_po ",locus_number," ",error_rate)
  system(INLINE2)
  tmp = read.delim("po.txt",header = F)
  
  genotype = cbind(tmp$V1,tmp$V2)
  snp.id = seq(1, locus_number)
  sample.id = c("1", "2")
  snp.position = rep(0, locus_number)
  snp.chromosome = rep(1, locus_number)
  
  #file.remove("hsp.gds")
  snpgdsCreateGeno(gds.fn = "po.gds", genmat = genotype, sample.id=sample.id, snp.id=snp.id, snp.rs.id=NULL,
                   snp.chromosome=snp.chromosome, snp.position=snp.position, snp.allele=NULL, snpfirstdim=TRUE,
                   compress.annotation="ZIP_RA.max", compress.geno="", other.vars=NULL)
  
  
  genofile <- snpgdsOpen("~/git/IBDplot/po.gds")
  tmp = read.delim("AF.txt",header = F)
  AF = tmp$V1
  ibd <- snpgdsIBDMLE(genofile,sample.id=sample.id, snp.id=snp.id, missing.rate=0.0, allele.freq = AF, num.thread = 12)
  
  ibd.coeff <- snpgdsIBDSelection(ibd)
  K0_tmp[i] = ibd.coeff$k0
  K1_tmp[i] = ibd.coeff$k1
  
  snpgdsClose(genofile)
  
}
save(list=ls(),file = paste0("res_po.Rdata"))
#par(pty="s")
#plot(NULL, xlim = c(0,1), ylim = c(0,1), xlab="k0", ylab="k1")
#lines(c(0,1), c(1,0), col=1, lty=2)
points(K0_tmp, K1_tmp, xlim=c(0,1), ylim=c(0,1), pch = 20, col="red")
#title(main = paste0("locus_number = ",locus_number,"\n error_rate = ",error_rate))
#hist(K1_tmp)
#######################
# HS
#######################
REP = 100
K0_tmp = K1_tmp = rep(0,REP)
for (i in 1: REP){
  cat("rep =",i,"\n")
  INLINE1 = "gcc -DHAVE_INLINE -lgsl -lm -lgslcblas  -Wall -o kinfer_hsp kinfer_hsp.c"
  system(INLINE1)
  INLINE2 = paste0("./kinfer_hsp ",locus_number," ",error_rate)
  system(INLINE2)
  tmp = read.delim("hsp.txt",header = F)
  
  genotype = cbind(tmp$V1,tmp$V2)
  snp.id = seq(1, locus_number)
  sample.id = c("1", "2")
  snp.position = rep(0, locus_number)
  snp.chromosome = rep(1, locus_number)
  
  #file.remove("hsp.gds")
  snpgdsCreateGeno(gds.fn = "hsp.gds", genmat = genotype, sample.id=sample.id, snp.id=snp.id, snp.rs.id=NULL,
                   snp.chromosome=snp.chromosome, snp.position=snp.position, snp.allele=NULL, snpfirstdim=TRUE,
                   compress.annotation="ZIP_RA.max", compress.geno="", other.vars=NULL)
  
  
  genofile <- snpgdsOpen("~/git/IBDplot/hsp.gds")
  tmp = read.delim("AF.txt",header = F)
  AF = tmp$V1
  ibd <- snpgdsIBDMLE(genofile,sample.id=sample.id, snp.id=snp.id, missing.rate=0.0, allele.freq = AF, num.thread = 12)
  
  ibd.coeff <- snpgdsIBDSelection(ibd)
  K0_tmp[i] = ibd.coeff$k0
  K1_tmp[i] = ibd.coeff$k1
  
  snpgdsClose(genofile)
  
}
save(list=ls(),file = paste0("res_hsp.Rdata"))
#par(pty="s")
#plot(NULL, xlim = c(0,1), ylim = c(0,1), xlab="k0", ylab="k1")
#lines(c(0,1), c(1,0), col=1, lty=2)
points(K0_tmp, K1_tmp, xlim=c(0,1), ylim=c(0,1), pch = 20, col="blue")
#title(main = paste0("locus_number = ",locus_number,"\n error_rate = ",error_rate))
#hist(K1_tmp)
#######################