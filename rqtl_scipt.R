library(qtl)
#Step 1: import your the data from the working directory
#step 2: fill your directory in the dir =""
#step 3: execute the block, you can change the jittermap as per required, here i used default
sug <- read.cross("csv", "C:/Users/Honey/Desktop/ASSIGN", "sug1.csv", genotypes=c("A", "B", "X"), alleles=c("A", "B"))
sug2 <- jittermap(object = sug, amount = 1e-6)
plot(sug2)
#step 4: single qtl analysis
sug2 <- calc.genoprob(sug2, step=1)
out.em <- scanone(sug2)
summary(out.em)
summary(out.em, threshold=3)
plot(out.em)
#step 5: genome scan via haley- knott regression
out.hk <- scanone(sug2, method="hk")
plot(out.em, out.hk, col=c("blue", "red"))
#step 6: LOD curves for just 7 and 19 chromosomes
plot(out.em, out.hk, col=c("blue", "red"), chr=c(7,19))
#step 7: plot the differences
plot(out.hk - out.em, ylim=c(-0.3, 0.3), ylab="LOD(HK)-LOD(EM)")
#step 8: permutation tests
operm <- scanone(sug2, method="hk", n.perm=1000)
plot(operm)
#step 9: interval estimates of qtl location 7
lodint(out.hk, chr=7, drop=2) 
bayesint(out.hk, chr=7, prob=0.99)
#step 10: QTL effects
max(out.hk) 
mar <- find.marker(sug2, chr=7, pos=0) 
plotPXG(sug2, marker=mar)
#step 11: 2d qtl scan
sug2 <- calc.genoprob(sug2, step=2)
out2 <- scantwo(sug2, method="hk")
plot(out2)
#step 12: multiple QTL analysis
sug2 <- calc.genoprob(sug2, step=1)
qtl <- makeqtl(sug2, chr=c(7,19), pos=c(0, 12), what="prob")
out.fq <- fitqtl(sug2, qtl=qtl, method="hk") 
rqtl <- refineqtl(sug2, qtl=qtl, method= "hk")
plotLodProfile(rqtl)
