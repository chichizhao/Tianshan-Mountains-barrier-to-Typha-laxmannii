#!/r
# tutorial: http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
# set the working directory
setwd("/home/chichi/data2/chichi/typha_lax/pca")
rm(list=ls())

# load the packages
library(gdsfmt)
library(SNPRelate)
vcf.fn <- "/home/chichi/data2/chichi/typha_lax/pca/126_raw_filter_variants_maf0.01_miss0.9.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "pac_pruned_126_5.gds", method="biallelic.only")
genofile <- snpgdsOpen("pac_pruned_126_5.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
str(snpset)
names(snpset)
snpset.id <- unlist(unname(snpset))
snpset.id
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=16)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],
    EV2 = pca$eigenvect[,2],
    stringsAsFactors = FALSE)
head(tab)



# plot with base graphics
png("pac_pruned_126_test.png", width=1800, height=1800, res=300)
plot(tab$EV1, tab$EV2, xlab="PC1", ylab="PC2", pch=20, cex=2,xlim=c(-0.3,0.1), ylim=c(-0.5,1))
text(tab$EV1, tab$EV2, labels=tab$sample.id, cex=0.8, pos=4)
dev.off()
# base on the west middle and east geographical location of the samples, we can color the samples
# get the sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.table("sample_pop.txt", header=FALSE, sep="\t")
pop_code <- pop_code$V2
head(cbind(sample.id, pop_code))
tab <- data.frame(sample.id = sample.id,
    EV1 = pca$eigenvect[,1],
    EV2 = pca$eigenvect[,2],
    pop = factor(pop_code)[match(sample.id, sample.id)],
    stringsAsFactors = FALSE)
head(tab)

png("pac_126.png", width=3000, height=3000, res=600)
par(family="Arial", cex=1)
# here we want to set the color of groups in specific color
colors <- c("#2c71bfc6", "#69b039c2")
# Create the plot
plot(tab$EV1, tab$EV2, xlab="PC1 2.38%", ylab="PC2 1.89%", 
     pch=20, cex=1.5, col=colors[tab$pop], 
     xlim=c(-0.25, 0.1), ylim=c(-0.12, 0.45), 
     cex.lab=1.2, cex.axis=1.2)  # Increase label and axis text size

legend("topleft", legend=c("North", "South"), col=colors[1:length(levels(tab$pop))], pch=20, cex=2)
dev.off()

