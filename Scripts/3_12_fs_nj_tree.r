
 RpairwiseFST.R
# aim: pairwise FST by population table and plot starting from vcf file
# needs: "SNPRelate" Rpackage, vcf input file, population metadata tab separated input file
# notes: Negative Fst are technical artifact of the computation (see Roesti el al. 2012) and are automatically replaced with zero.

##################################################

# closing previous genofile is just to assure a fresh start 
# if this give you an error it's ok, just means you already closed it
closefn.gds(genofile)
snpgdsClose(genofile)
library("SNPRelate")

#########
# Set input file names and parameters
#########
setwd("/home/chichi/data2/chichi/typha_lax/treemix_Fst/")

vcf.infile<-"/home/chichi/data2/chichi/typha_lax/treemix_Fst/merged.vcf.gz"
metadata.infile<- "/home/chichi/data2/chichi/typha_lax/treemix_Fst/pop_txt" # note: metadata file has to be a tab separated input text file with a "samples" column and a "pop" column
transformNegativeFSTtozero<- TRUE
keepLowerTriangularMatrixOnly<- TRUE

#########
# Parsing VCF file 
########

snpgdsVCF2GDS(vcf.infile, "ccm2.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm2.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

#########
# Parsing Population metadata 
#   needs a tab separated input text file,
#   with a "samples" column and a "pop" column
########

metadata=read.table(file = metadata.infile,header = T,sep = "\t",stringsAsFactors = F)
metadata=metadata[order(factor(metadata$samples, levels = sample.id)),]
pop_code=metadata$pop
poplevels=levels(as.factor(pop_code))

#####################################################
################## Pairwise FST  ####################
#####################################################



sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# pairwise populations matrix creation
res= outer(X= poplevels , Y= poplevels, 
           FUN = function(h,k){
             paste(h,k,sep = "/")
           }
)
colnames(res)=poplevels
rownames(res)=poplevels

as.data.frame(res)


# pairwise population matrix FST calculation

for(i in poplevels) {
  for(j in poplevels) {
    popelem= unlist(strsplit(res[i,j],"/"))
    
    #takes selection of samples an population to use for each pair
    flag<- pop_code %in% c(popelem[1],popelem[2])
    samp.sel<- sample.id[flag]
    pop.sel<- pop_code[flag]
    
    if (popelem[1]==popelem[2]){result="0"}else{
        result = snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), 
                         autosome.only=FALSE, method="W&C84")
        result = result$MeanFst
    }
    res[i,j]=as.character(result)
  }
}

########
# final edits and prints
########

# cell transformation from characters to numerics

res=data.frame(apply(res, 2, function(x) as.numeric(as.character(x))))
rownames(res)=poplevels

# trasforming negative FST into zero
# Negative Fst are technical artifact of the computation (see Roesti el al. 2012) and will be automatically replaced with zero inside this function.
if (transformNegativeFSTtozero==TRUE){
  for(i in poplevels) {
    for(j in poplevels) {
      if (res[i,j]<0){res[i,j]<-0} 
    }
  }
}

# keep only the lower triangular matrix
# (set the upper triangular to zero)
if(keepLowerTriangularMatrixOnly == TRUE){
  res[upper.tri(res)]<-0
}


########
# PRINT pairwise FST matrix to text file
########
timewithspaces= Sys.time()
timeAttr= gsub(pattern = " ",replacement = "_",x =  timewithspaces)
outfile <- paste("pairFstMatrix_", timeAttr, ".txt", sep="")
write.table(x = res, file = outfile, sep = "\t", dec = ".", 
            quote = F, row.names = T,col.names = NA)


# buid the NJ tree
# data
dist_matrix <- matrix(c(
  0, 0.0220845818406098, 0.0295167634493603, 0.017275170966017, 0.0906391542748899, 0.180873761016613,
  0.0220845818406098, 0, 0.0235132491599617, 0.0108138055111894, 0.0707225123213458, 0.146383299791471,
  0.0295167634493603, 0.0235132491599617, 0, 0.0194485323227546, 0.0864256934036283, 0.164764264405614,
  0.017275170966017, 0.0108138055111894, 0.0194485323227546, 0, 0.0688559944853521, 0.1329025916805,
  0.0906391542748899, 0.0707225123213458, 0.0864256934036283, 0.0688559944853521, 0, 0.201311091748804,
  0.180873761016613, 0.146383299791471, 0.164764264405614, 0.1329025916805, 0.201311091748804, 0
), nrow = 6, byrow = TRUE)

rownames(dist_matrix) <- colnames(dist_matrix) <- c("N1", "N2", "N3", "N4", "S", "T_angustifolia")

# Load required library

library(ape)
dist_obj <- as.dist(dist_matrix)
bionj_tree <- bionj(dist_obj)
min_positive <- min(bionj_tree$edge.length[bionj_tree$edge.length > 0])
bionj_tree$edge.length[bionj_tree$edge.length < 0] <- min_positive 
bionj_tree_rooted <- root(bionj_tree, outgroup = "T_angustifolia", resolve.root = TRUE)

distance_bootstrap <- function(dist_matrix, method = "bionj", B = 100) {
  n_taxa <- nrow(dist_matrix)
  trees <- vector("list", B)
  
  # get original upper triangular values
  tri_ind <- upper.tri(dist_matrix)
  orig_values <- dist_matrix[tri_ind]
  
  for (i in 1:B) {
    # resample the distance values with replacement
    boot_values <- sample(orig_values, replace = TRUE)
    boot_matrix <- matrix(0, n_taxa, n_taxa)
    boot_matrix[tri_ind] <- boot_values
    boot_matrix <- boot_matrix + t(boot_matrix)
    diag(boot_matrix) <- 0
    colnames(boot_matrix) <- rownames(boot_matrix) <- colnames(dist_matrix)
    
    # build tree
    if (method == "bionj") {
      boot_tree <- bionj(as.dist(boot_matrix))
      # deal with negative branch lengths
      boot_tree$edge.length[boot_tree$edge.length < 0] <- min_positive / 10
    } else {
      boot_tree <- fastme.bal(as.dist(boot_matrix))
      boot_tree$edge.length[boot_tree$edge.length < 0] <- min_positive / 10
    }
    
    # reroot the tree
    boot_tree <- root(boot_tree, outgroup = "T_angustifolia", resolve.root = TRUE)
    trees[[i]] <- boot_tree
  }
  return(trees)
}
set.seed(123)
bionj_boot_trees <- distance_bootstrap(dist_matrix, method = "bionj", B = 1000)
# Calculate bootstrap support values
bionj_support <- prop.clades(bionj_tree_rooted, bionj_boot_trees)
bionj_tree_rooted$node.label <- bionj_support
# Plot the BIONJ tree with bootstrap support values
plot(bionj_tree_rooted, type = "phylogram", 
     main = "BIONJ Tree with Bootstrap Support",
     edge.width = 2, cex = 1.0)
nodelabels(bionj_support, adj = c(1.2, -0.3), frame = "none", cex = 0.8, bg = "white")
add.scale.bar()

# save the tree to a text file
write.tree(bionj_tree_rooted, file = "bionj_tree_with_bootstrap.txt")

