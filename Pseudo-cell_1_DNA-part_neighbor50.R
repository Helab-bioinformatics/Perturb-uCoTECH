# ===========================================================================================================
# R codes of Pseudo-cell calculation method of perturbation effects, developed by Chen Li in August, 2023.
# ===========================================================================================================

setwd("/direction/to/your/triomics/rds")
library(Seurat)
library(dplyr)

# load rds files ------------------------------------------------------------------------------------------------------
triomics.rds <- readRDS("./your_file.rds")
triomics.rds.perturb <- subset(triomics.rds, group %in% c("perturb") & cellpopulation %in% c("HSPC") & nFeature_sgRNA>0)

# load the Cicero pair information ------------------------------------------------------------------------------------
pair.cicero <- read.csv("./Cicero_pairs.csv", header = T)
pair.near <- read.csv("./Nearest_pairs.csv", header = T)
common <- intersect(pair.cicero$enhancer, pair.near$enhancer)

enh.add <- setdiff(rownames(triomics.rds.perturb@assays$sgRNAonCRE@data), unique(pair.cicero$enhancer))
pair.cicero.add <- pair.near[pair.near$enhancer %in% enh.add, ]
pair.cicero <- pair.cicero[,-3]
pair.cicero.add <- pair.cicero.add[,-3]
pair.cicero.final <- rbind(pair.cicero, pair.cicero.add)

gene <- unique(pair.cicero.final$gene)
enh <- unique(pair.cicero.final$enhancer)

setdiff(unique(pair.cicero.final$gene), rownames(triomics.rds.perturb@assays$RNA@data))
setdiff(rownames(triomics.rds.perturb@assays$RNA@data), unique(pair.cicero.final$gene))
intersect(pair.cicero.final$gene, rownames(triomics.rds.perturb@assays$RNA@data))

pair.cicero.final.sub <- pair.cicero.final[pair.cicero.final$gene %in% rownames(triomics.rds.perturb@assays$RNA@data), ]
pair <- pair.cicero.final.sub
pair <- unique(pair)
pair <- pair[pair$enhancer %in% rownames(triomics.rds.perturb@assays$sgRNAonCRE@counts),]
rownames(pair) <- paste(pair$enhancer, pair$gene, sep = '.')

# calculate the nearest 50 neighbors ------------------
dt <- dist(triomics.rds.perturb@reductions$umap@cell.embeddings) %>% as.matrix()
cn <- colnames(dt)
cn <- apply(dt, MARGIN = 1, FUN = function(x) cn[order(x)[1:50]])
dim(cn)

# get sgRNAonCRE counts, cell names and feature names -----
sgRNAonCRE.mat <- triomics.rds.perturb@assays$sgRNAonCRE@counts

# add row and update row names ------------------------
sgRNAonCRE.mat <- sgRNAonCRE.mat[pair$enhancer,]
rownames(sgRNAonCRE.mat) <- rownames(pair)

cell.name <- colnames(triomics.rds.perturb@assays$sgRNAonCRE@counts)
feature.name <- rownames(sgRNAonCRE.mat)
feature.peak <- rownames(triomics.rds.perturb@assays$DNA@data)

sum(colnames(cn) != cell.name)

# calculate H3K27ac signals in enhancers --------------
num.feat <- 1
DNA.mean <- apply(sgRNAonCRE.mat, MARGIN = 1, FUN = function(x){
  num.cell <- 1
  tmp <- sapply(x, FUN = function(y) {
    if(y!=0){
      neighbor <- cn[, num.cell]
      neighbor0 <- neighbor[x[match(neighbor, cell.name)]==0]
      neighbor1 <- neighbor[x[match(neighbor, cell.name)]!=0]
      if(exists("neighbor1") & pair[feature.name[num.feat], "enhancer"] %in% feature.peak) {trueM <- triomics.rds.perturb@assays$DNA@data[pair[feature.name[num.feat],'enhancer'], neighbor1] %>% mean()} else {trueM <-0}
      if(exists("neighbor0") & pair[feature.name[num.feat], "enhancer"] %in% feature.peak) {falseM <- triomics.rds.perturb@assays$DNA@data[pair[feature.name[num.feat], "enhancer"], neighbor0] %>% mean()} else {falseM <-0}
      k27 <- data.frame(trueM = trueM, falseM = falseM)
    }
    num.cell <<- num.cell + 1
    if(exists("k27")) {return(k27)}
  })
  tmp <- tmp[!sapply(tmp, is.null)]
  num.feat <<- num.feat +1
  return(tmp)
})

DNA.mean <- DNA.mean[!sapply(DNA.mean, FUN = function(x){length(x)==0})]
summary <- sapply(DNA.mean, FUN = function(x){
  DNA.mean.unls <- unlist(x)
  k27.enh.df <- data.frame(trueM = DNA.mean.unls[seq(1,length(DNA.mean.unls), 2)], falseM=DNA.mean.unls[seq(2,length(DNA.mean.unls), 2)])
  k27.enh.df[is.na(k27.enh.df)]<-0
  k27.enh.df.rmT0 <- k27.enh.df[k27.enh.df$trueM!=0,]
  k27.enh.df.rmT0 <- k27.enh.df.rmT0[k27.enh.df.rmT0$falseM!=0,]
  mean <- mean(k27.enh.df.rmT0$trueM/k27.enh.df.rmT0$falseM)
  if (nrow(k27.enh.df.rmT0) >1) {
    wilcox_test <- wilcox.test(k27.enh.df.rmT0$trueM, k27.enh.df.rmT0$falseM, alternative = "greater",paired=TRUE, exact = F)
    w.p.value <- wilcox_test$p.value
  }
  else {
    w.p.value <- 1
  }
  data.frame(mean = mean, 
             w.p.value = w.p.value,
             trueM = mean(k27.enh.df.rmT0$trueM),
             falseM = mean(k27.enh.df.rmT0$falseM))
})
summary <- t(summary)
write.csv(summary, "Pseudo-cell_DNA-part.csv", quote = F)

