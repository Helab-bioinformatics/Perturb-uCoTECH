# ===================================================================================================================
# R codes of pseudocell grouping to calculate perturbation effects, developed and updated by Chen Li in June, 2024.
# ===================================================================================================================

setwd("/path/to/your/file")
library(Seurat)
library(dplyr)

# --- (1) load the files ----------------------------------------------
yourRds <- readRDS("your-Seurat-object.rds")
yourRds.perturb <- subset(yourRds, cellType %in% c("perturbed"))
DefaultAssay(yourRds.perturb) <- "RNA"

# --- (2) find neighboring 50 cells -----------------------------------
dt <- dist(yourRds.perturb@reductions$umap@cell.embeddings) %>% as.matrix()
cn <- colnames(dt)
cn <- apply(dt, MARGIN = 1, FUN = function(x) cn[order(x)[1:50]])
dim(cn)

sgrEnh.mat.raw <- yourRds.perturb@assays$sgrEnh@counts
sgrEnh.mat <- sgrEnh.mat.raw[rownames(sgrEnh.mat.raw) %in% rownames(yourRds.perturb@assays$DNA@data),]
cell.name <- colnames(yourRds.perturb@assays$sgrEnh@counts) 
feature.name <- rownames(yourRds.perturb@assays$sgrEnh@counts) 
feature.name <- feature.name[feature.name %in% rownames(yourRds.perturb@assays$DNA@data)] 
sum(colnames(cn) != cell.name) 

num.feat <- 1
neighbors <- apply(sgrEnh.mat, MARGIN = 1, FUN = function(x){
  num.cell <- 1
  tmp <- sapply(x, FUN = function(y) {
    if(y!=0){
      neighbor <- cn[, num.cell]
      neighbor0 <- neighbor[x[match(neighbor, cell.name)]==0]
      neighbor1 <- neighbor[x[match(neighbor, cell.name)]!=0]
      if(exists("neighbor1")) {
        trueM <- yourRds.perturb@assays$DNA@data[feature.name[num.feat], neighbor1] %>% mean()
        num.neighbor1 <- length(neighbor1)
      } else {trueM <-0; num.neighbor1 <- 0}
      if(exists("neighbor0")) {
        falseM <- yourRds.perturb@assays$DNA@data[feature.name[num.feat], neighbor0] %>% mean()
        num.neighbor0 <- length(neighbor0)
      } else {falseM <-0; num.neighbor0 <- 0}
      k27 <- data.frame(trueM = trueM, falseM = falseM, num.neighbor1 = num.neighbor1, num.neighbor0 = num.neighbor0)
    }
    num.cell <<- num.cell + 1
    if(exists("k27")) {return(k27)}
  })
  tmp <- tmp[!sapply(tmp, is.null)]
  num.feat <<- num.feat +1
  return(tmp)
})


