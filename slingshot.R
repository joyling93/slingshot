library(Seurat)
library(slingshot,quietly = T)
#library(Matrix)
#slingPseudotime
#slingCurveWeights
rds.path <- '/PROJ/development/xiezhuoming/proj/velocyto/input/demo/seekone_demo/demo_seurat.rds'
sce <- readRDS(rds.path)

cl<-sce$seurat_clusters
rd<-Embeddings(sce, "umap")

sds <- slingshot(rd, clusterLabels = cl, reducedDim = 'umap',start.clus =5,end.clus =8)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sce)))
names(x = ident.colors) <- levels(x = sce)
cell.colors <- ident.colors[Idents(object = sce)]
names(x = cell.colors) <- colnames(x = sce)

png('curv.png')
plot(Embeddings(sce, reduction = 'umap'), col = cell.colors, asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black')
dev.off()

png('line.png')
plot(Embeddings(sce, reduction = 'umap'), col = cell.colors, asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black',type='l')
dev.off()