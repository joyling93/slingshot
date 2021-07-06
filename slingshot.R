library(Seurat)
library(slingshot)
library(Matrix)
sce <- readRDS()

cl<-sce$seurat_clusters
rd<-Embeddings(sce, "umap")

sds <- slingshot(rd, clusterLabels = cl)

png('line.png')
plot(SlingshotDataSet(sds), col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black')
dev.off()

png('curv.png')
lines(SlingshotDataSet(sds), lwd=2, col='black',type='l')
dev.off()