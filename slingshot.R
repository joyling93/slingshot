library(Seurat)
library(slingshot,quietly = T)
library(argparse)
#library(Matrix)
#slingPseudotime
#slingCurveWeights


parser = ArgumentParser()
parser$add_argument("--seuart_obj", help="seurat 对象的RDS文件。",required=TRUE)
parser$add_argument("--inference_clusters)", help="设置用于轨迹推断的分组变量，默认是seurat_clusters"
                    ,required=TRUE,default='seurat_clusters')
parser$add_argument("--reduction_space", help="设置需要推断轨迹的降维空间，默认是umap",default='umap')
parser$add_argument("--start_cluster", help="设置轨迹推断起始 cluster")
parser$add_argument("--end_cluster", help="设置轨迹推断结束 cluster")
parser$add_argument("--outdir", help='默认是./output',required=TRUE,default="./output")


args <- parser$parse_args()
str(args)

seuart_obj <- args$seuart_obj
inference_clusters <- args$inference_clusters
reduction_space <- args$reduction_space
start_cluster <- args$start_cluster
end_cluster <- args$end_cluster
outdir <- args$outdir
#rds.path <- '/PROJ/development/xiezhuoming/proj/velocyto/input/demo/seekone_demo/demo_seurat.rds'


sce <- readRDS(seurat_onj)

cl<-sce[[inference_clusters]][,1]
rd<-Embeddings(sce, reduction_space)

sds <- slingshot(rd, clusterLabels = cl, reducedDim = 'umap',start.clus =5,end.clus =8)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sce)))
names(x = ident.colors) <- levels(x = sce)
cell.colors <- ident.colors[Idents(object = sce)]
names(x = cell.colors) <- colnames(x = sce)

png(file.path(outdir,'curv.png'))
plot(rd, col = cell.colors, asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black')
dev.off()

png(file.path(outdir,'line.png'))
plot(rd, col = cell.colors, asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black',type='l')
dev.off()