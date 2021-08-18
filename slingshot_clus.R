library(Seurat)
library(slingshot,quietly = T)
library(argparse)
#library(Matrix)
#slingPseudotime
#slingCurveWeights


parser = ArgumentParser()
parser$add_argument("--seurat_obj", help="seurat 对象的RDS文件。",required=TRUE)
parser$add_argument("--clus)", help="subclus file"
                    ,required=TRUE)
parser$add_argument("--group)", help="plot group"
                    ,required=TRUE,default='Sample')
parser$add_argument("--reduction_space", help="设置需要推断轨迹的降维空间，默认是umap",default='umap')
parser$add_argument("--start_cluster", help="设置轨迹推断起始 cluster")
parser$add_argument("--end_cluster", help="设置轨迹推断结束 cluster")
parser$add_argument("--outdir", help='默认是./output',required=TRUE,default="./output")


args <- parser$parse_args()
str(args)

seurat_obj <- args$seurat_obj
reduction_space <- args$reduction_space
start_cluster <- args$start_cluster
end_cluster <- args$end_cluster
outdir <- args$outdir
group<- args$group
#rds.path <- '/PROJ/development/xiezhuoming/proj/velocyto/input/demo/seekone_demo/demo_seurat.rds'

clus<-read.table(args$clus,header=T,sep=',',row.names=1)
obj <- readRDS(seurat_obj)
DefaultAssay(obj) <- "RNA"
obj<-subset(obj,cells=c(rownames(clus)))

sce <- 
        ScaleData(obj, verbose = FALSE) %>% 
        RunPCA(npcs = 30, verbose = FALSE) %>% 
        RunUMAP( reduction = "pca", dims = 1:30) %>% 
        RunTSNE( reduction = "pca", dims = 1:30) %>% 
        FindNeighbors( dims = 1:30) %>% 
        FindClusters(resolution = 0.6)

cl<-sce$seurat_clusters
rd<-Embeddings(sce, reduction_space)

sds <- slingshot(rd, clusterLabels = cl, reducedDim = reduction_space,start.clus =start_cluster,end.clus =end_cluster)
saveRDS(sds,file.path(outdir,'slingshot_obj.rds'))

group<-as.factor(clus[[group]])
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = group)))
names(x = ident.colors) <- levels(x = group)
cell.colors <- ident.colors[group]
names(x = cell.colors) <- group

png(file.path(outdir,'curv.png'))
plot(rd, col = cell.colors, asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black')
legend("topright",legend=names(ident.colors),col=ident.colors,cex=0.5,pch = 16)
dev.off()

png(file.path(outdir,'line.png'))
plot(rd, col = cell.colors, asp = 1, pch = 16)
lines(SlingshotDataSet(sds), lwd=2, col='black',type='l')
legend("topright",legend=names(ident.colors),col=ident.colors,cex=0.5,pch = 16)
dev.off()


