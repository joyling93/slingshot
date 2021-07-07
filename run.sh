#!bin/bash
source /PROJ/home/xiezhuoming/miniconda2/bin/activate /PROJ/home/xiezhuoming/miniconda2/envs/Slingshot
Rscript slingshot.R \
--seurat_obj /PROJ/development/xiezhuoming/proj/velocyto/input/demo/seekone_demo/demo_seurat.rds \
--inference_cluster seurat_clusters \
--outdir .