软件简介
Slingshot 旨在模拟单细胞 RNA 测序数据中的发育轨迹，并可以处理任意多个分支事件。

简要使用方法
修改seurat_obj、inference_cluster、outdir 三个参数，运行 run.sh。
Slingshot 在运行时需要提供一个分组变量（即 inference_cluster ），但不必须提供轨迹推断的起始。

参数说明
--seurat_obj	经过降维的 seurat 对象，RDS文件。
--inference_clusters	设置用于轨迹推断的分组变量，默认是seurat_clusters。
--reduction_space	设置需要推断轨迹的降维空间，默认是umap。
--start_cluster	设置轨迹推断起始 cluster，非必须。
--end_cluster	设置轨迹推断结束 cluster，非必须。
--outdi	默认是./output。

结果文件说明
slingshot_obj.rds 	slingshot 对象，包含轨迹推断结果和画图坐标，可用函数 slingPseudotime 提取 Pseudotime value。
curv.png	轨迹推断图，拟合曲线，更改轨迹推断开始和结束 cluster 可以修正拟合结果。
line.png	轨迹推断图，折线。