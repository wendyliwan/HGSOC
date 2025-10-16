setwd("F:/毕设/data3/finial/F5--免疫细胞/新的immune")
table(ovar3@meta.data$celltype.main)
immune_cell.ids1 <- which(ovar3@meta.data$celltype.main.22 == c("CD8+T","CD4+T"))
immune_cell.ids2 <- which(ovar3@meta.data$celltype.main.22 == "B_Plasma")
immune_cell.ids3 <- which(ovar3@meta.data$celltype.main.22 == c("Myeloid","NK","Tregs"))

immune_cells <- subset(ovar3, cells = c(immune_cell.ids1,immune_cell.ids2,immune_cell.ids3))

immuDimPlot(immune_cells,reduction = "umap",shuffle = F,label=F,group.by = "tissue_type",cols = c("#db5a6b","#0077ce"))
ggsave("immune_cells_umap.pdf",width=7,height=6)
DimPlot(immune_cells,reduction = "tsne",label=F,shuffle = F,group.by = "tissue_type",cols = c("#db5a6b","#0077ce"))
ggsave("immune_cells_tsne.pdf",width=7,height=6)    
                        
Immu_sce=CreateSeuratObject(counts = GetAssayData(immune_cells,assay = "RNA",layer = 'counts'),
                            meta.data = immune_cells@meta.data)
Immu_sce=NormalizeData(Immu_sce)%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA(verbose = F)
Immu_sce=RunHarmony(Immu_sce,group.by.var="orig.ident",assay.use="RNA",max.iter.harmony=20)
ElbowPlot(Immu_sce,ndims = 50,reduction = "pca")
Immu_sce=FindNeighbors(Immu_sce,reduction = "pca",dims = 1:20)
Immu_sce=FindClusters(Immu_sce,resolution = 0.5)
# setwd("F:/毕设/data3/immune cell/IMMU")
setwd("F:/毕设/data3/finial/F5--免疫细胞/新的immune")
table(Immu_sce@meta.data$seurat_clusters)
#    0    1    2    3    4    5    6    7    8    9   10 
# 1360 1172 1108  791  754  742  346  278   68   36   19 
#    0    1    2    3    4    5    6    7    8    9   10   11   12 
# 1169 1126  984  794  655  599  452  447  371  280  151  101   68
Immu_sce=RunUMAP(Immu_sce,reduction = "pca",dims = 1:20,min.dist = 0.85)%>%RunTSNE(dims=1:20,reduction = "pca")
col=pal_futurama()(13)
DimPlot(Immu_sce,reduction = "umap",label = T,shuffle = F,raster = F,cols = col)
ggsave(file="Immu_umap.pdf",width = 6,height = 6)
DimPlot(Immu_sce,reduction = "tsne",label = T,shuffle = F,raster = F)
ggsave(file="Immu_tsne.pdf",width = 6,height = 6)

DimPlot(Immu_sce,reduction = "umap",group.by = 'orig.ident',shuffle = F,raster = F,cols = Biocols)
#DimPlot(Immu_sce,reduction = "tsne",group.by = 'orig.ident',raster = F,cols = Biocols)
ggsave(file="Immu_orig_umap.pdf",width = 8,height = 6)

DimPlot(Immu_sce,reduction = "umap",split.by = "tissue_type",label = T,raster = F,cols = Biocols)
ggsave(file="Immu_tissue-map.pdf",width = 10,height = 6)
DimPlot(Immu_sce,reduction = "tsne",split.by = "tissue_type",label = T,raster = F,cols = Biocols)
ggsave(file="Immu_tissue-tsne.pdf",width = 10,height = 6)

DimPlot(Immu_sce,reduction = "umap",group.by = "tissue_type",shuffle = T,raster = F,cols = c("#db5a6b","#0077ce"))
ggsave(file="tissue_umap.pdf",width = 6,height = 6)
DimPlot(Immu_sce,reduction = "tsne",group.by = "tissue_type",shuffle = T,raster = F,cols = Biocols)
ggsave(file="tissue_tsne.pdf",width = 6,height = 6)

Immu_markers=FindAllMarkers(Immu_sce,test.use = "wilcox",
                            only.pos = T,
                            logfc.threshold = 0.25,
                            min.pct = 0.25)
write.csv(Immu_markers,file = "Immu_markers.csv")
Immu_markers=read.csv("Immu_markers.csv")

Immu.all.markers=Immu_markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)

Immu_top10=Immu.all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(Immu_top10,"Immu_top10.csv")

top10 <- head(VariableFeatures(Immu_sce), 10)
VariableFeaturePlot(Immu_sce)
LabelPoints(plot = plot3, points = top10, repel = TRUE,xnudge = 0, ynudge = 0)

Immu.markers1 <- FindMarkers(Immu_sce, ident.1 = 2, min.pct = 0.25)
head(Immu.markers1)
write.csv(Immu.markers1,file = "Immu.markers1.csv")
Immu.markers1=read.csv("Immu.markers1.csv")

Immu.markers2 <- FindMarkers(Immu_sce, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
write.csv(Immu.markers2,file = "Immu.markers2.csv")
Immu.markers2=read.csv("Immu.markers2.csv")




###删掉第10簇##
Immu_sce1=Immu_sce
library(forcats)
Immu_sce1 <- subset(Immu_sce, idents = 10, invert = TRUE)
Immu_sce1$seurat_clusters <- fct_recode(Immu_sce1$seurat_clusters, 
                                        "10" = "11", 
                                        "11" = "12") %>%
  droplevels()  # 删除未使用的level（原10）
table(Immu_sce1$seurat_clusters)

col=c("#FF6F00FF", "#C71000FF" ,"#008EA0FF" ,"#8A4198FF", "#FF69B4" ,"#7B68EE", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF", "#ADE2D0FF" ,"#4B0082" ,"#D8BFD8")
DimPlot(Immu_sce1,reduction = "umap",label = T,shuffle = F,raster = F,cols = col)
ggsave(file="Immu_umap.pdf",width = 6,height = 6)
DimPlot(Immu_sce,reduction = "tsne",label = T,shuffle = F,raster = F)
ggsave(file="Immu_tsne.pdf",width = 6,height = 6)

DimPlot(Immu_sce,reduction = "umap",group.by = 'orig.ident',shuffle = F,raster = F,cols = col)
#DimPlot(Immu_sce,reduction = "tsne",group.by = 'orig.ident',raster = F,cols = Biocols)
ggsave(file="Immu_orig_umap.pdf",width = 8,height = 6)

DimPlot(Immu_sce1,reduction = "umap",split.by = "tissue_type",label = T,raster = F,cols = col)
ggsave(file="Immu_tissue-map.pdf",width = 10,height = 6)
DimPlot(Immu_sce,reduction = "tsne",split.by = "tissue_type",label = T,raster = F,cols = Biocols)
ggsave(file="Immu_tissue-tsne.pdf",width = 10,height = 6)

DimPlot(Immu_sce1,reduction = "umap",group.by = "tissue_type",shuffle = T,raster = F,cols = c("#db5a6b","#0077ce"))
ggsave(file="tissue_umap.pdf",width = 6,height = 6)
DimPlot(Immu_sce,reduction = "tsne",group.by = "tissue_type",shuffle = T,raster = F,cols = Biocols)
ggsave(file="tissue_tsne.pdf",width = 6,height = 6)



# markers <- c(
#               # "CD3D","CD3E","CD8A","CTLA4","HAVCR2","PDCD1",#Tcell
#               # "JCHAIN","CD79A","CD27","CD38",#Plasma
#               # "IL1B","OLR1","C5AR1","CD4","CD14","C1QA",#Mono/Macro
#               # "CD19","CD22",#Bcell
#               # "PEG3","NID2","FOXL2","MEG8","WFDC1","ABCA8",#germ cell
#               # "CXCL2","EREG","CCL20",
#               "CCL4","CCL3L1","IL1RN"
#               )
#FeaturePlot(Immu_sce,features = markers,ncol=3)

#markers <- c("GZMK","CD8A","CD8B")#####CD8+T
#markers=c("IGHG1","IGHG4","IGHG3")#####Plasma
#markers=c("IL7R","CD40LG","AQP3")#####CD4+T
#markers=c("TRDC","KRT81","FGFBP2")#####NK
#markers=c("GPNMB","C1QC","TREM2")#####macro
#markers=c("FOXP3","TNFRSF4","IL2RA")#####Tregs
#markers=c("MS4A1","BANK1","VPREB3")#####B
#markers=c("FCN1","CFP","VCAN")####mono
#markers=c("LILRA4","CLEC4C","PTCRA")####DC
#markers=c("AURKB","MKI67","TYMS")####HPICs
FeaturePlot(Immu_sce,features = markers,ncol=3)
Idents(Immu_sce1) <- Immu_sce1$seurat_clusters
markers=c("IGHG1","IGHG4","IGHG3",
          "GZMK","CD8A","CD8B",
          "IL7R","CD40LG","AQP3",
          "TRDC","KRT81","FGFBP2",
          "GPNMB","C1QC","TREM2",
          "AURKB","MKI67","TYMS",
          "FOXP3","TNFRSF4","IL2RA",
          "MS4A1","BANK1","VPREB3",
          "FCN1","CFP","VCAN",
          "LILRA4","CLEC4C","PTCRA")
DotPlot(Immu_sce1,features = markers) + RotatedAxis()+ 
  theme(legend.position = "none")
###### 气泡图#####
Idents(Immu_sce1)="subcelltype"
DotPlot(Immu_sce1, features = markers) +
  RotatedAxis() +  # 旋转x轴标签
  theme(
    legend.position = "none",  # 移除图例
    axis.title.x = element_blank(),  # 移除x轴标题
    axis.title.y = element_blank(),  # 移除y轴标题
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  # 调整x轴标签角度和字体
    axis.text.y = element_text(size = 12, face = "bold"),  # 调整y轴标签字体
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 设置网格线
    panel.background = element_rect(fill = "white"),  # 设置背景颜色
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 设置标题
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # 设置颜色渐变
  labs(title = "Marker Gene Expression")  # 添加标题
DotPlot(Immu_sce1, features = markers) +
  RotatedAxis() +
  theme(
    axis.title.x = element_blank(),  # 移除x轴标题
    axis.title.y = element_blank(),  # 移除y轴标题
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  # 调整x轴标签
    axis.text.y = element_text(size = 12, face = "bold"),  # 调整y轴标签
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 设置网格线
    panel.background = element_rect(fill = "white"),  # 设置背景颜色
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 设置标题
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # 设置颜色渐变
  scale_size(range = c(2, 8)) +  # 调整点的大小范围
  labs(title = "Marker Gene Expression Bubble Plot")  # 添加标题
#################
Immu_sce1$subcelltype<- recode(Immu_sce1@meta.data$seurat_clusters,
                              "0"="Plasma",
                              "1"="CD8+ T",
                              "2"="CD4+ T",
                              "3"="NK",
                              "4"="Macrophage",
                              "5"="CD8+ T",
                              "6"="HPICs",
                              "7"="Tregs",
                              "8"="NK",
                              "9"="B",
                              "10"="Monocytes",
                              "11"="DC")

table(Immu_sce1@meta.data$subcelltype)
Immu_sce1@meta.data$subcelltype <- factor(Immu_sce1@meta.data$subcelltype, 
                                          levels = setdiff(levels(Immu_sce1@meta.data$subcelltype), "12"))
cols=pal_lancet()(9)
cols=c(cols,"#FF69B4")
DimPlot(Immu_sce1,reduction = "umap",label=T,group.by = "subcelltype",cols = cols)
ggsave("Immu_anno_umap2.pdf",width=7,height=6)

DimPlot(Immu_sce,reduction = "umap",label=F,group.by = "subcelltype",split.by="tissue_type",cols = col)
ggsave("Immu_anno_umap_tissue.pdf",width=12,height=6)

DimPlot(Immu_sce,reduction = "tsne",label=T,group.by = "subcelltype",cols = Biocols)
ggsave("Immu_anno_tsne.pdf",width=7,height=6)

DimPlot(Immu_sce,reduction = "tsne",label=F,group.by = "subcelltype",split.by="tissue_type",cols = Biocols)
ggsave("Immu_anno_tsne_tissue.pdf",width=12,height=6)

Immu_sce1$celltype<- recode(Immu_sce1@meta.data$seurat_clusters,
                               "0"="B_cells",
                               "1"="T_cells",
                               "2"="T_cells",
                               "3"="ILCs",
                               "4"="Myeloid",
                               "5"="T_cells",
                               "6"="HPICs",
                               "7"="T_cells",
                               "8"="ILCs",
                               "9"="B_cells",
                               "10"="Myeloid",
                               "11"="Myeloid")
colss=pal_startrek()(4)
colss=c(colss,"#8A4198FF")
DimPlot(Immu_sce1,reduction = "umap",label=T,group.by = "celltype",cols = colss)
Immu_sce1@meta.data$celltype <- factor(Immu_sce1@meta.data$celltype, 
                                          levels = setdiff(levels(Immu_sce1@meta.data$celltype), "12"))

table(Idents(Immu_sce1))
Idents(Immu_sce1)="celltype"
Immu_markers=FindAllMarkers(Immu_sce1,test.use = "wilcox",
                            only.pos = T,
                            logfc.threshold = 0.25,
                            min.pct = 0.25)
write.csv(Immu_markers,file = "Immu_markers_anno5.csv")
Immu_markers=read.csv("Immu_markers.csv")
Immu.all.markers=Immu_markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)

Immu_top10=Immu.all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(Immu_top10,"Immu2_top10.csv")
marker=c("FCRL5","IGHG2",##B
         "CD3D","TRAC",##T
         "TRDC","GNLY",##ilc
         "LYZ","IL1B",##my
         "RRM2","AURKB"##HPICs
         )
#marker=c("RRM2","AURKB")
###########marker#####
DotPlot(Immu_sce1, features = marker) +
    RotatedAxis() +  # 旋转x轴标签
    theme(
      legend.position = "none",  # 移除图例
      axis.title.x = element_blank(),  # 移除x轴标题
      axis.title.y = element_blank(),  # 移除y轴标题
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  # 调整x轴标签角度和字体
      axis.text.y = element_text(size = 12, face = "bold"),  # 调整y轴标签字体
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 设置网格线
      panel.background = element_rect(fill = "white"),  # 设置背景颜色
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 设置标题
    ) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # 设置颜色渐变
    labs(title = "Marker Gene Expression")  # 添加标题
DotPlot(Immu_sce1, features = marker) +
  RotatedAxis() +
  theme(
    axis.title.x = element_blank(),  # 移除x轴标题
    axis.title.y = element_blank(),  # 移除y轴标题
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  # 调整x轴标签
    axis.text.y = element_text(size = 12, face = "bold"),  # 调整y轴标签
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 设置网格线
    panel.background = element_rect(fill = "white"),  # 设置背景颜色
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 设置标题
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # 设置颜色渐变
  scale_size(range = c(2, 8)) +  # 调整点的大小范围
  labs(title = "Marker Gene Expression Bubble Plot")  # 添加标题
###############
table(Idents(Immu_sce1))
Idents(Immu_sce1)="celltype"
DotPlot(Immu_sce1,features = marker) + RotatedAxis()
ggsave("Markers_Dot.pdf",width=14,height=6)
Immu.markers3<-FindAllMarkers(Immu_sce1,only.pos = TRUE , logfc.threshold = 1 ,min.pct = 0.3)
write.csv(Immu.markers3,file="Immu.markers.celltype.csv")
Immu.markers3=read.csv("Immu.markers.celltype.csv")
###备份一下
Immu_sce1=Immu_sce
top10<-Immu.markers3%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
###对这些基因进行scaledata
markers=as.data.frame(top10[,"gene"])
Immu_sce1=ScaleData(Immu_sce1,features = as.character(unique(markers$gene)))
p = DoHeatmap(Immu_sce1,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
ggsave("heatmap.pdf",p,width = 10,height = 12)

##每个细胞亚群抽取最少的细胞类型
allCells=names(Idents(Immu_sce1))
allType=levels(Idents(Immu_sce1))
choose_Cells=unlist(lapply(allType,function(x){
  cgCells=allCells[Idents(Immu_sce1)==x]
  cg=sample(cgCells,min(table(Immu_sce1@meta.data$subcelltype)))
}))

##提取
cg_sce=Immu_sce1[,allCells %in% choose_Cells]
table(Idents(cg_sce))

p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p+theme(axis.text.x = element_text(angle = 0, size = 12))
ggsave("heatmap157.pdf",p,width = 10,height = 12)

#########堆积图
cell.prop<-as.data.frame(prop.table(table(Immu_sce1@meta.data$celltype,Immu_sce1@meta.data$tissue_type)))
colnames(cell.prop)<-c("cluster","group","proportion")
# 
# col=col<-pal_npg()(10)
# 绘制条形图
p <- ggplot(cell.prop, aes(group, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colss) +  # 使用自定义的蓝色填充色
  ggtitle("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整x轴标签角度
    legend.position = "right"  # 图例位置
  ) +
  guides(fill = guide_legend(title = NULL))  # 移除图例标题
p



save(Immu_sce,file = "Immune cells2.Rdata")

######CytoTRACE#########
exp1=as.matrix(GetAssayData(Immu_sce1,assay = "RNA",layer = "counts"))
exp1=exp1[apply(exp1>0,1,sum)>=5,]
results=CytoTRACE(exp1,ncores = 1)
phenot1=Immu_sce1$celltype
phenot1=as.character(phenot1)
names(phenot1)=rownames(Immu_sce1@meta.data)
emb=Immu_sce1@reductions[["umap"]]@cell.embeddings
par(cex.main = 0.8)  # 调整主标题字体大小
plotCytoTRACE(results,phenotype = phenot1,emb = emb,outputDir = "./CytoTRACE")
plotCytoGenes(results,numOfGenes = 30,outputDir = "./CytoTRACE")
save(Immu_sce,file = "Immune cells.Rdata")
load("Immune cells.Rdata")

###################Monocle2########################
DimPlot(Immu_sce1,reduction="umap",group.by = "celltype",label = T,cols = Biocols)+
  DimPlot(Immu_sce,reduction="umap",group.by = "seurat_clusters",label = T,cols = Biocols)
ggsave("./注释/anno+clusters.pdf",width = 14,height = 6)

Idents(Immu_sce)=Immu_sce$subcelltype
table(Immu_sce1$celltype)
Immu_sce1@meta.data$seurat_clusters <- factor(Immu_sce1@meta.data$seurat_clusters, 
                                       levels = setdiff(levels(Immu_sce1@meta.data$seurat_clusters), "12"))


colnames(Immu_sce1@meta.data)

data=GetAssayData(Immu_sce1,assay = "RNA",layer = "counts")###=exprmatrix
data=as.matrix(data)
# exp1=as.matrix(GetAssayData(Immu_sce,assay = "RNA",layer = "counts"))
pd=new('AnnotatedDataFrame',data=Immu_sce1@meta.data[,c(1,6,28,35)])
fData=data.frame(gene_short_name=row.names(data),row.names=row.names(data))
fd=new('AnnotatedDataFrame',data=fData)
dim(data)  # 21425  6674

# 检查 phenoData 的维度
dim(pd) #6655           4 
if (nrow(pd) != ncol(data)) {
  pd <- pd[colnames(data), ]
}
# 检查 featureData 的维度
dim(fd)  # 21425

mycds=newCellDataSet(data,
                     phenoData =pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
mycds=estimateSizeFactors(mycds)
mycds=estimateDispersions(mycds,cores=8)

save(mycds,file = "mycds.Rdata")
load("mycds.Rdata")

disp_table=dispersionTable(mycds)
order.genes=subset(disp_table,mean_expression>=0.005&dispersion_empirical>= +
                     1*dispersion_fit)%>%pull(gene_id)%>%as.character()
mycds=setOrderingFilter(mycds,order.genes)
plot_ordering_genes(mycds)

p=plot_ordering_genes(mycds)
ggsave("./拟时序/0OrderGenes.pdf",p,width = 8,height = 6)

# 检查 mycds 的细胞数量
ncol(exprs(mycds))

# 检查 phenoData 的样本名称
rownames(pData(mycds))

# 检查 phenoData 的列名
colnames(pData(mycds))

# 如果样本名称不一致，修复 phenoData
if (!identical(colnames(exprs(mycds)), rownames(pData(mycds)))) {
  rownames(pData(mycds)) <- colnames(exprs(mycds))
}

# 如果 orig.ident 不存在，手动添加
if (!"orig.ident" %in% colnames(pData(mycds))) {
  pData(mycds)$orig.ident <- "sample1"  # 假设所有细胞属于同一个样本
}

# 重新运行 reduceDimension
mycds <- reduceDimension(mycds,
                         max_components = 2,
                         reduction_method = 'DDRTree',
                         residualModelFormulaStr = "~orig.ident")

#mycds=reduceDimension(mycds,max_components = 2,reduction_method = 'DDRTree',
#                     residualModelFormulaStr = "~orig.ident")

mycds=orderCells(mycds)
##若报错
trace('project2MST',edit=T,where=asNamespace("monocle"))
#
##mycds=orderCells(mycds,root_state=5)

plot_cell_trajectory(mycds,color_by = "State")
ggsave("./拟时序/1Trajectory_State.pdf",width = 10,height = 6.5)

plot_cell_trajectory(mycds,color_by = "Pseudotime")
ggsave("./拟时序/2Trajectory_Pseudotime.pdf",width = 10,height = 6.5)

plot_cell_trajectory(mycds,color_by = "celltype")
ggsave("./拟时序/3Trajectory_tissue.pdf",width = 10,height = 6.5)

plot_cell_trajectory(mycds,color_by = "orig.ident")
ggsave("./拟时序/4Trajectory_Orig.pdf",width = 10,height = 6.5)

ggplot(pData(mycds),aes(Pseudotime,colour=celltype,fill=celltype))+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()
ggsave("./拟时序/5Trajectory_density.pdf",width = 10,height = 6.5)

marker=c("FCRL5","IGHG2",##B
         "CD3D","TRAC",##T
         "TRDC","GNLY",##ilc
         "LYZ","IL1B",##my
         "RRM2","AURKB"##HPICs
)
#genes=c("FOXP3","PAX5","ID2","CSF1R")
genes=c("FCRL5","CD3D","TRDC","CSF1R","RRM2")
#p1=plot_genes_in_pseudotime(mycds[genes],color_by = "State")
p2=plot_genes_in_pseudotime(mycds[genes],color_by = "Pseudotime")
p3=plot_genes_in_pseudotime(mycds[genes],color_by = "celltype")
p2|p3
ggsave("./拟时序/6Trajectory_pseudotime.pdf",p1|p2|p3,width = 10,height = 6.5)

p1=plot_genes_jitter(mycds[genes],color_by = "celltype")
p2=plot_genes_violin(mycds[genes],grouping = "State",color_by = "State")
p3=plot_genes_in_pseudotime(mycds[genes],color_by = "State")
ggsave("./拟时序/7Trajectory_jitter.pdf",p1|p2|p3,width = 10,height = 6.5)

pData(mycds)$CD3D=log2(exprs(mycds)["CD3D",]+1)
p1=plot_cell_trajectory(mycds,color_by = "CD3D")+
  scale_color_continuous(type = "viridis")

pData(mycds)$CSF1R=log2(exprs(mycds)["CSF1R",]+1)
p2=plot_cell_trajectory(mycds,color_by = "CSF1R")+
  scale_color_continuous(type = "viridis")

pData(mycds)$FCRL5=log2(exprs(mycds)["FCRL5",]+1)
p3=plot_cell_trajectory(mycds,color_by = "FCRL5")+
  scale_color_continuous(type = "viridis")

pData(mycds)$RRM2=log2(exprs(mycds)["RRM2",]+1)
p4=plot_cell_trajectory(mycds,color_by = "RRM2")+
  scale_color_continuous(type = "viridis")

pData(mycds)$TRDC=log2(exprs(mycds)["TRDC",]+1)
p5=plot_cell_trajectory(mycds,color_by = "TRDC")+
  scale_color_continuous(type = "viridis")

ggsave("./拟时序/8Trajectory_experssion.pdf",p1|p2|p3|p4|p5,width = 10,height = 4)

###寻找拟时序差异基因——————monocle
Time_diff=differentialGeneTest(mycds,cores = 8,
                               fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"Time_diff.csv",row.names = F)
Time_diff=read.csv("./拟时序/DEG--monocle/Time_diff.csv")

Time_genes=Time_diff[order(Time_diff$qval),"gene_short_name"][1:100]
p=plot_pseudotime_heatmap(mycds[Time_genes,],num_clusters = 3,
                          show_rownames = T,return_heatmap = T)
ggsave("./拟时序/DEG--monocle/9Time_heatmap.pdf",p,width=5,height=10)
hp.genes=p$tree_row$labels[p$tree_row$order]
Time_diff_sig=Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"./拟时序/DEG--monocle/Time_diff_sig.csv",row.names = F)

####寻找拟时序差异基因——————beam
beam_res=BEAM(mycds,branch_point = 2,cores = 8,progenitor_method="duplicate")
trace('buildBranchCellDataSet',edit=T,where=asNamespace("monocle"))

write.csv(beam_res,"./拟时序/DEG--BEAM/Beam_all.csv",row.names = F)
beam_res=read.csv("./拟时序/DEG--BEAM/Beam_all.csv")
beam_genes=beam_res[order(beam_res$qval),"gene_short_name"][1:100]
p=plot_genes_branched_heatmap(mycds[beam_genes,],branch_point = 1,num_clusters = 3,show_rownames = T,return_heatmap = T)
ggsave("./拟时序/DEG--BEAM/10Beam_heatmap.pdf",p$ph_res,width=6.5,height=10)

hp.genes=p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
Beam_sig=beam_res[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Beam_sig,"./拟时序/DEG--BEAM/Beam_sig.csv",row.names = F)

Idents(Immu_sce)


###########差异分析#####
Idents(Immu_sce1)
Immu_subcell_markers=FindAllMarkers(Immu_sce1,logfc.threshold = 0.25,
                                    min.pct = 0.2,
                                    only.pos = FALSE,
                                    # ident.1="Fibroblast",ident.2="Endothelial",
                                    # ident.3="myoFibroblast",ident.4="Mesothelial",
                                    group.by="celltype",
)%>%mutate(gene=rownames(.))
write.csv(Immu_subcell_markers,"差异分析Immu_subcell_markers.csv",row.names = F)
jjVolcano(diffData=Immu_subcell_markers,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          #tile.col=Biocols,
          col.type="adjustP",
          topGeneN=5)

ggsave("./差异分析/subcell-volcane2.pdf",width = 22,height = 15)

# 创建一个逻辑向量来标识不包含小数点的基因名
no_decimal_genes <- !grepl("\\.", Immu_subcell_markers$gene)

# 过滤掉包含小数点的基因名
filtered_markers <- Immu_subcell_markers[no_decimal_genes, ]

Immu_subcell_markers_2=read.csv("差异分析Immu_subcell_markers.csv")
rownames(Immu_subcell_markers) <- Immu_subcell_markers$gene
deg_top20=filtered_markers%>%
  dplyr::group_by(cluster)%>%
  dplyr::top_n(n=20,wt=avg_log2FC)
col<-pal_npg()(10)
#Immu_sce1=JoinLayers(Immu_sce)
plot=averageHeatmap(Immu_sce1,
                    markerGene = deg_top20$gene,
                    group.by = "celltype",
                    gene.order = deg_top20$gene,
                    annoCol = T
)
ggsave("./差异分析/DEG-Heatmap.pdf",plot,width = 22,height = 15)

jjVolcano(diffData=filtered_markers,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          #tile.col=Biocols,
          col.type="adjustP",
          topGeneN=5)
ggsave("./差异分析/volcano.pdf",width = 15,height = 10)

rownames(Immu_subcell_markers) <- Immu_subcell_markers$gene
ids=bitr(filtered_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Immu_subcell_markers=merge(filtered_markers,ids,by.x='gene',by.y='SYMBOL')
Immu_subcell_markers=Immu_subcell_markers[order(Immu_subcell_markers$avg_log2FC,decreasing = T),]
Immu_subcell_markers.list=as.numeric(Immu_subcell_markers$avg_log2FC)
names(Immu_subcell_markers.list)=Immu_subcell_markers$ENTREZID
Immu_ana.de=names(Immu_subcell_markers.list)[abs(Immu_subcell_markers.list)>1]
##GO
Go.ana=enrichGO(Immu_ana.de,OrgDb = "org.Hs.eg.db",ont="BP",readable = T)
dotplot(Go.ana,showCategory=10,title="GO-BP")
ggsave("./差异分析/GO.BP.png",width = 8,height = 7)
Go.ana2=enrichGO(Immu_ana.de,OrgDb = "org.Hs.eg.db",ont="MF",readable = T)
dotplot(Go.ana2,showCategory=10,title="GO-MF")
ggsave("./差异分析/GO.MF.png",width = 8,height = 7)
Go.ana2=enrichGO(Immu_ana.de,OrgDb = "org.Hs.eg.db",ont="ALL",readable = T)
dotplot(Go.ana2,showCategory=10)
##KEGG
KEGG.ana=enrichKEGG(Immu_ana.de,organism = "hsa",pvalueCutoff = 0.05)
dotplot(KEGG.ana,showCategory=15,title="KEGG")
ggsave("./差异分析/KEGG.png",width = 8,height = 7)
##GSEA
GSEA.ana=gseKEGG(Immu_subcell_markers.list,organism = "hsa",pvalueCutoff = 0.05)
GSEA.ana.arrange=arrange(GSEA.ana,desc(abs(NES)))

gsekp1=gseaplot2(GSEA.ana.arrange,1:5,pvalue_table=F,base_size=14)
gsekp2=upsetplot(GSEA.ana.arrange,n=5)
cowplot::plot_grid(gsekp1,gsekp2,rel_widths = c(1,.6),labels = c("A","B"))
ggsave("./差异分析/GSEA.png",gsekp1|gsekp2,width = 14,height = 7)



##################################cnv#############
setwd("F:/毕设/data3/immune cell/IMMU/infercnv")
DimPlot(Immu_sce,reduction = "umap",label=T,split.by = "tissue_type")
DimPlot(Immu_sce,reduction = "umap",label=F,group.by = "orig.ident")

table(Immu_sce$seurat_clusters)
table(Immu_sce$tissue_type)
###抽样
###ovar2_part=ovar2[,sample(1:ncol(ovar2),round(ncol(ovar2)/5))]

Immu_matrix_counts=as.matrix(GetAssayData(Immu_sce,assay = "RNA" , layer = "counts"))
###创建新的分组
Immu_sce$Immu_cnv_type=paste0(Immu_sce$tissue_type,"_",Immu_sce$subcelltype)
unique(Immu_sce$Immu_cnv_type)
####根据cluster分析
write.table(Immu_sce$Immu_cnv_type,"./infercnv/2_Immu_cnv_type.txt",sep = "\t",quote = F , col.names = F)

gencode=read_tsv("F:/毕设/data3/hg38_gencode_v27.txt",col_names = c("gene","chr","start","end"))
gencode=gencode[!duplicated(gencode$gene),]
###提取基因交集
common_genes=intersect(gencode$gene,rownames(Immu_matrix_counts))
sort(unique(Immu_sce$Immu_cnv_type))####$后改

###创建Infercnv的object####cluster
infercnv_Immuobj=CreateInfercnvObject(raw_counts_matrix = Immu_matrix_counts[common_genes,],
                                      annotations_file = "./infercnv/2_Immu_cnv_type.txt",
                                      delim = "\t",
                                      gene_order_file = "F:/毕设/data3/hg38_gencode_v27.txt",
                                      ref_group_names = c("Normal_B cell","Normal_germ cell","Normal_Mono/Macro" ,
                                                          "Normal_Plasma cell","Normal_T cell"))
save(infercnv_Immuobj,file = "./infercnv/infercnv_Immuobj.Rdata")
load("infercnv_Immuobj.Rdata")

infercnv_obj=infercnv::run(infercnv_Immuobj,
                           cutoff = 0.1,
                           out_dir = "./infercnv/",
                           no_prelim_plot = T,
                           cluster_by_groups = T,
                           denoise = T,
                           HMM = F,
                           min_cells_per_gene = 10,
                           num_threads = 10,
                           write_expr_matrix = T)
save(infercnv_obj,file = "infercnv_obj.Rdata")
load("infercnv_Immuobj.Rdata")


grp=read.table("./infercnv/infercnv.observation_groupings.txt",sep = "",header = T,row.names = NULL)
obs=data.table::fread("./infercnv/infercnv.observations.txt",data.table = F)%>%column_to_rownames(var = 'V1')

cnvScore=function(data){
  data=data%>%
    as.matrix()%>%
    t()%>%
    scale()%>%
    rescale(to=c(-1,1))%>%
    t()
  cnv_score=as.data.frame(colSums(data*data))
  return(cnv_score)
}
cnv_score=cnvScore(obs)

grp$cluster=sapply(strsplit(grp$Dendrogram.Group,"_"),function(x) paste(x[1],x[2],sep = "_"))
sort(unique(grp$cluster))
grp$cluster=factor(grp$cluster,levels=c("Cancer_B cell","Cancer_germ cell" ,"Cancer_Mono/Macro","Cancer_Plasma cell","Cancer_T cell"))
row.names(grp)=row.names(cnv_score)
grp=grp[-1]
identical(row.names(cnv_score),row.names(grp))
cnv_score$cluster=grp$cluster
colnames(cnv_score)=c("Score","Cluster")
p=ggboxplot(cnv_score,"Cluster","Score",fill="Cluster")+
  coord_cartesian(ylim = c(1,1000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(cnv_score$Cluster)

ggplot(cnv_score,aes(x=Cluster,y=Score,fill=Cluster))+
  geom_violin()+
  coord_cartesian(ylim = c(0,500))+
  theme_minimal()+
  labs(x="Cluster",y="CNV Score")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))




#####################拟bulk差异分析################
bs=split(colnames(Immu_sce1),Immu_sce1$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(Immu_sce1)%in%bs[[x]]
    data=GetAssayData(Immu_sce1,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(Immu_sce1@meta.data[,c('orig.ident',"tissue_type")])
group_list=phe[match(names(bs),phe$orig.ident),'tissue_type']
exprset=exprset[apply(exprset,1, function(x) sum(x>1)>1),]
colData=data.frame(row.names = colnames(exprset),group_list=group_list)
dds=DESeqDataSetFromMatrix(countData=round(exprset),
                           colData=colData,
                           design=~group_list)
dds2=DESeq(dds)
tmp=results(dds2,contrast = c("group_list","Normal","Cancer"))
DEG_DESeq2=as.data.frame(tmp[order(tmp$padj),])
DEG_DESeq2=na.omit(DEG_DESeq2)

DEG_DESeq2=DEG_DESeq2%>%
  mutate(Type=if_else(pvalue>0.05,"ns",
                      if_else(abs(log2FoldChange)<0.5,"ns",
                              if_else(log2FoldChange>=0.5,"up","down"))))%>%
  arrange(desc(abs(log2FoldChange)))%>%rownames_to_column("Gene_Sambol")
write.csv(DEG_DESeq2,"DEG_DESeq2.csv")
table(DEG_DESeq2$Type)
DEG_DESeq22<- na.omit(DEG_DESeq2)

rows_to_keep <- !grepl("-", DEG_DESeq22[[1]])  # 注意这里使用df[[1]]来引用第一列

# 使用逻辑索引来选择不匹配（即不包含'-'）的行
DEG_DESeq22_cleaned <- DEG_DESeq22[rows_to_keep, ]

p=ggplot(DEG_DESeq2,aes(log2FoldChange,-log10(pvalue)))+
  geom_point(size=1.5,
             alpha=0.8,
             aes(color=Type),
             show.legend = T)+
  scale_color_manual(values = c("#00468B","gray","#E64B35"))+
  ylim(0,15)+
  xlim(-10,10)+
  labs(x="Log2(fold change)",y="-log10(padj)")+
  geom_hline(yintercept = -log10(0.05),
             linetype=2,
             color='black',lwd=0.8)+
  geom_vline(xintercept = c(-1,1),
             linetype=2,
             color='black',lwd=0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave("Myeloid_volcano_unlabel.png",p,width = 8,height = 7)

top_20 <-bind_rows(
  DEG_DESeq2 %>%
    filter(Type == 'up') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(14),
  DEG_DESeq2 %>%
    filter(Type == 'down') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(24)
)
#write.csv(top_20,"Myeloid_DEG_top10.csv")
rownames(top_20)=top_20$Gene_Sambol

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
volc_plot2
ggsave("Myeloid_volcano.png",volc_plot2,width = 8,height = 7)

##-------------------Go
sig_DEG=subset(DEG_DESeq2,Type!='ns')
genelist=bitr(sig_DEG$Gene_Sambol,fromType = "SYMBOL",
              toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
ego1=enrichGO(gene=genelist$ENTREZID,
              OrgDb =org.Hs.eg.db,
              ont = "ALL",
              minGSSize = 1,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.2,
              readable = T)
ego1_res=ego1@result
ego1_res$Group="Cancer vs Normal"
ego1_res$LogP=-log(ego1_res$p.adjust)
# ###选取感兴趣的基因集
# select_ego1_res=ego1_res%>%
#   dplyr::filter(grepl("pathway|signal|cancer",Description))%>%
#   dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
#   mutate(Description=forcats::fct_inorder(Description))

set.seed(123)

# 随机选择10个BP、MF、CC的Description
selected_bp <- sample(ego1_res$Description[ego1_res$ONTOLOGY == "BP"],5 )
selected_mf <- sample(ego1_res$Description[ego1_res$ONTOLOGY == "MF"], 5)
selected_cc <- sample(ego1_res$Description[ego1_res$ONTOLOGY == "CC"], 5)

# 合并选择的Description
selected_descriptions <- c(selected_bp, selected_mf, selected_cc)

# 筛选数据
filtered_go_results <- ego1_res[ego1_res$Description %in% selected_descriptions, ]
ggplot(filtered_go_results, aes(x = reorder(Description, RichFactor), y = RichFactor, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +  # 按ONTOLOGY分面，纵向排列
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "GO Enrichment Analysis", x = "Description", y = "Rich Factor") +
  scale_fill_manual(values = c("BP" = "blue", "MF" = "red", "CC" = "green")) +
  coord_flip()  # 将条形图改为横向排列
# 绘制条形图，按logp设置渐变颜色

ggplot(filtered_go_results, aes(x = reorder(Description,Count), y = Count, fill = LogP)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1)+  # 按ONTOLOGY分面，纵向排列
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "GO Enrichment Analysis", 
       x = "Description", 
       y = "Count", 
       fill = "LogP") +  # 图例标题
  # scale_fill_gradient(low = "grey", high = "purple") +  # 设置渐变颜色
  coord_flip()  # 将条形图改为横向排列
# 

ontology_colors <- list(
  BP = c("lightblue", "darkblue")
  CC = c("pink", "red"),
  MF = c("lightgreen", "darkgreen"),  # 分子功能（MF）
    # 细胞组分（CC）
 # 生物过程（BP）
)
p <- ggplot(filtered_go_results) +
  geom_bar(
    data = subset(filtered_go_results, ONTOLOGY == "BP"),
    aes(x = reorder(Description, Count), y = Count, fill = LogP),
    stat = "identity"
  ) +
  scale_fill_gradientn(
    colors = ontology_colors[["BP"]],  # BP颜色
    name = "LogP (BP)"
  ) +
  new_scale_fill() +  # 添加新的填充尺度
  geom_bar(
    data = subset(filtered_go_results, ONTOLOGY == "MF"),
    aes(x = reorder(Description, Count), y = Count, fill = LogP),
    stat = "identity"
  ) +
  scale_fill_gradientn(
    colors = ontology_colors[["MF"]],  # MF颜色
    name = "LogP (MF)"
  ) +
  new_scale_fill() +  # 添加新的填充尺度
  geom_bar(
    data = subset(filtered_go_results, ONTOLOGY == "CC"),
    aes(x = reorder(Description, Count), y = Count, fill = LogP),
    stat = "identity"
  ) +
  scale_fill_gradientn(
    colors = ontology_colors[["CC"]],  # CC颜色
    name = "LogP (CC)"
  ) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "GO Enrichment Analysis", x = "Description", y = "Count") +
  coord_flip()# ggplot(ego1_res[1:40,],aes(Count,Description))+
#   geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
#   #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
#   ggtitle("GO-BP_all")+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
# 
# ggplot(select_ego1_res[1:40,],aes(Count,Description))+
#   geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
#   #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
#   ggtitle("GO-BP_select")+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
p
#-----------------------------kegg
kk1=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk1_res=kk1@result
kk1_res$Group="Cancer vs Normal"
kk1_res$LogP=-log(kk1_res$p.adjust)

top20_res <- kk1_res %>%
  arrange(p.adjust) %>%  # 按p.adjust排序
  head(20)  # 取前20个

# 计算GeneRatio
top20_res <- top20_res %>%
  mutate(GeneRatio = Count / as.numeric(sub("/\\d+", "", BgRatio)))
p <- ggplot(top20_res, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    title = "Top 20 KEGG Pathways",
    color = "-Log10(p.adjust)",
    size = "Gene Count"
  ) +
  scale_color_gradient(low = "blue", high = "red")  # 颜色渐变
p

ggplot(kk1_res[40:60,],aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  ggtitle("KEGG")+
  #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
##---------------------------------------------------------gsea
geneSet=msigdbr(species="Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
###geneSet=subset(geneSet,gs_name%in%)

geneList=DEG_DESeq2$log2FoldChange
names(geneList)=DEG_DESeq2$Gene_Sambol
geneList=sort(geneList,decreasing = T)

GSEA_enrichment=GSEA(geneList,
                     TERM2GENE = geneSet,
                     pvalueCutoff = 1,
                     minGSSize = 20,
                     maxGSSize = 1000,
                     eps=0,
                     pAdjustMethod = "BH")
result=data.frame(GSEA_enrichment)
gseaplot2(GSEA_enrichment,c(1:5),pvalue_table = T)
dotplot(GSEA_enrichment,showCategory=10,split=".sign")+
  facet_grid(~.sign)+
  theme(plot.title = element_text(size = 10,color = "black",hjust=0.5),
        axis.title = element_text(size = 10,color = "black"),
        axis.text = element_text(size = 5,color = "red"),
        axis.text.x = element_text(angle = 0,hjust = 1),
        legend.position = "right",
        legend.text = element_text(size=10),
        legend.title=element_text(size=10))+
  ggtitle("GSEA")
getwd()
setwd("F:/毕设/data3/immune cell/IMMU")
save(Immu_sce,file = "Immune cells2.Rdata")

#######################################all#######################
######################cellchat######################all
data1=GetAssayData(Immu_sce,assay = "RNA",layer = "data")
table(Immu_sce@meta.data$subcelltype)
table(Immu_sce@meta.data$orig.ident)
cellchat=createCellChat(data1,meta = Immu_sce@meta.data,
                        group.by = "subcelltype")
CellChatDB=CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB=CellChatDB
cellchat=subsetData(cellchat)
future::plan("multisession",workers=10)
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat=computeCommunProb(cellchat,raw.use = F)
cellchat=filterCommunication(cellchat,min.cells=10)
cellchat=computeCommunProbPathway(cellchat)
df.net=subsetCommunication(cellchat)
write.csv(df.net,"Gene.csv",row.names = F)
df.netP=subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netP,"Pathway.csv",row.names = F)
cellchat=aggregateNet(cellchat)
save(cellchat,file="Immu_cellchat.Rdata")
groupSize=as.numeric(table(cellchat@idents))

pdf("NetVisual_overview_all.pdf",width = 8,height = 6)
par(xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Interaction weights/strength")
dev.off()

pdf("NetVisual_overview_split2.pdf",width = 6,height = 5)
mat=cellchat2@net$weight
for(i in 1:nrow(mat)){
  mat2=matrix(0,nrow = nrow(mat),ncol=ncol(mat),dimnames = dimnames(mat))
  mat2[i,]=mat[i,]
  par(xpd=TRUE)
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,
                   edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()
mypathways=cellchat@netP$pathways
mypathways1=mypathways[1:5]
mypathways1

pdf("NetVisual_pathways_circle.pdf",width=6,height=5)
for(pathways.show in mypathways1){
  par(xpd=T)
  netVisual_aggregate(cellchat,signaling = pathways.show,layout="circle")
}
dev.off()

pdf("NetVisual_pathways_chord.pdf",width = 10,height = 8)
for(pathways.show in mypathways1){
  netVisual_aggregate(cellchat,signaling = pathways.show,layout="chord")
}
dev.off()

dir.create("Pathways")
for(pathways.show in mypathways1){
  pdf(paste0("Pathways/",pathways.show,".pdf"),width = 8,height=6.5)
  netAnalysis_contribution(cellchat,signaling = pathways.show)
  pairLR=extractEnrichedLR(cellchat,signaling = pathways.show,geneLR.return = FALSE)$interaction_name
  for(LR.show in pairLR){
    netVisual_individual(cellchat,signaling = pathways.show,pairLR.use = LR.show,layout = "circle")
  }
  for(LR.show in pairLR){
    netVisual_individual(cellchat,signaling = pathways.show,pairLR.use = LR.show,layout = "chord")
  }
  dev.off()
  }
p=netVisual_bubble(cellchat,sources.use =1:length(levels(cellchat@idents)),
                   targets.use = 1:length(levels(cellchat@idents)),remove.isolate = F)
ggsave("CCI_all.pdf",p,width = 15,height=40,limitsize = F)

p=netVisual_bubble(cellchat,sources.use = 1:5,targets.use = 6:10,remove.isolate = F)
ggsave("CCI_part.pdf",p,width = 15,height=40,limitsize = F)

p=netVisual_bubble(cellchat,sources.use ="T",targets.use = c("CD8+ T","CD4+ T","Tregs"),remove.isolate = F)
ggsave("CCI_T.pdf",p,width = 8,height=20,limitsize = F)

p=netVisual_bubble(cellchat,sources.use =c("T","B","NK"),targets.use = c("CD8+ T","CD4+ T","Tregs","Plasma"),signaling=c("CCL","TNF"),remove.isolate = F)
ggsave("CCI_T_pathway.pdf",p,width = 8,height=6,limitsize = F)

pairLR.use=c("CCL3_CCR1","CCL4_CCR5","CCL5_CCR3","TNF_TNFRSF1A","TNF_TNFRSF1b")
pairLR.use=data.frame(interaction_name=pairLR.use)
p=netVisual_bubble(cellchat,sources.use =1:length(levels(cellchat@idents)),
                   targets.use = 1:length(levels(cellchat@idents)),
                   pairLR.use = pairLR.use,remove.isolate = F)
ggsave("CCI_LR.pdf",p,width = 16,height=6,limitsize = F)


p=plotGeneExpression(cellchat,signaling = "CCL")
ggsave("GeneExpression_violin_sig.pdf",p,width = 10,height=9,limitsize = F)
p=plotGeneExpression(cellchat,signaling = "CCL",enriched.only = FALSE)
ggsave("GeneExpression_violin_sig2.pdf",p,width = 10,height=9,limitsize = F)

cellchat2=netAnalysis_computeCentrality(cellchat,slot.name = "netP")
pdf("signalingRole.pdf",width = 6,height = 4.5)
for(pathways.show in mypathways1){
  netAnalysis_signalingRole_network(cellchat2,signaling = pathways.show,
                                    width=8,height = 2.5,font.size = 10)
}
dev.off()
netAnalysis_signalingRole_scatter(cellchat2)
ggsave("function.pdf",width = 8,height=8)
netAnalysis_signalingRole_heatmap(cellchat2,pattern = "outgoing",signaling = c("CCL","CXCL","LCK","TRAIL","MHC-I","GALECTIN"))
netAnalysis_signalingRole_heatmap(cellchat2,pattern = "incoming",signaling = c("CCL","CXCL","LCK","TRAIL","MHC-I","GALECTIN"))


###########NMF聚类#######信号流出通讯模式的鉴定和可视化
selectK(cellchat2,pattern = "outgoing")
nPatterna=3
pdf("NMF2.pdf",width = 10,height = 20)
cellchat2=identifyCommunicationPatterns(cellchat2,pattern = "outgoing",
                                        k=nPatterna,width=6,height=20,font.size = 6)
dev.off()
netAnalysis_river(cellchat2,pattern = "outgoing")

netAnalysis_dot(cellchat2,pattern = "outgoing")


selectK(cellchat2,pattern = "incoming")
nPatterna=3
pdf("NMF2.pdf",width = 10,height = 20)
cellchat2=identifyCommunicationPatterns(cellchat2,pattern = "incoming",
                                        k=nPatterna,width=6,height=20,font.size = 6)
dev.off()
netAnalysis_river(cellchat2,pattern = "incoming")

netAnalysis_dot(cellchat2,pattern = "incoming")
save(cellchat,file="cellchat1.Rdata")
load("cellchat1.Rdata")
save(cellchat2,file="cellchat2.Rdata")
##############################cancer
##################################Immu_cancer####################
cell.ids1 <- which(Immu_sce@meta.data$tissue_type == "Cancer")
Immu_Cancer<- subset(Immu_sce, cells =cell.ids1)

data1=GetAssayData(Immu_Cancer,assay = "RNA",layer = "data")
table(Immu_Cancer@meta.data$subcelltype)
table(Immu_Cancer@meta.data$orig.ident)
cellchat=createCellChat(data1,meta = Immu_Cancer@meta.data,
                        group.by = "subcelltype")
CellChatDB=CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB=CellChatDB
cellchat=subsetData(cellchat)
future::plan("multisession",workers=10)
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat=computeCommunProb(cellchat,raw.use = F)

cellchat_cancer=filterCommunication(cellchat,min.cells=10)
cellchat_cancer=computeCommunProbPathway(cellchat_cancer)
df.net.cancer=subsetCommunication(cellchat_cancer)
write.csv(df.net.cancer,"Gene.csv",row.names = F)
df.netP.cancer=subsetCommunication(cellchat_cancer,slot.name = "netP")
write.csv(df.netP.cancer,"Pathway.csv",row.names = F)
cellchat_cancer=aggregateNet(cellchat_cancer)
save(cellchat_cancer,file="Immu_cancer_cellchat.Rdata")
groupSize=as.numeric(table(cellchat_cancer@idents))

pdf("Cancer_NetVisual_overview_all.pdf",width = 8,height = 6)
par(xpd=TRUE)
netVisual_circle(cellchat_cancer@net$count,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Number of interactions")
netVisual_circle(cellchat_cancer@net$weight,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Interaction weights/strength")
dev.off()

pdf("Cancer_NetVisual_overview_split.pdf",width = 6,height = 5)
mat=cellchat_cancer@net$weight
for(i in 1:nrow(mat)){
  mat2=matrix(0,nrow = nrow(mat),ncol=ncol(mat),dimnames = dimnames(mat))
  mat2[i,]=mat[i,]
  par(xpd=TRUE)
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,
                   edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()
mypathways=cellchat_cancer@netP$pathways
mypathways1=mypathways[6:10]
mypathways1

pdf("Cancer_NetVisual_pathways_circle.pdf",width=6,height=5)
for(pathways.show in mypathways1){
  par(xpd=T)
  netVisual_aggregate(cellchat_cancer,signaling = pathways.show,layout="circle")
}
dev.off()

pdf("Cancer_NetVisual_pathways_chord.pdf",width = 10,height = 8)
for(pathways.show in mypathways1){
  netVisual_aggregate(cellchat_cancer,signaling = pathways.show,layout="chord")
}
dev.off()

dir.create("Pathways")
for(pathways.show in mypathways1){
  pdf(paste0("Pathways/",pathways.show,".pdf"),width = 8,height=6.5)
  netAnalysis_contribution(cellchat_cancer,signaling = pathways.show)
  pairLR=extractEnrichedLR(cellchat_cancer,signaling = pathways.show,geneLR.return = FALSE)$interaction_name
  for(LR.show in pairLR){
    netVisual_individual(cellchat_cancer,signaling = pathways.show,pairLR.use = LR.show,layout = "circle")
  }
  for(LR.show in pairLR){
    netVisual_individual(cellchat_cancer,signaling = pathways.show,pairLR.use = LR.show,layout = "chord")
  }
  dev.off()
}
p=netVisual_bubble(cellchat_cancer,sources.use =1:length(levels(cellchat_cancer@idents)),
                   targets.use = 1:length(levels(cellchat_cancer@idents)),remove.isolate = F)
ggsave("CCI_all.pdf",p,width = 15,height=40,limitsize = F)

p=netVisual_bubble(cellchat_cancer,sources.use = 1:5,targets.use = 6:10,remove.isolate = F)
ggsave("CCI_part.pdf",p,width = 15,height=40,limitsize = F)

p=netVisual_bubble(cellchat_cancer,sources.use ="T",targets.use = c("CD8+ T","CD4+ T","Tregs"),remove.isolate = F)
ggsave("CCI_T.pdf",p,width = 8,height=20,limitsize = F)

p=netVisual_bubble(cellchat_cancer,sources.use =c("T","B","NK"),targets.use = c("CD8+ T","CD4+ T","Tregs","Plasma"),signaling=c("CCL","TNF"),remove.isolate = F)
ggsave("CCI_T_pathway.pdf",p,width = 8,height=6,limitsize = F)

pairLR.use=c("CCL3_CCR1","CCL4_CCR5","CCL5_CCR3","TNF_TNFRSF1A","TNF_TNFRSF1b")
pairLR.use=data.frame(interaction_name=pairLR.use)
p=netVisual_bubble(cellchat_cancer,sources.use =1:length(levels(cellchat_cancer@idents)),
                   targets.use = 1:length(levels(cellchat_cancer@idents)),
                   pairLR.use = pairLR.use,remove.isolate = F)
ggsave("CCI_LR.pdf",p,width = 16,height=6,limitsize = F)


p=plotGeneExpression(cellchat_cancer,signaling = "CCL")
ggsave("GeneExpression_violin_sig.pdf",p,width = 10,height=9,limitsize = F)
p=plotGeneExpression(cellchat_cancer,signaling = "CCL",enriched.only = FALSE)
ggsave("GeneExpression_violin_sig2.pdf",p,width = 10,height=9,limitsize = F)

cellchat_cancer2=netAnalysis_computeCentrality(cellchat_cancer,slot.name = "netP")
pdf("signalingRole.pdf",width = 6,height = 4.5)
for(pathways.show in mypathways1){
  netAnalysis_signalingRole_network(cellchat_cancer2,signaling = pathways.show,
                                    width=8,height = 2.5,font.size = 10)
}
dev.off()
netAnalysis_signalingRole_scatter(cellchat_cancer2)
ggsave("function.pdf",width = 8,height=8)
netAnalysis_signalingRole_heatmap(cellchat_cancer2,pattern = "outgoing",signaling = c("CCL","CXCL","LCK","TRAIL","MHC-I","GALECTIN"))
netAnalysis_signalingRole_heatmap(cellchat_cancer2,pattern = "incoming",signaling = c("CCL","CXCL","LCK","TRAIL","MHC-I","GALECTIN"))


###########NMF聚类#######信号流出通讯模式的鉴定和可视化
selectK(cellchat_cancer2,pattern = "outgoing")
nPatterna=4
pdf("NMF2.pdf",width = 10,height = 20)
cellchat_cancer2=identifyCommunicationPatterns(cellchat_cancer2,pattern = "outgoing",
                                        k=nPatterna,width=6,height=20,font.size = 6)
dev.off()
netAnalysis_river(cellchat_cancer2,pattern = "outgoing")

netAnalysis_dot(cellchat_cancer2,pattern = "outgoing")


selectK(cellchat_cancer2,pattern = "incoming")
nPatterna=3
pdf("incoming_NMF2.pdf",width = 10,height = 20)
cellchat_cancer2=identifyCommunicationPatterns(cellchat_cancer2,pattern = "incoming",
                                        k=nPatterna,width=6,height=20,font.size = 6)
dev.off()
netAnalysis_river(cellchat_cancer2,pattern = "incoming")

netAnalysis_dot(cellchat_cancer2,pattern = "incoming")

save(cellchat_cancer,file="cellchat_cancer1.Rdata")
load("cellchat1.Rdata")
save(cellchat_cancer2,file="cellchat_cancer2.Rdata")


##################################Immu_normal#####
cell.ids1 <- which(Immu_sce@meta.data$tissue_type == "Normal")
Immu_Normal<- subset(Immu_sce, cells =cell.ids1)

data1=GetAssayData(Immu_Normal,assay = "RNA",layer = "data")
table(Immu_Normal@meta.data$subcelltype)
table(Immu_Normal@meta.data$orig.ident)
cellchat=createCellChat(data1,meta = Immu_Normal@meta.data,
                        group.by = "subcelltype")
CellChatDB=CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB=CellChatDB
cellchat=subsetData(cellchat)
future::plan("multisession",workers=10)
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat=computeCommunProb(cellchat,raw.use = F)

cellchat_normal=filterCommunication(cellchat,min.cells=10)
cellchat_normal=computeCommunProbPathway(cellchat_normal)
df.net.cancer=subsetCommunication(cellchat_normal)
write.csv(df.net.cancer,"Gene.csv",row.names = F)
df.netP.cancer=subsetCommunication(cellchat_normal,slot.name = "netP")
write.csv(df.netP.cancer,"Pathway.csv",row.names = F)
cellchat_normal=aggregateNet(cellchat_normal)
save(cellchat_normal,file="Immu_Normal_cellchat.Rdata")
groupSize=as.numeric(table(cellchat_normal@idents))

pdf("Cancer_NetVisual_overview_all.pdf",width = 8,height = 6)
par(xpd=TRUE)
netVisual_circle(cellchat_normal@net$count,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Number of interactions")
netVisual_circle(cellchat_normal@net$weight,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Interaction weights/strength")
dev.off()

pdf("Cancer_NetVisual_overview_split.pdf",width = 6,height = 5)
mat=cellchat_normal@net$weight
for(i in 1:nrow(mat)){
  mat2=matrix(0,nrow = nrow(mat),ncol=ncol(mat),dimnames = dimnames(mat))
  mat2[i,]=mat[i,]
  par(xpd=TRUE)
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,
                   edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()
mypathways=cellchat_normal@netP$pathways
mypathways1=mypathways[6:10]
mypathways1

pdf("Cancer_NetVisual_pathways_circle.pdf",width=6,height=5)
for(pathways.show in mypathways1){
  par(xpd=T)
  netVisual_aggregate(cellchat_normal,signaling = pathways.show,layout="circle")
}
dev.off()

pdf("Cancer_NetVisual_pathways_chord.pdf",width = 10,height = 8)
for(pathways.show in mypathways1){
  netVisual_aggregate(cellchat_normal,signaling = pathways.show,layout="chord")
}
dev.off()

dir.create("Pathways")
for(pathways.show in mypathways1){
  pdf(paste0("Pathways/",pathways.show,".pdf"),width = 8,height=6.5)
  netAnalysis_contribution(cellchat_normal,signaling = pathways.show)
  pairLR=extractEnrichedLR(cellchat_normal,signaling = pathways.show,geneLR.return = FALSE)$interaction_name
  for(LR.show in pairLR){
    netVisual_individual(cellchat_normal,signaling = pathways.show,pairLR.use = LR.show,layout = "circle")
  }
  for(LR.show in pairLR){
    netVisual_individual(cellchat_normal,signaling = pathways.show,pairLR.use = LR.show,layout = "chord")
  }
  dev.off()
}
p=netVisual_bubble(cellchat_normal,sources.use =1:length(levels(cellchat_normal@idents)),
                   targets.use = 1:length(levels(cellchat_normal@idents)),remove.isolate = F)
ggsave("CCI_all.pdf",p,width = 15,height=40,limitsize = F)

p=netVisual_bubble(cellchat_normal,sources.use = 1:5,targets.use = 6:10,remove.isolate = F)
ggsave("CCI_part.pdf",p,width = 15,height=40,limitsize = F)

p=netVisual_bubble(cellchat_normal,sources.use =c("NK","CD8+ T","CD4+ T","B","Plasma"),targets.use = c("CD8+ T","CD4+ T","B","Plasma","NK"),remove.isolate = F)
ggsave("CCI_LB.pdf",p,width = 8,height=20,limitsize = F)

p=netVisual_bubble(cellchat_normal,sources.use =c("NK","CD8+ T","CD4+ T","B","Plasma"),targets.use = c("CD8+ T","CD4+ T","B","Plasma","NK"),
                   signaling=c("CCL","CXCL","IL4"),remove.isolate = F)
ggsave("CCI_LB_pathway.pdf",p,width = 8,height=6,limitsize = F)

pairLR.use=c("CCL3_CCR1","CCL4_CCR5","CCL5_CCR3","TNF_TNFRSF1A","TNF_TNFRSF1b")
pairLR.use=data.frame(interaction_name=pairLR.use)
p=netVisual_bubble(cellchat_normal,sources.use =1:length(levels(cellchat_normal@idents)),
                   targets.use = 1:length(levels(cellchat_normal@idents)),
                   pairLR.use = pairLR.use,remove.isolate = F)
ggsave("CCI_LR.pdf",p,width = 16,height=6,limitsize = F)


p=plotGeneExpression(cellchat_normal,signaling = "CCL")
ggsave("GeneExpression_violin_sig.pdf",p,width = 10,height=9,limitsize = F)
p=plotGeneExpression(cellchat_normal,signaling = "CCL",enriched.only = FALSE)
ggsave("GeneExpression_violin_sig2.pdf",p,width = 10,height=9,limitsize = F)

cellchat_normal2=netAnalysis_computeCentrality(cellchat_normal,slot.name = "netP")
pdf("signalingRole.pdf",width = 6,height = 4.5)
for(pathways.show in mypathways1){
  netAnalysis_signalingRole_network(cellchat_normal2,signaling = pathways.show,
                                    width=8,height = 2.5,font.size = 10)
}
dev.off()
netAnalysis_signalingRole_scatter(cellchat_normal2)
ggsave("function.pdf",width = 8,height=8)
netAnalysis_signalingRole_heatmap(cellchat_normal2,pattern = "outgoing",signaling = c("CCL","CXCL","LCK","TRAIL","MHC-I","GALECTIN"))
netAnalysis_signalingRole_heatmap(cellchat_normal2,pattern = "incoming",signaling = c("CCL","CXCL","LCK","TRAIL","MHC-I","GALECTIN"))


###########NMF聚类#######信号流出通讯模式的鉴定和可视化
selectK(cellchat_normal2,pattern = "outgoing")
nPatterna=5
pdf("outgoing_NMF2.pdf",width = 10,height = 20)
cellchat_normal2=identifyCommunicationPatterns(cellchat_normal2,pattern = "outgoing",
                                               k=nPatterna,width=6,height=20,font.size = 6)
dev.off()
netAnalysis_river(cellchat_normal2,pattern = "outgoing")

netAnalysis_dot(cellchat_normal2,pattern = "outgoing")


selectK(cellchat_normal2,pattern = "incoming")
nPatterna=4
pdf("incoming_NMF2.pdf",width = 10,height = 20)
cellchat_normal2=identifyCommunicationPatterns(cellchat_normal2,pattern = "incoming",
                                               k=nPatterna,width=6,height=20,font.size = 6)
dev.off()
netAnalysis_river(cellchat_normal2,pattern = "incoming")

netAnalysis_dot(cellchat_normal2,pattern = "incoming")

save(cellchat_normal,file="cellchat_normal1.Rdata")
save(cellchat_normal2,file="cellchat_normal2.Rdata")

##################################两组合并#############
Immu_all_cellchat=list(normal=cellchat_normal2,cancer=cellchat_cancer2)
Immu_all_cellchat=mergeCellChat(Immu_all_cellchat,add.names = names(Immu_all_cellchat))
###互作次数bar图
p1=compareInteractions(Immu_all_cellchat,show.legend = F,group = c(1,2))
p2=compareInteractions(Immu_all_cellchat,show.legend = F,group = c(1,2),measure = "weight")
p1+p2
ggsave("num_interactions.png",p1+p2,width = 8,height=8)

#####互作次数网络图
par(mfrow=c(1,2),xpd=T)
netVisual_diffInteraction(Immu_all_cellchat,weight.scale = T)
netVisual_diffInteraction(Immu_all_cellchat,weight.scale = T,measure = "weight")
dev.off()

###数量与强度差异热图
par(mfrow=c(1,1))
p1=netVisual_heatmap(Immu_all_cellchat)
p2=netVisual_heatmap(Immu_all_cellchat,measure = "weight")
p1+p2

###保守和特异信号通路的识别和可视化
p1=rankNet(Immu_all_cellchat,mode = "comparison",stacked = T,do.stat = T)
p2=rankNet(Immu_all_cellchat,mode = "comparison",stacked = F,do.stat = T)
p1+p2


diff.count=Immu_all_cellchat@net$cancer$count-Immu_all_cellchat@net$normal$count
write.csv(Immu_all_cellchat@net$cancer$count,"output.Cancer.count.csv",quote=F)
write.csv(Immu_all_cellchat@net$normal$count,"output.Normal.count.csv",quote=F)

pheatmap(diff.count,
         treeheight_row = "0",
         treeheight_col = "0",
         cluster_rows = T,
         cluster_cols = T)

view(Immu_all_cellchat)
view(Immu_all_cellchat@LR[["cancer"]])
