library("dplyr")
library("patchwork")
library("Seurat")
library("ggplot2")
library("ggrepel")
library("reshape2")
library("data.table")
library("stringr")
library("tibble")
library("harmony")
library("clustree")
library("infercnv")
library("readr")
library("RColorBrewer")
library("scales")
library("ggpubr")
library("devtools")
library("gplots")
library("sva")
#devtools::install_local("D:/下载/CytoTRACE_0.3.3.tar.gz")
library("CytoTRACE")
#devtools::install_github('junjunlab/scRNAtoolVis')
library("scRNAtoolVis")
library("clusterProfiler")
#BiocManager::install('AnnotationDbi',dependencies=TRUE)
library("AnnotationDbi")
#install.packages("D:/下载/org.Hs.eg.db_3.19.1.tar.gz",repos = NULL,type = "source")
library("org.Hs.eg.db")
library("enrichplot")
library("ggupset")
library("monocle")
library("igraph")
library("ggsci")
library("scales")
library("msigdbr")
library("DESeq2")
library("NMF")
library("circlize")
library("CellChat")
library("assorthead")
library("BiocNeighbors")
library("GEOquery")
library(reactome.db)
library(ReactomePA)
library(ggnewscale)
library(survival) 
library(survminer)
# 查看当前工作目录
getwd()
# 设置工作目录（将工作目录切换到指定路径下）
setwd("F:/毕设/data3")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files("./GSE184880")

# 创建一个空的列表来存储Seurat对象
seurat_list <- list()
# 读取每个样本的10x数据并创建Seurat对象
for (sample in samples) {
  # 拼接文件路径
  data.path <- paste0("./GSE184880/", sample)
  
  # 读取10x数据，data.dir参数指定存放文件的路径
  seurat_data <- Read10X(data.dir = data.path)
  
  # 创建Seurat对象，并指定项目名称为样本文件名
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample,
                                   min.features = 1000,
                                   min.cells = 5)
  
  # 将Seurat对象添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

#合并对象
ovar <- merge(seurat_list[[1]], 
              y = seurat_list[-1],
              add.cell.ids = samples)
ovar<-JoinLayers(ovar)
#table(ovar@meta.data$orig.ident)
#table(str_split(colnames(ovar),'-',simplify = T)[,2])

# 使用PercentageFeatureSet()函数计算线粒体 QC 指标，该函数计算来自一组功能的计数百分比
ovar[["percent.mt"]] <- PercentageFeatureSet(ovar, pattern = "^MT-")

#红细胞比例
HB.genes<-c("HBA1","HBA2","HBB","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ","HBD")
HB.genes<-CaseMatch(HB.genes,rownames(ovar))
ovar[["percent.HB"]] <- PercentageFeatureSet(ovar, features = HB.genes)

# 查看相关性
FeatureScatter(ovar,"nCount_RNA","percent.mt",group.by = "orig.ident")
FeatureScatter(ovar,"nCount_RNA","nFeature_RNA",group.by = "orig.ident")
#查看指控指标
#设置绘图元素
theme.set2=theme(axis.title = element_blank())
plot.features=c("nFeature_RNA","nCount_RNA","percent.mt")
group="orig.ident"
#指控前小提琴图
plots=list()
for(i in c(1:length(plot.features))){
  plots[[i]]=VlnPlot(ovar,group.by=group,pt.size=0,
                     features=plot.features[i]) + theme.set2 + NoLegend()
}
violin<-wrap_plots(plots=plots,ncol=3)
ggsave("22vln_before_qc.pdf",plot = violin,width = 14,height = 8)

#设置质控指标
quantile(ovar$nFeature_RNA,seq(0.01,0.1,0.01))
quantile(ovar$nFeature_RNA,seq(0.9,1,0.01))

quantile(ovar$nCount_RNA,seq(0.01,0.1,0.01))
quantile(ovar$nCount_RNA,seq(0.9,1,0.01))

quantile(ovar$percent.mt,seq(0.01,0.1,0.01))
quantile(ovar$percent.mt,seq(0.9,1,0.01))

#质控
ovar <- subset(ovar, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & nCount_RNA > 1500 & percent.mt < 10 & percent.mt > 0)
print(ovar)

plots=list()
for(i in seq_along(plot.features)){
  plots[[i]]=VlnPlot(ovar,group.by=group,pt.size=0,
                     features=plot.features[i]) + theme.set2 + NoLegend()
}
violin<-wrap_plots(plots=plots,ncol=3)
ggsave("vln_after_qc.pdf",width = 14,height = 8)

saveRDS(ovar,"1ovar_qc.rds")

ovar=readRDS("1ovar_qc.rds")

#标准化
ovar <- NormalizeData(ovar)
# 高变异基因的选择
ovar <- FindVariableFeatures(ovar, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ovar)
ovar <- ScaleData(ovar, features = all.genes)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ovar), 10)
top10
plot3 <- VariableFeaturePlot(ovar)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE,xnudge = 0, ynudge = 0)

#ovar<-SCTransform(ovar)
plot4
ggsave("top10.pdf",width = 14,height = 13)

##PCA分析
ovar <- RunPCA(ovar, verbose = F)
#ovar1 <- RunPCA(ovar, features = VariableFeatures(object = ovar))
# Examine and visualize PCA results a few different ways
print(ovar[["pca"]], dims = 1:5, nfeatures = 5)
##主成分分析图
pdf(file = "pca.pdf",width = 7,height = 6)
DimPlot(ovar, reduction = 'pca', pt.size = 0.5)
dev.new()
##绘制每个PCA成分的相关基因
pdf(file = "PAC_rel.pdf",width = 10,height = 9)
VizDimLoadings(ovar, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()
#主成分分析热图
DimHeatmap(ovar, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30 , ncol = 2)
pdf("PCA_Heat.pdf",width = 10,height = 9)
dev.off()
#
##DimHeatmap(ovar, dims = 1:12, cells = 500, balanced = TRUE , nfeatures = 30 , ncol = 3)

##选取合适的PC
##主成分积累贡献大于90%选择拐点
pdf(file = "PC.pdf",width = 7,height = 6)
ElbowPlot(ovar,ndims = 50)
dev.off()

#确定与每个PC的累计百分比
pct<-ovar[["pca"]]@stdev/sum(ovar[["pca"]]@stdev)*100
#计算与每个PC的累计百分比
cumu<-cumsum(pct)
#设置pc
pcs=1:40
ovar<-RunHarmony(ovar,group.by.vars="orig.ident",assay.use="RNA",max.iter.harmony=20)
table(ovar@meta.data$orig.ident)
#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
seq=seq(0.1,2,by=0.1)
ovar=FindNeighbors(ovar,dims=pcs)
for(res in seq){
  ovar=FindClusters(ovar,resolution = res)
}

#画图
p1=clustree(ovar,prefix="RNA_snn_res.")+coord_flip()
p=p1+plot_layout(widths = c(3,1))
ggsave("RNA_snn_res.png",p,width = 30,height = 14)

#降维聚类h
ovar<-FindNeighbors(ovar,reduction = "harmony",dims=pcs)%>%FindClusters(resolution = 1)
ovar<-RunUMAP(ovar,reduction = "harmony",dims=pcs)%>%RunTSNE(dims=pcs,reduction = "harmony")

DimPlot(ovar,reduction = "umap",label = T)
ggsave("umap.png",width = 8,height = 6)

DimPlot(ovar,reduction = "umap",label = F,group.by="orig.ident")
ggsave("umap_oi.png",width = 8,height = 6)

DimPlot(ovar,reduction = "tsne",label = T)
ggsave("tsne.png",width = 8,height = 6)

DimPlot(ovar,reduction = "tsne",label = F,group.by="orig.ident")
ggsave("tsne_oi.png",width = 8,height = 6)
###wuh
#ovar1<-FindNeighbors(ovar,reduction = "pca",dims=pcs)%>%FindClusters(resolution = 1)
#ovar1<-RunUMAP(ovar,reduction = "pca",dims=pcs)%>%RunTSNE(dims=pcs,reduction = "pca")

#DimPlot(ovar1,reduction = "umap",label = T)
#ggsave("umap.png",width = 8,height = 6)
 
#DimPlot(ovar1,reduction = "umap",label = F,group.by="orig.ident")
#ggsave("umap_oi.png",width = 8,height = 6)
 
#DimPlot(ovar1,reduction = "tsne",label = T)
#ggsave("tsne.png",width = 8,height = 6)
 
#DimPlot(ovar1,reduction = "tsne",label = F,group.by="orig.ident")
#ggsave("tsne_oi.png",width = 8,height = 6)

saveRDS(ovar,"5ovar_Umap_Tsne.rds")

ovar=readRDS("5ovar_Umap_Tsne.rds")



#寻找样本marker
# ovar.markers1 <- FindMarkers(ovar, ident.1 = 2, min.pct = 0.25)
# write.csv(ovar.markers1,file = "markers1.RNA.csv")
ovar.markers1=read.csv("./细胞注释/markers1.RNA.csv")
# ovar.markers2 <- FindMarkers(ovar, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# write.csv(ovar.markers2,file = "markers2.RNA.csv")
ovar.markers2=read.csv("./细胞注释/markers2.RNA.csv")

##第一次注释
markers <- c("PTPRC",#immune
             "EPCAM",#epithelial
             "PECAM1","MME")#stromal

p=FeaturePlot(ovar,features = markers,ncol=3)
ggsave("1_map.pdf",p,width=15,height=10)

p=DotPlot(ovar,features = markers) + RotatedAxis()
ggsave("2_Dot.pdf",p,width=14,height=6)

p=VlnPlot(ovar,features = markers,stack = T , flip = T)
ggsave("3_Vln.pdf",p,width=14,height=6)

p=DimPlot(ovar,reduction = "umap",label=T)
DimPlot(ovar,reduction = "tsne",label=T)
ggsave("./细胞注释/第一次注释/4_tsne.pdf",width=7,height=6)

# ovar$celltype.1<- recode(ovar@meta.data$seurat_clusters,
#                          "0"="0",
#                          "1"="immune cell",
#                          "2"="2",
#                          "3"="immune cell",
#                          "4"="immune cell",
#                          "5"="Epithelial",
#                          "6"="immune cell",
#                          "7"="immune cell",
#                          "8"="8",
#                          "9"="Endothelial",
#                          "10"="immune cell",
#                          "11"="immune cell",
#                          "12"="12",
#                          "13"="13",
#                          "14"="Fibroblast",
#                          "15"="immune cell",
#                          "16"="Epithelial",
#                          "17"="17",
#                          "18"="Fibroblast",
#                          "19"="immune cell",
#                          "20"="immune cell",
#                          "21"="21",
#                          "22"="22",
#                          "23"="immune cell",
#                          "24"="immune cell",
#                          "25"="immune cell",
#                          "26"="Epithelial",
#                          "27"="Epithelial",
#                          "28"="28",
#                          "29"="immune cell",
#                          "30"="Endothelial",
#                          "31"="immune cell",
#                          "32"="0",
#                          "33"="immune cell" )

ovar2$celltype.1<- recode(ovar2@meta.data$seurat_clusters,
                         "0"="Stromal cell",
                         "1"="immune cell",
                         "2"="Stromal cell",
                         "3"="immune cell",
                         "4"="immune cell",
                         "5"="Epithelial",
                         "6"="immune cell",
                         "7"="immune cell",
                         "8"="Stromal cell",
                         "9"="Stromal cell",
                         "10"="immune cell",
                         "11"="immune cell",
                         "12"="Stromal cell",
                         "13"="Stromal cell",
                         "14"="Stromal cell",
                         "15"="immune cell",
                         "16"="Stromal cell",
                         "17"="immune cell",
                         "18"="Stromal cell",
                         "19"="immune cell",
                         "20"="immune cell",
                         "21"="Stromal cell",
                         "22"="Stromal cell",
                         "23"="immune cell",
                         "24"="immune cell",
                         "25"="immune cell",
                         "26"="Epithelial",
                         "27"="Stromal cell",
                         "28"="Stromal cell",
                         "29"="immune cell",
                         "30"="Stromal cell",
                         "31"="immune cell",
                         "32"="Stromal cell",
                         "33"="immune cell" )

table(ovar2@meta.data$celltype.1)

Biocols=c("#9b409b","#d5cabd","#009d99","#a34c08","#e07f3e","#6a5fbe","#737bae",
          "#0077ce","#0088ca","#0095b6","#b6a6b5","#8487aa","#474b6b","#50434f")
p=DimPlot(ovar,reduction = "umap",label=T,group.by = "celltype.1",cols = Biocols)
ggsave("anno1.pdf",p,width=7,height=6)
DimPlot(ovar,reduction = "tsne",label=T,group.by = "celltype.1",cols = Biocols)
ggsave("anno1.1.pdf",width=7,height=6)

##第二次注释
markers <- c("IGHG1","MZB1","SDC1","CD19","CD79A","MS4A1", #B_Plasma cell
               "OGN","DCN","COL1A1", #Fibroblast
               "PECAM1","CLDN5","VWF", #Endothelia
               "CD4","CD14","C1QA", #Monocytic
               "KRT18","EPCAM","CD24","KRT19", #Epithelia
              "NCR1","GNLY","NKG7","CD3D","CD3E","CD8A") #T cell


#p=FeaturePlot(ovar,features = markers[1:9],ncol=3)
#ggsave("2.1_map.pdf",p,width=15,height=10)
p=DotPlot(ovar,features = markers) + RotatedAxis()
ggsave("2.2_Dot.pdf",p,width=14,height=6)

ovar$celltype.main<- recode(ovar@meta.data$seurat_clusters,
                         "0"="Fibroblast",
                         "1"="T cell",
                         "2"="Fibroblast",
                         "3"="T cell",
                         "4"="Monocytic",
                         "5"="Epithelia",
                         "6"="T cell",
                         "7"="T cell",
                         "8"="Fibroblast",
                         "9"="Endothelia",
                         "10"="B_Plasma cell",
                         "11"="T cell",
                         "12"="Fibroblast",
                         "13"="Fibroblast",
                         "14"="Fibroblast",
                         "15"="T cell",
                         "16"="Fibroblast",
                         "17"="Fibroblast",
                         "18"="Fibroblast",
                         "19"="B_Plasma cell",
                         "20"="Monocytic",
                         "21"="Fibroblast",
                         "22"="Fibroblast",
                         "23"="Monocytic",
                         "24"="T cell",
                         "25"="T cell",
                         "26"="Epithelia",
                         "27"="Fibroblast",
                         "28"="Fibroblast",
                         "29"="Monocytic",
                         "30"="Endothelia",
                         "31"="T cell",
                         "32"="Fibroblast",
                         "33"="B_Plasma cell" )

table(ovar@meta.data$celltype.main)

p=DimPlot(ovar,reduction = "umap",label=T,group.by = "celltype.main",cols = Biocols)
p
ggsave("22_anno_umap.pdf",p,width=7,height=6)


p=DimPlot(ovar,reduction = "tsne",label=T,group.by = "celltype.main",cols = Biocols)
ggsave("22_anno_tsne.pdf",p,width=7,height=6)
###cluster分组统计
tab.1=table(ovar@meta.data$seurat_clusters,ovar@meta.data$celltype.main)
balloonplot(tab.1)

###tissue分组统计
tab.2=table(ovar2@meta.data$tissue_type,ovar2@meta.data$celltype.main)
balloonplot(tab.2)

###样本分组统计
tab.3=table(ovar@meta.data$orig.ident,ovar2@meta.data$celltype.main)
balloonplot(tab.3,label.size=1)


table(Idents(ovar))
ovar1=ovar2
########备份数据
table(Idents(ovar1))
Idents(ovar1)="celltype.main"
p=DotPlot(ovar1,features = markers) + RotatedAxis()
ggsave("2.3_Dot.pdf",p,width=14,height=6)

###寻找每个细胞类型的markers
# ovar.markers3<-FindAllMarkers(ovar,only.pos = TRUE , logfc.threshold = 1 ,min.pct = 0.3)
# write.csv(ovar.markers3,file="markers.celltype.RNA.csv")
ovar.markers3=read.csv("F:/毕设/data3/质控前/markers.celltype.RNA.csv")

top10<-ovar.markers3%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
###对这些基因进行scaledata

markers=as.data.frame(top10[,"gene"])
ovar1=ScaleData(ovar1,features = as.character(unique(markers$gene)))


p = DoHeatmap(ovar2,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
ggsave("heatmap.pdf",p,width = 10,height = 9)

##每个细胞亚群抽取最少的细胞类型
allCells=names(Idents(ovar1))
allType=levels(Idents(ovar1))
choose_Cells=unlist(lapply(allType,function(x){
  cgCells=allCells[Idents(ovar1)==x]
  cg=sample(cgCells,min(table(ovar1@meta.data$celltype.main)))
}))

##提取
cg_sce=ovar1[,allCells %in% choose_Cells]
table(Idents(cg_sce))

p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
ggsave("heatmap1279.pdf",p,width = 10,height = 9)
ovar2=ovar
##堆叠柱状图
cell.prop<-as.data.frame(prop.table(table(ovar2@meta.data$celltype.main,ovar2@meta.data$tissue_type)))
colnames(cell.prop)<-c("cluster","group","proportion")

p<-ggplot(cell.prop,aes(group,proportion,fill = cluster))+
  geom_bar(stat="identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(title = NULL))

ggsave("./细胞注释/堆积图4.pdf",p,width = 14,height = 9)
p
saveRDS(ovar,"ovar_annotation.rds")
load("ovar_annotation.rds")
####第一次分型柱状图
cell.prop<-as.data.frame(prop.table(table(ovar2@meta.data$celltype.1,ovar2@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")

p<-ggplot(cell.prop,aes(group,proportion,fill = cluster))+
  geom_bar(stat="identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(title = NULL))
ggsave("./细胞注释/堆积图2.pdf",p,width = 14,height = 9)

###tissue堆积图
cell.prop<-as.data.frame(prop.table(table(ovar2@meta.data$celltype.1,ovar2@meta.data$tissue_type)))
colnames(cell.prop)<-c("cluster","group","proportion")

p<-ggplot(cell.prop,aes(group,proportion,fill = cluster))+
  geom_bar(stat="identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("./细胞注释/堆积图3.pdf",p,width = 14,height = 9)
#####

p<-ggplot(cell.prop,aes(group,proportion,fill = cluster))+
  geom_bar(stat="identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(title = NULL))+
  scale_fill_manual(values = c("Stromal cell" = "#9b409b", "immune cell" = "#d5cabd", "Epithelial" = "#009d99"))
####差异分析
Idents(ovar)
Idents(ovar)="celltype.main"
cell.markers=FindAllMarkers(ovar,
                            only.pos = FALSE,
                            test.use = "wilcox",
                            slot = "data",
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
write.csv(cell.markers,file="all.cell.markers.csv")
cell.markers=read.csv("./质控前/all.cell.markers.csv")
##根据p值进行筛选
# cell.markers.fil=cell.markers%>%filter(pct.1>0.5 & p_val_adj <0.2)%>%
#   filter(abs(avg_log2FC)>1)
# colnames(cell.markers.fil)
# table(abs(cell.markers.fil$avg_log2FC)>2)
# plotdt=cell.markers.fil%>%mutate(gene=ifelse(abs(avg_log2FC)>=2,gene,NA))
# ggplot(plotdt,aes(x=avg_log2FC,y=-log10(p_val_adj),
#                   size = pct.1,
#                   color=avg_log2FC))+
#   geom_point()+
#   ggtitle(label = "DCT")+
#   geom_text_repel(aes(label = gene),size=3,color="black")+
#   theme_bw()+
#   theme(plot.title = element_text(face = "bold",hjust = 0.5),
#         plot.background = element_rect(fill="transparent",colour = NA))+
#   scale_color_gradient2(low='olivedrab',high = 'salmon2',
#                         mid = 'grey',midpoint = 0)+
#   scale_size(range = c(1,3))

rownames(cell.markers) <- cell.markers$X
cell.markers=cell.markers[,-1]
ids=bitr(cell.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
cell.markers=merge(cell.markers,ids,by.x='gene',by.y='SYMBOL')
cell.markers=cell.markers[order(cell.markers$avg_log2FC,decreasing = T),]
cell.markers.list=as.numeric(cell.markers$avg_log2FC)
names(cell.markers.list)=cell.markers$ENTREZID
##GO
#筛选差异较大的基因集
ana.de=names(cell.markers.list)[abs(cell.markers.list)>1]
Go.ana=enrichGO(ana.de,OrgDb = "org.Hs.eg.db",ont="BP",readable = T)
dotplot(Go.ana,showCategory=10,title="GO-BP")
ggsave("./差异分析/GO.BP.png",width = 8,height = 7)
Go.ana2=enrichGO(ana.de,OrgDb = "org.Hs.eg.db",ont="MF",readable = T)
dotplot(Go.ana2,showCategory=10,title="GO-MF")
ggsave("./差异分析/GO.MF.png",width = 8,height = 7)
##KEGG
KEGG.ana=enrichKEGG(ana.de,organism = "hsa",pvalueCutoff = 0.05)
dotplot(KEGG.ana,showCategory=10,title="KEGG")
ggsave("./差异分析/KEGG.png",width = 8,height = 7)
##GSEA
GSEA.ana=gseKEGG(cell.markers.list,organism = "hsa",pvalueCutoff = 0.05)
GSEA.ana.arrange=arrange(GSEA.ana,desc(abs(NES)))

gsekp1=gseaplot2(GSEA.ana.arrange,1:5,color=Biocols,pvalue_table=F,base_size=14)
gsekp2=upsetplot(GSEA.ana.arrange,n=5)
cowplot::plot_grid(gsekp1,gsekp2,rel_widths = c(1,.6),labels = c("A","B"))
ggsave("./差异分析/GSEA.png",gsekp1+gsekp2,width = 14,height = 7)










