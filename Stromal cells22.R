setwd("F:/毕设/data3/Stromal cell")
Stromal_cell.ids <- which(ovar3@meta.data$celltype.main.2 == c("Fibroblast","Endothelia"))
Stromal_cells <- subset(ovar3, cells = Stromal_cell.ids)
#########################
# test=Stromal_cells
# test=FindNeighbors(test,reduction = "pca",dims = 1:20)
# test=FindClusters(test,resolution = 0.05)
# table(test@meta.data$seurat_clusters)
# #    0    1    2    3    4    5    6    7    8    9 
# # 1605 1414 1304  672  646  635  295   70   58   42 
# test=RunUMAP(test,reduction = "pca",dims = 1:20,min.dist = 0.85)%>%RunTSNE(dims=1:20,reduction = "pca")
# DimPlot(test,reduction = "umap",label = T,shuffle = F,raster = F,cols = Biocols)
# DimPlot(test,reduction = "tsne",label = T,shuffle = F,raster = F,cols = Biocols)
# 
# DimPlot(test,reduction = "umap",group.by = 'orig.ident',shuffle = F,raster = F,cols = Biocols)
# 
# DimPlot(test,reduction = "umap",split.by = "tissue_type",label = T,raster = F,cols = Biocols)
# 
# DimPlot(test,reduction = "umap",group.by = "tissue_type",shuffle = T,raster = F,cols = Biocols)
# 
# markers <- c("PECAM1","CLDN5","EGFL7",###内皮细胞
#              "DCN","OGN","LUM","LAMB1",##上皮细胞
#              "ACTA2","MYH11","TAGLN","PDGFRB",###上皮-平滑肌细胞
#              "UPK3B")###间皮细胞
# 
# FeaturePlot(test,features = markers,ncol=3)
#######################
DimPlot(Stromal_cells,reduction = "umap",shuffle = F,label=F,group.by = "tissue_type",cols = c("#db5a6b","#0077ce"))
ggsave("Stromal_cellS_umap.pdf",width=7,height=6)
DimPlot(Stromal_cells,reduction = "tsne",label=F,shuffle = F,group.by = "tissue_type",cols = c("#db5a6b","#0077ce"))
ggsave("Stromal_cellS_tsne.pdf",width=7,height=6)

Stro_sce=CreateSeuratObject(counts = GetAssayData(Stromal_cells,assay = "RNA",layer = 'counts'),
                            meta.data = Stromal_cells@meta.data)
Stro_sce=NormalizeData(Stro_sce)%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA(verbose = F)
Stro_sce=RunHarmony(Stro_sce,group.by.var="orig.ident",assay.use="RNA",max.iter.harmony=20)
ElbowPlot(Stro_sce,ndims = 50,reduction = "pca")
Stro_sce=FindNeighbors(Stro_sce,reduction = "pca",dims = 1:20)
Stro_sce=FindClusters(Stro_sce,resolution = 0.05)
table(Stro_sce@meta.data$seurat_clusters)
#    0    1    2    3    4    5    6    7    8    9 
# 1605 1414 1304  672  646  635  295   70   58   42 
Stro_sce=RunUMAP(Stro_sce,reduction = "pca",dims = 1:20,min.dist = 0.85)%>%RunTSNE(dims=1:20,reduction = "pca")

setwd("F:/毕设/data3/Stromal cell")
setwd("F:/毕设/data3/Stromal cell/Stromalcell2")

cols=pal_igv()(12)
#"#5050FFFF" "#CE3D32FF" "#749B58FF" "#F0E685FF" "#466983FF" "#BA6338FF" "#5DB1DDFF" "#802268FF" "#6BD76BFF"
#[10] "#D595A7FF"
DimPlot(Stro_sce,reduction = "umap",label = T,shuffle = F,raster = F,cols = col)
ggsave(file="Stro_umap.pdf",width = 6,height = 6)
DimPlot(Stro_sce,reduction = "tsne",label = T,shuffle = F,raster = F,cols = col)
ggsave(file="Stro_tsne.pdf",width = 6,height = 6)

DimPlot(Stro_sce,reduction = "umap",group.by = 'orig.ident',shuffle = F,label=F,raster = F,cols = Biocols)
DimPlot(Stro_sce,reduction = "tsne",group.by = 'orig.ident',raster = F,cols = col)
ggsave(file="Stro_orig_umap.pdf",width = 8,height = 6)

DimPlot(Stro_sce,reduction = "umap",split.by = "tissue_type",label = T,raster = F,cols = Biocols)
ggsave(file="Stro_tissue.pdf",width = 10,height = 6)
#DimPlot(Stro_sce,reduction = "tsne",split.by = "tissue_type",label = T,raster = F,cols = Biocols)

DimPlot(Stro_sce,reduction = "umap",group.by = "tissue_type",shuffle = T,raster = F,cols = c("#db5a6b","#0077ce"))
ggsave(file="tissue_umap.pdf",width = 6,height = 6)
DimPlot(Stro_sce,reduction = "tsne",group.by = "tissue_type",shuffle = T,raster = F,cols = Biocols)

####样本marker
Stro_markers=FindAllMarkers(Stro_sce,test.use = "wilcox",
                            only.pos = T,
                            logfc.threshold = 0.25,
                            min.pct = 0.25)
write.csv(Stro_markers,file = "Stro_markers.csv")
Stro_markers=read.csv("F:/毕设/data3/Stromal cell/markers/Stro_markers.csv")
Stro.all.markers=Stro_markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)

Stro_top10=Stro.all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(Stro_top10,"./markers/Stro_top10.csv")

# Stro_markers2=FindAllMarkers(Stro_sce,
#                             only.pos = FALSE,
#                             test.use = "wilcox",
#                             slot = "data",
#                             min.pct = 0.25,
#                             logfc.threshold = 0.25)
# write.csv(Stro_markers2,file = "Stro_markers2.csv")
#寻找样本marker
Stro.markers1 <- FindMarkers(Stro_sce, ident.1 = 2, min.pct = 0.25)
head(Stro.markers1)
write.csv(Stro.markers1,file = "Stro.markers1.csv")
Stro.markers1=read.csv("F:/毕设/data3/Stromal cell/markers/Stro.markers1.orig.csv")

Stro.markers2 <- FindMarkers(Stro_sce, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
write.csv(Stro.markers2,file = "Stro.markers2.csv")
Stro.markers2=read.csv("Stro.markers2.csv")

markers <- c("PECAM1","CLDN5","EGFL7",###内皮细胞
             "DCN","OGN","LAMB1",##上皮细胞
             "ACTA2","MYH11","TAGLN",###上皮-平滑肌细胞
             "UPK3B")###间皮细胞
markers=c("PECAM1","DCN","MYH11","UPK3B")
FeaturePlot(Stro_sce,features = markers,ncol=4)
ggsave("1markers_map.pdf",width=15,height=20)

DotPlot(Stro_sce,features = markers) + RotatedAxis()+ 
  theme(legend.position = "none")
ggsave("2markers_Dot.pdf",width=14,height=6)

VlnPlot(Stro_sce,features = markers,stack = T , flip = T)
ggsave("3_Vln.pdf",width=14,height=6)

DimPlot(Stro_sce,reduction = "umap",label=T)
DimPlot(Stro_sce,reduction = "tsne",label=T)
#ggsave("./细胞注释/第一次注释/4_tsne.pdf",width=7,height=6)

Stro_sce$subcelltype<- recode(Stro_sce@meta.data$seurat_clusters,
                              "0"="Mesothelial",
                              "1"="Mesothelial",
                              "2"="Fibroblast",
                              "3"="Fibroblast",
                              "4"="Fibroblast",
                              "5"="Fibroblast",
                              "6"="myoFibroblast",
                              "7"="Fibroblast",
                              "8"="Mesothelial",
                              "9"="Endothelial")

table(Stro_sce@meta.data$subcelltype)
col=pal_d3()(4)
DimPlot(Stro_sce,reduction = "umap",label=F,group.by = "subcelltype",cols = col)+
 labs(title = "Stromal cells")
ggsave("Stro_anno_umap.pdf",width=7,height=6)

DimPlot(Stro_sce,reduction = "umap",label=F,group.by = "subcelltype",split.by="tissue_type",cols = col)+
  labs(title = "Stromal cells")
ggsave("Stro_anno_umap_tissue.pdf",width=12,height=6)

DimPlot(Stro_sce,reduction = "tsne",label=T,group.by = "subcelltype",cols = col)+
  labs(title = "Stromal cells")
ggsave("Stro_anno_tsne.pdf",width=7,height=6)

DimPlot(Stro_sce,reduction = "tsne",label=F,group.by = "subcelltype",split.by="tissue_type",cols = col)
ggsave("Stro_anno_tsne_tissue.pdf",width=12,height=6)

####细胞类型marker
table(Idents(Stro_sce))
Idents(Stro_sce)="subcelltype"
DotPlot(Stro_sce,features = markers) + RotatedAxis()
ggsave("2.3_Dot.pdf",p,width=14,height=6)
Stro.markers3<-FindAllMarkers(Stro_sce,only.pos = TRUE , logfc.threshold = 1 ,min.pct = 0.3)
write.csv(Stro.markers3,file="./markers/Stro.markers.celltype.csv")
Stro.markers3=read.csv("F:/毕设/data3/Stromal cell/markers/Stro.markers.celltype.csv")
###备份一下
Stro_sce1=Stro_sce
top10<-Stro.markers3%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
###对这些基因进行scaledata
markers=as.data.frame(top10[,"gene"])
Stro_sce1=ScaleData(Stro_sce1,features = as.character(unique(markers$gene)))
DoHeatmap(Stro_sce1,
              features = as.character(unique(markers$gene)),
              group.by = "subcelltype")+
  scale_fill_viridis()
DoHeatmap(Stro_sce1,
          features = markers$gene,
          group.by = "subcelltype", draw.lines = F)+
  scale_fill_viridis()
ggsave("heatmap.pdf",p,width = 10,height = 9)

##每个细胞亚群抽取最少的细胞类型
allCells=names(Idents(Stro_sce1))
allType=levels(Idents(Stro_sce1))
choose_Cells=unlist(lapply(allType,function(x){
  cgCells=allCells[Idents(Stro_sce1)==x]
  cg=sample(cgCells,min(table(Stro_sce1@meta.data$subcelltype)))
}))

##提取
cg_sce=Stro_sce1[,allCells %in% choose_Cells]
table(Idents(cg_sce))

DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "subcelltype")+
  scale_fill_viridis()
p+theme(axis.text.x = element_text(angle = 0, size = 12))
ggsave("heatmap42.pdf",p,width = 10,height = 9)

#########堆积图
cell.prop<-as.data.frame(prop.table(table(Stro_sce1@meta.data$subcelltype,Stro_sce1@meta.data$tissue_type)))
colnames(cell.prop)<-c("cluster","group","proportion")

p<-ggplot(cell.prop,aes(group,proportion,fill = cluster))+
  geom_bar(stat="identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(title = NULL))
ggsave("堆积图1.pdf",p,width = 14,height = 9)

#######################CytoTRACE#############
exp1=as.matrix(GetAssayData(Stro_sce,assay = "RNA",layer = "counts"))
exp1=exp1[apply(exp1>0,1,sum)>=5,]
results=CytoTRACE(exp1,ncores = 1)
phenot1=Stro_sce$subcelltype
phenot1=as.character(phenot1)
names(phenot1)=rownames(Stro_sce@meta.data)
emb=Stro_sce@reductions[["umap"]]@cell.embeddings
par(cex.main = 0.8)  # 调整主标题字体大小
plotCytoTRACE(results,phenotype = phenot1,emb = emb,outputDir = "./CytoTRACE")
plotCytoGenes(results,numOfGenes = 30,outputDir = "./CytoTRACE")
save(Stro_sce,file = "Stromal cells.Rdata")
load("Stromal cells.Rdata")


#######################Monocle2################
DimPlot(Stro_sce,reduction="umap",group.by = "subcelltype",label = T,cols = col)+
  DimPlot(Stro_sce,reduction="umap",group.by = "seurat_clusters",label = T,cols = col)
ggsave("./注释/anno+clusters.pdf",width = 14,height = 6)
Idents(Stro_sce)=Stro_sce$subcelltype
table(Stro_sce$subcelltype)

colnames(Stro_sce@meta.data)
data=GetAssayData(Stro_sce,assay = "RNA",slot = "counts")###=exprmatrix
data=as.matrix(data)
# exp1=as.matrix(GetAssayData(Stro_sce,assay = "RNA",layer = "counts"))
pd=new('AnnotatedDataFrame',data=Stro_sce@meta.data[,c(1,6,27,32)])
fData=data.frame(gene_short_name=row.names(data),row.names=row.names(data))
fd=new('AnnotatedDataFrame',data=fData)
mycds=newCellDataSet(as.matrix(data),
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
#ggsave("./拟时序/0OrderGenes.pdf",p,width = 8,height = 6)
mycds=reduceDimension(mycds,max_components = 2,reduction_method = 'DDRTree',
                      residualModelFormulaStr = "~orig.ident")

mycds=orderCells(mycds)
##若报错
trace('project2MST',edit=T,where=asNamespace("monocle"))
#
##mycds=orderCells(mycds,root_state=5)

plot_cell_trajectory(mycds,color_by = "State")
ggsave("./拟时序/1Trajectory_State.pdf",width = 10,height = 6.5)

plot_cell_trajectory(mycds,color_by = "Pseudotime")
ggsave("./拟时序/2Trajectory_Pseudotime.pdf",p2,width = 10,height = 6.5)

plot_cell_trajectory(mycds,color_by = "subcelltype")
ggsave("./拟时序/3Trajectory_tissue.pdf",p3,width = 10,height = 6.5)

plot_cell_trajectory(mycds,color_by = "orig.ident")
ggsave("./拟时序/4Trajectory_Orig.pdf",width = 10,height = 6.5)

ggplot(pData(mycds),aes(Pseudotime,colour=subcelltype,fill=subcelltype))+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()
ggsave("./拟时序/5Trajectory_density.pdf",width = 10,height = 6.5)

genes=c("PECAM1","DCN","MYH11","UPK3B")
p1=plot_genes_in_pseudotime(mycds[genes],color_by = "State")
p1=plot_genes_in_pseudotime(mycds[genes],color_by = "Pseudotime")
p2=plot_genes_in_pseudotime(mycds[genes],color_by = "subcelltype")
p1|p2|p3
ggsave("./拟时序/6Trajectory_pseudotime.pdf",p1|p2|p3,width = 10,height = 6.5)

p1=plot_genes_jitter(mycds[genes],grouping = "State",color_by = "State")
p2=plot_genes_violin(mycds[genes],grouping = "State",color_by = "State")
p3=plot_genes_in_pseudotime(mycds[genes],color_by = "State")
ggsave("./拟时序/7Trajectory_jitter.pdf",p1|p2|p3,width = 10,height = 6.5)

pData(mycds)$PECAM1=log2(exprs(mycds)["PECAM1",]+1)
p1=plot_cell_trajectory(mycds,color_by = "PECAM1")+
  scale_color_continuous(type = "viridis")

pData(mycds)$DCN=log2(exprs(mycds)["DCN",]+1)
p2=plot_cell_trajectory(mycds,color_by = "DCN")+
  scale_color_continuous(type = "viridis")

pData(mycds)$MYH11=log2(exprs(mycds)["MYH11",]+1)
p3=plot_cell_trajectory(mycds,color_by = "MYH11")+
  scale_color_continuous(type = "viridis")

pData(mycds)$UPK3B=log2(exprs(mycds)["UPK3B",]+1)
p4=plot_cell_trajectory(mycds,color_by = "UPK3B")+
  scale_color_continuous(type = "viridis")
ggsave("./拟时序/8Trajectory_experssion.pdf",p1|p2|p3|p4,width = 10,height = 4)

###寻找拟时序差异基因——————monocle
Time_diff=differentialGeneTest(mycds,cores = 10,
                               fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"./拟时序/DEG--monocle/Time_diff.csv",row.names = F)
Time_diff=read.csv("./拟时序/DEG--monocle/Time_diff.csv")

Time_genes=Time_diff[order(Time_diff$qval),"gene_short_name"][1:100]
plot_pseudotime_heatmap(mycds[Time_genes,],num_clusters = 3,
                          show_rownames = T,return_heatmap = T)
ggsave("./拟时序/DEG--monocle/9Time_heatmap.pdf",p,width=5,height=10)
hp.genes=p$tree_row$labels[p$tree_row$order]
Time_diff_sig=Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"./拟时序/DEG--monocle/Time_diff_sig.csv",row.names = F)

####寻找拟时序差异基因——————beam
beam_res=BEAM(mycds,branch_point = 2,cores = 10,progenitor_method="duplicate")
write.csv(beam_res,"./拟时序/DEG--BEAM/Beam_all.csv",row.names = F)
beam_res=read.csv("F:/毕设/data3/Stromal cell/拟时序/DEG--BEAM/Beam_all.csv")
beam_genes=beam_res[order(beam_res$qval),"gene_short_name"][1:100]
plot_genes_branched_heatmap(mycds[beam_genes,],branch_point = 1,num_clusters = 3,show_rownames = T,return_heatmap = T)
ggsave("./拟时序/DEG--BEAM/10Beam_heatmap.pdf",p$ph_res,width=6.5,height=10)

hp.genes=p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
Beam_sig=beam_res[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Beam_sig,"F:/毕设/data3/Stromal cell/markers/拟时序/DEG--BEAM/Beam_sig.csv",row.names = F)

Idents(Stro_sce)



##################差异分析celltype############
Stro_subcell_markers=FindAllMarkers(Stro_sce,logfc.threshold = 0.25,
                                    min.pct = 0.2,
                                    only.pos = FALSE,
                                    # ident.1="Fibroblast",ident.2="Endothelial",
                                    # ident.3="myoFibroblast",ident.4="Mesothelial",
                                    group.by="subcelltype",
)%>%mutate(gene=rownames(.))
write.csv(Stro_subcell_markers,"./差异分析/Stro_subcell_markers.csv",row.names = F)
Stro_subcell_markers=read.csv("F:/毕设/data3/Stromal cell/差异分析/Stro_subcell_markers.csv")
jjVolcano(diffData=Stro_subcell_markers,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          tile.col=col,
          col.type="adjustP",
          topGeneN=10,
          flip=T)
ggsave("火山图差异分析all.png",width=22,height=15)
ggsave("./差异分析/subcell-volcane.pdf",width = 22,height = 15)

# 创建一个逻辑向量来标识不包含小数点的基因名
no_decimal_genes <- !grepl("\\.", Stro_subcell_markers$gene)

# 过滤掉包含小数点的基因名
filtered_markers <- Stro_subcell_markers[no_decimal_genes, ]

Stro_subcell_markers_2=read.csv("F:/毕设/data3/Stromal cell/差异分析/Stro_subcell_markers.csv")
rownames(Stro_subcell_markers_2) <- Stro_subcell_markers_2$gene
deg_top20=filtered_markers%>%
  dplyr::group_by(cluster)%>%
  dplyr::top_n(n=20,wt=avg_log2FC)
col<-pal_npg()(4)
#Stro_sce1=JoinLayers(Stro_sce)
plot=averageHeatmap(Stro_sce1,
                    markerGene = deg_top20$gene,
                    group.by = "subcelltype",
                    gene.order = deg_top20$gene,
                    annoCol = T,
                    myanCol = col
)
ggsave("./差异分析/DEG-Heatmap.pdf",plot,width = 10,height = 15)
jjVolcano(diffData=filtered_markers,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          tile.col=col,
          col.type="adjustP",
          topGeneN=10,legend.position = c(0.9,0.9))
ggsave("筛选后volcano.png",width = 22,height = 15)
# #根据p值进行筛选
# Stro_subcell_markers_fil=Stro_subcell_markers%>%filter(pct.1>0.5 & p_val_adj <0.2)%>%
#   filter(abs(avg_log2FC)>1)
# colnames(Stro_subcell_markers_fil)
# table(abs(Stro_subcell_markers_fil$avg_log2FC)>2)
# plotdt=Stro_subcell_markers_fil%>%mutate(gene=ifelse(abs(avg_log2FC)>=2,gene,NA))
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

# rownames(Stro_subcell_markers) <- Stro_subcell_markers$gene
# ids=bitr(Stro_subcell_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
# Stro_subcell_markers=merge(Stro_subcell_markers,ids,by.x='gene',by.y='SYMBOL')
# Stro_subcell_markers=Stro_subcell_markers[order(Stro_subcell_markers$avg_log2FC,decreasing = T),]
# Stro_subcell_markers.list=as.numeric(Stro_subcell_markers$avg_log2FC)
# names(Stro_subcell_markers.list)=Stro_subcell_markers$ENTREZID
# Stro_ana.de=names(Stro_subcell_markers.list)[abs(Stro_subcell_markers.list)>1]


rownames(filtered_markers) <- filtered_markers$gene
ids=bitr(filtered_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Stro_subcell_markers_dif=merge(filtered_markers,ids,by.x='gene',by.y='SYMBOL')
Stro_subcell_markers_dif=Stro_subcell_markers_dif[order(filtered_markers$avg_log2FC,decreasing = T),]
Stro_subcell_markers.list=as.numeric(Stro_subcell_markers_dif$avg_log2FC)
names(Stro_subcell_markers.list)=Stro_subcell_markers_dif$ENTREZID
Stro_ana.de=names(Stro_subcell_markers.list)[abs(Stro_subcell_markers.list)>1]



##GO
Go.ana=enrichGO(Stro_ana.de,OrgDb = "org.Hs.eg.db",ont="BP",readable = T)
dotplot(Go.ana,showCategory=10,title="GO-BP")
ggsave("./差异分析/GO.BP.png",width = 8,height = 7)
Go.ana2=enrichGO(Stro_ana.de,OrgDb = "org.Hs.eg.db",ont="MF",readable = T)
dotplot(Go.ana2,showCategory=10,title="GO-MF")
Go.ana=enrichGO(Stro_ana.de,OrgDb = "org.Hs.eg.db",ont="ALL",readable = T)
dotplot(Go.ana,showCategory=10,title="GO-ALL")
ggsave("./差异分析/GO.MF.png",width = 8,height = 7)
##KEGG
KEGG.ana=enrichKEGG(Stro_ana.de,organism = "hsa",pvalueCutoff = 0.05)
dotplot(KEGG.ana,showCategory=10,title="KEGG")
ggsave("./差异分析/KEGG.png",width = 8,height = 7)
##GSEA
Stro_subcell_geneList <- sort(Stro_subcell_markers.list, decreasing = TRUE)
GSEA.ana=gseKEGG(Stro_subcell_geneList,organism = "hsa",pvalueCutoff = 0.05)
GSEA.ana.arrange=arrange(GSEA.ana,desc(abs(NES)))

gsekp1=gseaplot2(GSEA.ana.arrange,1:5,color=cols,pvalue_table=F,base_size=14)
gsekp2=upsetplot(GSEA.ana.arrange,n=5)
cowplot::plot_grid(gsekp1,gsekp2,rel_widths = c(1,.6),labels = c("A","B"))
ggsave("./差异分析/GSEA.png",gsekp1|gsekp2,width = 14,height = 7)






#####################CNV##########
DimPlot(Stro_sce,reduction = "umap",label=T,split.by = "tissue_type")
DimPlot(Stro_sce,reduction = "umap",label=F,group.by = "orig.ident")

table(Stro_sce$seurat_clusters)
table(Stro_sce$tissue_type)
###抽样
###ovar2_part=ovar2[,sample(1:ncol(ovar2),round(ncol(ovar2)/5))]

Stro_matrix_counts=as.matrix(GetAssayData(Stro_sce,assay = "RNA" , layer = "counts"))
###创建新的分组
Stro_sce$Stro_cnv_type=paste0(Stro_sce$tissue_type,"_",Stro_sce$subcelltype)
unique(Stro_sce$Stro_cnv_type)
####根据cluster分析
write.table(Stro_sce$Stro_cnv_type,"./infercnv/2_Stro_cnv_type.txt",sep = "\t",quote = F , col.names = F)

gencode=read_tsv("F:/毕设/data3/hg38_gencode_v27.txt",col_names = c("gene","chr","start","end"))
gencode=gencode[!duplicated(gencode$gene),]
###提取基因交集
common_genes=intersect(gencode$gene,rownames(Stro_matrix_counts))
sort(unique(Stro_sce$Stro_cnv_type))####$后改

###创建Infercnv的object####cluster
infercnv_Stroobj=CreateInfercnvObject(raw_counts_matrix = Stro_matrix_counts[common_genes,],
                                      annotations_file = "./infercnv/2_Stro_cnv_type.txt",
                                      delim = "\t",
                                      gene_order_file = "F:/毕设/data3/hg38_gencode_v27.txt",
                                      ref_group_names = c("Normal_Endothelial","Normal_Fibroblast" ,"Normal_Mesothelial","Normal_myoFibroblast"))
save(infercnv_Stroobj,file = "./infercnv/infercnv_Stroobj.Rdata")
load("infercnv_Stroobj.Rdata")

infercnv_obj=infercnv::run(infercnv_Stroobj,
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
load("infercnv_Stroobj.Rdata")

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
cnv_score<- na.omit(cnv_score)

grp$cluster=sapply(strsplit(grp$Dendrogram.Group,"_"),function(x) paste(x[1],x[2],sep = "_"))
grp<- na.omit(grp)
sort(unique(grp$cluster))
grp$cluster=factor(grp$cluster,levels=c("Cancer_Endothelial","Cancer_Fibroblast","Cancer_myoFibroblast"))
rownames(grp)=grp[,1]
grp=grp[-1]

cnv_score$row.name=rownames(cnv_score)
grp$row.name=rownames(grp)
# 执行内连接，只保留共同的行名
cnv_score_clean <- inner_join(cnv_score, grp, by = "row.name")
cnv_score_clean=cnv_score_clean[,-c(3:7)]
rownames(cnv_score_clean)=cnv_score_clean[,2]
cnv_score_clean=cnv_score_clean[1]


identical(row.names(cnv_score_clean),row.names(grp))
cnv_score_clean$cluster=grp$cluster
colnames(cnv_score_clean)=c("Score","Cluster")
ggboxplot(cnv_score_clean,"Cluster","Score",fill="Cluster")+
  coord_cartesian(ylim = c(1,1000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(cnv_score_clean$Cluster)

ggplot(cnv_score_clean,aes(x=Cluster,y=Score,fill=Cluster))+
  geom_violin()+
  coord_cartesian(ylim = c(0,500))+
  theme_minimal()+
  labs(x="Cluster",y="CNV Score")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))








#################差异分析（正常vs疾病）#####
table(Idents(Stro_sce))
Idents(Stro_sce)="tissue_type"
Stro_subcell_markers=FindAllMarkers(Stro_sce,logfc.threshold = 0.25,
                                    min.pct = 0.2,
                                    only.pos = FALSE,
                                    # ident.1="Fibroblast",ident.2="Endothelial",
                                    # ident.3="myoFibroblast",ident.4="Mesothelial",
                                    group.by="tissue_type",
)%>%mutate(gene=rownames(.))
write.csv(Stro_subcell_markers,"F:/毕设/data3/Stromal cell/差异分析/Stro_subcell_markers_nc.csv",row.names = F)
#Stro_subcell_markers=read.csv("F:/毕设/data3/Stromal cell/差异分析/Stro_subcell_markers.csv")
jjVolcano(diffData=Stro_subcell_markers,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          tile.col=c("#db5a6b","#0077ce"),
          col.type="adjustP",
          topGeneN=10,legend.position = c(0.9,0.95))
ggsave("火山图差异分析all.png",width=22,height=15)
ggsave("./差异分析/subcell-volcane.pdf",width = 22,height = 15)

# 创建一个逻辑向量来标识不包含小数点的基因名
no_decimal_genes <- !grepl("\\.", Stro_subcell_markers$gene)

# 过滤掉包含小数点的基因名
filtered_markers <- Stro_subcell_markers[no_decimal_genes, ]

# Stro_subcell_markers_fil=read.csv("F:/毕设/data3/Stromal cell/差异分析/Stro_subcell_markers.csv")
# rownames(Stro_subcell_markers_2) <- Stro_subcell_markers_2$gene
deg_top20=Stro_subcell_markers%>%
  dplyr::group_by(cluster)%>%
  dplyr::top_n(n=20,wt=avg_log2FC)
col<-pal_npg()(4)
#Stro_sce1=JoinLayers(Stro_sce)
plot=averageHeatmap(Stro_sce1,
                    markerGene = deg_top20$gene,
                    group.by = "tissue_type",
                    gene.order = deg_top20$gene,
                    annoCol = T,
                    myanCol = c("#db5a6b","#0077ce")
)
ggsave("./差异分析/DEG-Heatmap.pdf",plot,width = 10,height = 15)
jjVolcano(diffData=filtered_markers,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          tile.col=c("#db5a6b","#0077ce"),
          col.type="adjustP",
          topGeneN=10,legend.position = c(0.9,0.9))
ggsave("筛选后volcano.png",width = 22,height = 15)
# #根据p值进行筛选
# Stro_subcell_markers_fil=Stro_subcell_markers%>%filter(pct.1>0.5 & p_val_adj <0.2)%>%
#   filter(abs(avg_log2FC)>1)
# colnames(Stro_subcell_markers_fil)
# table(abs(Stro_subcell_markers_fil$avg_log2FC)>2)
# plotdt=Stro_subcell_markers_fil%>%mutate(gene=ifelse(abs(avg_log2FC)>=2,gene,NA))
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

# rownames(Stro_subcell_markers) <- Stro_subcell_markers$gene
# ids=bitr(Stro_subcell_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
# Stro_subcell_markers=merge(Stro_subcell_markers,ids,by.x='gene',by.y='SYMBOL')
# Stro_subcell_markers=Stro_subcell_markers[order(Stro_subcell_markers$avg_log2FC,decreasing = T),]
# Stro_subcell_markers.list=as.numeric(Stro_subcell_markers$avg_log2FC)
# names(Stro_subcell_markers.list)=Stro_subcell_markers$ENTREZID
# Stro_ana.de=names(Stro_subcell_markers.list)[abs(Stro_subcell_markers.list)>1]


rownames(filtered_markers) <- filtered_markers$gene
ids=bitr(filtered_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Stro_subcell_markers_dif=merge(filtered_markers,ids,by.x='gene',by.y='SYMBOL')
Stro_subcell_markers_dif=Stro_subcell_markers_dif[order(filtered_markers$avg_log2FC,decreasing = T),]
Stro_subcell_markers.list=as.numeric(Stro_subcell_markers_dif$avg_log2FC)
names(Stro_subcell_markers.list)=Stro_subcell_markers_dif$ENTREZID
Stro_ana.de=names(Stro_subcell_markers.list)[abs(Stro_subcell_markers.list)>1]



##GO
Go.ana=enrichGO(Stro_ana.de,OrgDb = "org.Hs.eg.db",ont="BP",readable = T)
dotplot(Go.ana,showCategory=10,title="GO-BP")
ggsave("./差异分析/GO.BP.png",width = 8,height = 7)
Go.ana2=enrichGO(Stro_ana.de,OrgDb = "org.Hs.eg.db",ont="MF",readable = T)
dotplot(Go.ana2,showCategory=10,title="GO-MF")
Go.ana=enrichGO(Stro_ana.de,OrgDb = "org.Hs.eg.db",ont="ALL",readable = T)
dotplot(Go.ana,showCategory=10,title="GO-ALL")
ggsave("./差异分析/GO.MF.png",width = 8,height = 7)
##KEGG
KEGG.ana=enrichKEGG(Stro_ana.de,organism = "hsa",pvalueCutoff = 0.05)
dotplot(KEGG.ana,showCategory=10,title="KEGG")
ggsave("./差异分析/KEGG.png",width = 8,height = 7)
##GSEA
Stro_subcell_geneList <- sort(Stro_subcell_markers.list, decreasing = TRUE)
GSEA.ana=gseKEGG(Stro_subcell_geneList,organism = "hsa",pvalueCutoff = 0.05)
GSEA.ana.arrange=arrange(GSEA.ana,desc(abs(NES)))

gsekp1=gseaplot2(GSEA.ana.arrange,1:5,color=cols,pvalue_table=F,base_size=14)
gsekp2=upsetplot(GSEA.ana.arrange,n=5)
cowplot::plot_grid(gsekp1,gsekp2,rel_widths = c(1,.6),labels = c("A","B"))
ggsave("./差异分析/GSEA.png",gsekp1|gsekp2,width = 14,height = 7)

######################cellchat######
data1=GetAssayData(Stro_sce,assay = "RNA",layer = "data")
table(Stro_sce@meta.data$subcelltype)
table(Stro_sce@meta.data$orig.ident)
cellchat=createCellChat(data1,meta = Stro_sce@meta.data,
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