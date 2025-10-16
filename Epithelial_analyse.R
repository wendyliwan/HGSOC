####Epithelial cells
setwd("F:/毕设/data3/finial/F3----上皮细胞")
load("F:/毕设/data3/Epithelial/Epithelial_cells_annotation.Rdata")
Epithelial_ids <- which(ovar3@meta.data$celltype.main.2 == "Epithelia")
Epithelial_cells <- subset(ovar3, cells = Epithelial_ids)

DimPlot(Epithelial_cells,reduction = "umap",label=F,group.by = "tissue_type",cols = c("#db5a6b","#0077ce"))
ggsave("Epithelial_cells_umap.pdf",p,width=7,height=6)
DimPlot(Epithelial_cells,reduction = "tsne",label=F,group.by = "tissue_type",cols = c("#db5a6b","#0077ce"))
ggsave("Epithelial_cells_tsne.pdf",p,width=7,height=6)

Idents(Epithelial_cells)="celltype.main.2"
DimPlot(Epithelial_cells,reduction = "tsne",group.by = "celltype.main.2",label = T,repel = T,pt.size = 0.1)
save(Epithelial_cells,file = "./Epithelial/Epithelial_cells_annotation.Rdata")

Epi_sce=CreateSeuratObject(counts = GetAssayData(Epithelial_cells,assay = "RNA",layer = 'counts'),
                           meta.data = Epithelial_cells@meta.data)
Epi_sce=NormalizeData(Epi_sce)%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA(verbose = F)


Epi_sce=RunHarmony(Epi_sce,group.by.var="orig.ident",assay.use="RNA",max.iter.harmony=20)
ElbowPlot(Epi_sce,ndims = 50,reduction = "pca")
Epi_sce=FindNeighbors(Epi_sce,reduction = "pca",dims = 1:20)
Epi_sce=FindClusters(Epi_sce,resolution = 0.25)

table(Epi_sce@meta.data$seurat_clusters)
###0   1   2   3   4   5   6   7 
#653 416 359 153  77  34  32  15 
Epi_sce=RunUMAP(Epi_sce,reduction = "pca",dims = 1:20)%>%RunTSNE(dims=1:20,reduction = "pca")
col=pal_cosmic()(8)
DimPlot(Epi_sce,reduction = "umap",group.by="seurat_clusters",label = F,shuffle = F,raster = F,cols = col)
##p2=DimPlot(Epi_sce,reduction = "tsne",label = F,raster = F,cols = Biocols)
p1
ggsave(file="./Epithelial/Epi_umap.pdf",width = 6,height = 6)

DimPlot(Epi_sce,reduction = "umap",group.by = 'orig.ident',shuffle = F,raster = F,cols = col)
#DimPlot(Epi_sce,reduction = "tsne",group.by = 'orig.ident',raster = F,cols = Biocols)
ggsave(file="./Epithelial/Epi_orig_umap.pdf",width = 8,height = 6)

DimPlot(Epi_sce,reduction = "umap",split.by = "tissue_type",label = T,raster = F,cols = col)
ggsave(file="./Epithelial/Epi_tissue.pdf",width = 10,height = 6)
#DimPlot(Epi_sce,reduction = "tsne",split.by = "tissue_type",label = T,raster = F,cols = Biocols)

DimPlot(Epi_sce,reduction = "umap",group.by = "tissue_type",shuffle = T,raster = F,cols = c("#db5a6b","#0077ce"))
ggsave(file="./Epithelial/tissue_umap.pdf",width = 6,height = 6)
#DimPlot(Epi_sce,reduction = "tsne",group.by = "tissue_type",shuffle = T,raster = F,cols = Biocols)

DimPlot(Epi_sce,reduction = "umap",group.by = "tumor_stage",shuffle = T,raster = F,cols = col)


###寻找上皮细胞marker
##Epithelial_cells.markers1 <- FindMarkers(Epithelial_cells, ident.1 = 2, min.pct = 0.25)
write.csv(Epithelial_cells.markers1,file = "./Epithelial/Epithelial_cells.markers1.RNA.csv")

Epi_markers=FindAllMarkers(Epi_sce,test.use = "wilcox",
                           only.pos = T,
                           logfc.threshold = 0.25,
                           min.pct = 0.25)
write.csv(Epi_markers,file = "./Epithelial/Epi_markers.csv")

Epi.all.markers=Epi_markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)

Epi_top10=Epi.all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(Epi_top10,"./Epithelial/Epi_top10.csv")
###气泡图看marker
markers_Epi<-c("IGF2","BARX1",
               "SPINK1","TFF3","CST1",
               "PROK2","KIF20A","IGLC2","IGHG3",
               "MUC5B","CDH1",
               "FLG",
               "ITGAL","XCL2",
               "GRP","SFRP2")
DotPlot(Epi_sce,features = markers_Epi)+coord_flip()

FeaturePlot(Epi_sce,features = markers_Epi,ncol=3)
ggsave("./Epithelial/Epi_Features_umap.pdf",width=9,height=18)

VlnPlot(Epi_sce,features = markers_Epi,stack = T , flip = T)
##########CytoTRACE
Idents(Epi_sce)=Epi_sce$tumor_stage
exp1=as.matrix(GetAssayData(Epi_sce,assay = "RNA",layer = "counts"))
exp1=exp1[apply(exp1>0,1,sum)>=5,]
results=CytoTRACE(exp1,ncores = 1)
phenot1=Epi_sce$tumor_stage
phenot1=as.character(phenot1)
names(phenot1)=rownames(Epi_sce@meta.data)
emb=Epi_sce@reductions[["umap"]]@cell.embeddings
par(cex.main = 0.8)  # 调整主标题字体大小
plotCytoTRACE(results,phenotype = phenot1,emb = emb,outputDir = "./CytoTRACE")
plotCytoGenes(results,numOfGenes = 30,outputDir = "./CytoTRACE")

########Monocle2
DimPlot(Epi_sce,reduction="tsne",group.by = "tissue_type",label = T)+DimPlot(Epi_sce,reduction="tsne",group.by = "seurat_clusters",label = T)
head(Epi_sce@meta.data$celltype.1)
colnames(Epi_sce@meta.data)
data=GetAssayData(Epi_sce,assay = "RNA",slot = "counts")
pd=new('AnnotatedDataFrame',data=Epi_sce@meta.data[,c(1,6,28,29)])
fData=data.frame(gene_short_name=row.names(data),row.names=row.names(data))
fd=new('AnnotatedDataFrame',data=fData)
mycds=newCellDataSet(as.matrix(data),
                     phenoData =pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
mycds=estimateSizeFactors(mycds)
mycds=estimateDispersions(mycds,cores=8)

disp_table=dispersionTable(mycds)
order.genes=subset(disp_table,mean_expression>=0.005&dispersion_empirical>= +
                     1*dispersion_fit)%>%pull(gene_id)%>%as.character()
mycds=setOrderingFilter(mycds,order.genes)
plot_ordering_genes(mycds)

p=plot_ordering_genes(mycds)
ggsave("./Epithelial/拟时序/0OrderGenes.pdf",p,width = 8,height = 6)
mycds=reduceDimension(mycds,max_components = 2,reduction_method = 'DDRTree',
                      residualModelFormulaStr = "~orig.ident")

mycds=orderCells(mycds)
##若报错
trace('project2MST',edit=T,where=asNamespace("monocle"))
#
mycds=orderCells(mycds,root_state=1)

p1=plot_cell_trajectory(mycds,color_by = "State")+
  theme(legend.position = "right")
ggsave("./Epithelial/拟时序/1Trajectory_State.pdf",p1,width = 10,height = 6.5)

p2=plot_cell_trajectory(mycds,color_by = "Pseudotime")+
  theme(legend.position = "right")
ggsave("./Epithelial/拟时序/2Trajectory_Pseudotime.pdf",p2,width = 10,height = 6.5)

p3=plot_cell_trajectory(mycds,color_by = "tumor_stage")+
  theme(legend.position = "right")+
  scale_color_manual(values = col)
  
ggsave("./Epithelial/拟时序/3Trajectory_tissue.pdf",p3,width = 10,height = 6.5)

p4=plot_cell_trajectory(mycds,color_by = "orig.ident")+
  theme(legend.position = "right")+
  scale_color_manual(values = col)
pggsave("./Epithelial/拟时序/4Trajectory_Orig.pdf",p4,width = 10,height = 6.5)

p5=plot_cell_trajectory(mycds,color_by = "seurat_clusters",backbone_color = col)+
  scale_color_manual(values = col)+
  theme(legend.position = "right")
p2|p3|p5

ggplot(pData(mycds),aes(Pseudotime,colour=tumor_stage,fill=tumor_stage))+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()
ggsave("./Epithelial/拟时序/5Trajectory_density.pdf",p5,width = 10,height = 6.5)

genes=c(order.genes)[1:4]
genes=c("WSB1","IGF2","APOL1")
p1=plot_genes_in_pseudotime(mycds[genes],color_by = "seurat_clusters")
p2=plot_genes_in_pseudotime(mycds[genes],color_by = "Pseudotime")
p3=plot_genes_in_pseudotime(mycds[genes],color_by = "tumor_stage")
p1|p2|p3
ggsave("./Epithelial/拟时序/6Trajectory_pseudotime.pdf",p1|p2|p3,width = 10,height = 6.5)

p1=plot_genes_jitter(mycds[genes],grouping = "tumor_stage",color_by = "tumor_stage")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(legend.position = "none")
  
p2=plot_genes_violin(mycds[genes],grouping = "tumor_stage",color_by = "tumor_stage")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
p3=plot_genes_in_pseudotime(mycds[genes],color_by = "tumor_stage")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
p1+p2+p3
ggsave("./Epithelial/拟时序/7Trajectory_jitter.pdf",p1|p2|p3,width = 10,height = 6.5)

pData(mycds)$IGF2=log2(exprs(mycds)["IGF2",]+1)
p1=plot_cell_trajectory(mycds,color_by = "IGF2")+
  scale_color_continuous(type = "viridis")

pData(mycds)$APOL1=log2(exprs(mycds)["APOL1",]+1)
p2=plot_cell_trajectory(mycds,color_by = "APOL1")+
  scale_color_continuous(type = "viridis")

pData(mycds)$WSB1=log2(exprs(mycds)["WSB1",]+1)
p3=plot_cell_trajectory(mycds,color_by = "WSB1")+
  scale_color_continuous(type = "viridis")
p2+p1+p3
pData(mycds)$NOC2L=log2(exprs(mycds)["NOC2L",]+1)
p4=plot_cell_trajectory(mycds,color_by = "NOC2L")+
  scale_color_continuous(type = "viridis")
ggsave("./Epithelial/拟时序/8Trajectory_experssion.pdf",p1|p2|p3|p4,width = 10,height = 4)


###寻找拟时序差异基因——————monocle
Time_diff=differentialGeneTest(mycds,cores = 10,
                               fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"./Epithelial/拟时序/Time_diff.csv",row.names = F)
Time_diff=read.csv("./Epithelial/拟时序/Time_diff.csv")
Time_genes=Time_diff[order(Time_diff$qval),"gene_short_name"][1:100]
p=plot_pseudotime_heatmap(mycds[Time_genes,],num_clusters = 3,
                          show_rownames = T,return_heatmap = T)
ggsave("./Epithelial/拟时序/9Time_heatmap.pdf",p,width=5,height=10)
hp.genes=p$tree_row$labels[p$tree_row$order]
Time_diff_sig=Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"./Epithelial/拟时序/Time_diff_sig.csv",row.names = F)

####寻找拟时序差异基因——————beam
trace('project2MST',edit=T,where=asNamespace("monocle"))
beam_res=BEAM(mycds,branch_point = 2,cores = 10,progenitor_method="duplicate")

write.csv(beam_res,"Beam_all.csv",row.names = F)
beam_res=read.csv("Beam_all.csv")
beam_genes=beam_res[order(beam_res$qval),"gene_short_name"][1:100]
p=plot_genes_branched_heatmap(mycds[beam_genes,],branch_point = 2,num_clusters = 3,show_rownames = T,return_heatmap = T)
ggsave("Beam_heatmap.22.pdf",p$ph_res,width=5,height=7.5)

hp.genes=p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
Beam_sig=beam_res[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Beam_sig,"Beam_sig.csv",row.names = F)
head(Idents(Epi_sce))
Idents(Epi_sce)="tissue_type"
####上皮细胞差异分析
Epi_markers2=FindAllMarkers(Epi_sce,logfc.threshold = 0.25,
                            #min.pct = 0.2,
                            only.pos = FALSE,
                            #ident.1="Normal",ident.2="Cancer",
                            group.by="tissue_type",
                            )%>%mutate(gene=rownames(.))
jjVolcano(diffData=Epi_markers2,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          tile.col=Biocols,
          col.type="adjustP",
          topGeneN=10)

Idents(Epi_sce)="orig.ident"

Epi_markers3=FindAllMarkers(Epi_sce,logfc.threshold = 0.25,
                            #min.pct = 0.2,
                            only.pos = FALSE,
                            #ident.1="Normal",ident.2="Cancer",
                            group.by="orig.ident",
                            )%>%mutate(gene=rownames(.))
jjVolcano(diffData=Epi_markers3,
          log2FC.cutoff=0.25,
          size=3.5,
          fontface='italic',
          tile.col=Biocols,
          col.type="adjustP",
          topGeneN=5)
####富集分析

library(org.Hs.eg.db)
ids=bitr(Epi_markers3$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
Epi_markers3.1=merge(Epi_markers3,ids,by.x='gene',by.y='SYMBOL')

##KEGG
Epi_markers3.1=split(Epi_markers3.1$ENTREZID, Epi_markers3.1$cluster) 
## KEGG
xx <- compareCluster(Epi_markers3.1,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("./Epithelial/富集分析/kegg.pdf",width = 6,height = 10)

##############################cnv##############
setwd("F:/毕设/data3/finial/")
getwd()
DimPlot(ovar2,reduction = "umap",label=T,split.by = "tissue_type")
DimPlot(Epithelial_cells,reduction = "tsne",label=F,group.by = "orig.ident")

table(Epithelial_cells$seurat_clusters)
table(Epithelial_cells$tissue_type)
###抽样
###ovar2_part=ovar2[,sample(1:ncol(ovar2),round(ncol(ovar2)/5))]

Epi_matrix_counts=as.matrix(GetAssayData(Epithelial_cells,assay = "RNA" , layer = "counts"))
###创建新的分组
###ovar2$cnv_type=paste0(ovar2$tissue_type,"_",ovar2$seurat_clusters)
unique(Epithelial_cells$cnv_type)
####根据cluster分析
write.table(Epithelial_cells$cnv_type,"./Epithelial/Epi_cnv_type.txt",sep = "\t",quote = F , col.names = F)

gencode=read_tsv("hg38_gencode_v27.txt",col_names = c("gene","chr","start","end"))
gencode=gencode[!duplicated(gencode$gene),]
###提取基因交集
common_genes=intersect(gencode$gene,rownames(Epi_matrix_counts))
sort(unique(Epithelial_cells$cnv_type))####$后改


DimPlot(Epi_sce,reduction = "umap",label=T,split.by = "tissue_type")
DimPlot(Epi_sce,reduction = "umap",label=F,group.by = "orig.ident")

table(Epi_sce$seurat_clusters)
table(Epi_sce$tissue_type)
###抽样
###ovar2_part=ovar2[,sample(1:ncol(ovar2),round(ncol(ovar2)/5))]

Epi_matrix_counts=as.matrix(GetAssayData(Epi_sce,assay = "RNA" , layer = "counts"))
###创建新的分组
Epi_sce$Epi_cnv_type=paste0(Epi_sce$tissue_type,"_",Epi_sce$seurat_clusters)
unique(Epi_sce$Epi_cnv_type)
####根据cluster分析
write.table(Epi_sce$Epi_cnv_type,"Epi_cnv_type.txt",sep = "\t",quote = F , col.names = F)

gencode=read_tsv("F:/毕设/data3/hg38_gencode_v27.txt",col_names = c("gene","chr","start","end"))
gencode=gencode[!duplicated(gencode$gene),]
###提取基因交集
common_genes=intersect(gencode$gene,rownames(Epi_matrix_counts))
sort(unique(Epi_sce$Epi_cnv_type))####$后改

###创建Infercnv的object####cluster
infercnv_Epiobj=CreateInfercnvObject(raw_counts_matrix = Epi_matrix_counts[common_genes,],
                                     annotations_file = "Epi_cnv_type.txt",
                                     delim = "\t",
                                     gene_order_file = "F:/毕设/data3/hg38_gencode_v27.txt",
                                     ref_group_names = c("Normal_5" ,"Normal_2" ,"Normal_3"))


infercnv_obj=infercnv::run(infercnv_Epiobj,
                           cutoff = 0.1,
                           out_dir = "./cnv/",
                           no_prelim_plot = T,
                           cluster_by_groups = T,
                           denoise = T,
                           HMM = F,
                           min_cells_per_gene = 10,
                           num_threads = 10,
                           write_expr_matrix = T)

grp=read.table("./Epithelial/Epi_infercnv/infercnv.observation_groupings.txt",sep = "",header = T,row.names = NULL)
obs=data.table::fread("./Epithelial/Epi_infercnv/infercnv.observations.txt",data.table = F)%>%column_to_rownames(var = 'V1')

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
grp$cluster=factor(grp$cluster,levels=c("Cancer_0" ,"Cancer_1" ,"Cancer_2" ,"Cancer_3" ,"Cancer_4" ,"Cancer_6" ,"Cancer_7"))
row.names(grp)=row.names(cnv_score)
grp=grp[-1]
identical(row.names(cnv_score),row.names(grp))
cnv_score$cluster=grp$cluster
colnames(cnv_score)=c("Score","Cluster")
ggboxplot(cnv_score,"Cluster","Score",fill="Cluster")+
  coord_cartesian(ylim = c(1,1000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(cnv_score$Cluster)

ggplot(cnv_score,aes(x=Cluster,y=Score,fill=Cluster))+
  geom_violin()+
  coord_cartesian(ylim = c(0,500))+
  theme_minimal()+
  labs(x="Cluster",y="CNV Score")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))








bs=split(colnames(Epi_sce),Epi_sce$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(Epi_sce)%in%bs[[x]]
    data=GetAssayData(Epi_sce,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(Epi_sce@meta.data[,c('orig.ident',"tissue_type")])
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
write.csv(DEG_DESeq2,"EPI_DEG_DESeq2.csv")
table(DEG_DESeq2$Type)

p=ggplot(DEG_DESeq2,aes(log2FoldChange,-log10(pvalue)))+
  geom_point(size=3.5,
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
ggsave("./差异分析/Stro_volcano_unlabel.png",p,width = 8,height = 7)

top_20 <-bind_rows(
  DEG_DESeq2 %>%
    filter(Type == 'up') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(10),
  DEG_DESeq2 %>%
    filter(Type == 'down') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(10)
)
write.csv(top_20,"./差异分析/Stro_DEG_top10.csv")
rownames(top_20)=top_20$Gene_Sambol

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("./差异分析/Stro_volcano.png",volc_plot2,width = 8,height = 7)

##-------------------Go
middle=read.csv("MIDDLE.csv")
sig_DEG=subset(middle,Type!='ns')
genelist=bitr(sig_DEG$Gene_Sambol,fromType = "SYMBOL",
              toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
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
###选取感兴趣的基因集
select_ego1_res=ego1_res%>%
  dplyr::filter(grepl("pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))

ggplot(ego1_res[1:40,],aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  ggtitle("GO-BP")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
#-----------------------------kegg
kk1=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk1_res=kk1
kk1_res$Group="Cancer vs Normal"
kk1_res$LogP=-log(kk1_res$p.adjust)

ggplot(kk1_res[1:40,],aes(Count,Description))+
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







