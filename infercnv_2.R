ovar2=ovar
ovar2$tissue_type<- recode(ovar2@meta.data$orig.ident,
                           "GSM5599220_Norm1"="Normal",
                           "GSM5599221_Norm2"="Normal",
                           "GSM5599222_Norm3"="Normal",
                           "GSM5599223_Norm4"="Normal",
                           "GSM5599224_Norm5"="Normal",
                           "GSM5599225_Cancer1"="Cancer",
                           "GSM5599226_Cancer2"="Cancer",
                           "GSM5599227_Cancer3"="Cancer",
                           "GSM5599228_Cancer4"="Cancer",
                           "GSM5599229_Cancer5"="Cancer",
                           "GSM5599230_Cancer6"="Cancer",
                           "GSM5599231_Cancer7"="Cancer")
p=DimPlot(ovar2,reduction = "umap",label=T,group.by = "tissue_type",cols = Biocols)
ggsave("Tissue_umap.pdf",p,width=7,height=6)
p

p=DimPlot(ovar2,reduction = "tsne",label=T,shuffle=T,group.by = "tissue_type",cols = Biocols)
ggsave("Tissue_tsne.pdf",p,width=7,height=6)
p

save(ovar2,file = "Ovar_fin.Rdata")
load("Ovar_fin.Rdata")
normal_ids <- which(ovar2@meta.data$tissue_type == "Normal")
cancer_ids2 <- which(ovar2@meta.data$tissue_type == "Cancer")
# 使用Seurat的Subset函数来创建一个新的Seurat对象，只包含normal_ids指定的细胞
ovar_Normal <- subset(ovar2, cells = normal_ids)
ovar_Cancer <- subset(ovar2, cells = cancer_ids2)

save(ovar_Normal,file = "ovar_Normal.Rdata")
load("ovar_Normal.Rdata")

save(ovar_Cancer,file = "ovar_Cancer.Rdata")
load("ovar_Cancer.Rdata")

p1=DimPlot(ovar_Normal,reduction = "tsne",label=T,group.by = "celltype.main",cols = Biocols)+
  labs(title = "Normal")+ 
  theme(legend.position = "none")
ggsave("Normal_tsne.pdf",p1,width=7,height=6)

p2=DimPlot(ovar_Cancer,reduction = "tsne",label=T,group.by = "celltype.main",cols = Biocols)+
  labs(title = "Cancer")
ggsave("Cancer_tsne.pdf",p2,width=7,height=6)

p1+p2
ggsave("Normal_Cancer_tsne.pdf",p1+p2,width=14,height=6)


##############################################拷贝数变异#########
BiocManager::install("infercnv")
DimPlot(ovar2,reduction = "umap",label=T,split.by = "tissue_type")
DimPlot(ovar2,reduction = "tsne",label=T,group.by = "orig.ident")

table(ovar3$seurat_clusters)
table(ovar3$tissue_type)
###抽样
ovar2_part=ovar2[,sample(1:ncol(ovar2),round(ncol(ovar2)/5))]

matrix_counts=as.matrix(GetAssayData(ovar3,assay = "RNA" , layer = "counts"))
###创建新的分组
ovar2$cnv_type=paste0(ovar2$tissue_type,"_",ovar2$seurat_clusters)
unique(ovar2$cnv_type)
####根据cluster分析
##write.table(ovar2_part$cnv_type,"cnv_type.txt",sep = "\t",quote = F , col.names = F)
####根据细胞注释分析
write.table(ovar2_part$celltype.main,"cnv_anno_type.txt",sep = "\t",quote = F , col.names = F)

gencode=read_tsv("hg38_gencode_v27.txt",col_names = c("gene","chr","start","end"))
gencode=gencode[!duplicated(gencode$gene),]
###提取基因交集
common_genes=intersect(gencode$gene,rownames(matrix_counts))
sort(unique(ovar2_part$celltype.main))####$后改

###创建Infercnv的object####cluster
#infercnv_obj=CreateInfercnvObject(raw_counts_matrix = matrix_counts[common_genes,],
#                                  annotations_file = "./cnv_type.txt",
#                                  delim = "\t",
#                                  gene_order_file = "./hg38_gencode_v27.txt",
#                                  ref_group_names = c("Normal_0", "Normal_1", "Normal_10" ,"Normal_11" ,"Normal_12", "Normal_13" ,"Normal_14",
#                                                      "Normal_15", "Normal_16", "Normal_17", "Normal_19", "Normal_2",  "Normal_20",
#                                                      "Normal_21", "Normal_22", "Normal_23", "Normal_24" ,"Normal_25" ,"Normal_26",
#                                                      "Normal_27", "Normal_28", "Normal_29", "Normal_3",  "Normal_30", "Normal_31",
#                                                      "Normal_32", "Normal_4",  "Normal_5",  "Normal_6",  "Normal_7",  "Normal_8", 
#                                                      "Normal_9"))


#infercnv_obj=infercnv::run(infercnv_obj,
#                           cutoff = 0.1,
#                           out_dir = "infercnv/",
#                           no_prelim_plot = T,
#                           cluster_by_groups = T,
#                           denoise = T,
#                           HMM = F,
#                           min_cells_per_gene = 10,
#                           num_threads = 10,
#                           write_expr_matrix = T)


save(infercnv_obj,file="infercnv_obj.Rdata    ,HMM = F")
infercnv::plot_cnv(infercnv_obj,
                   output_filename = "inferCNV_cluster_heatmap",
                   output_format = "pdf",
                   custom_color_pal = color.palette(c("#2067AE","white","#B11F2B")))

####cnvscore111111
setwd("F:/毕设/data3/finial/cnv")
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
grp$cluster=factor(grp$cluster,levels=c("B_Plasma","CD4+T_s1","CD4+T_s10","CD4+T_s11","CD4+T_s12","CD4+T_s13",    
                                        "CD4+T_s14","CD4+T_s15","CD4+T_s16","CD4+T_s17","CD4+T_s18","CD4+T_s19",    
                                        "CD4+T_s2","CD4+T_s20","CD4+T_s21","CD4+T_s22","CD4+T_s23","CD4+T_s24",   
                                        "CD4+T_s25","CD4+T_s26" ,    "CD4+T_s27" ,    "CD4+T_s28" ,    "CD4+T_s29",     "CD4+T_s3",  
                                        "CD4+T_s30","CD4+T_s31"  ,   "CD4+T_s32"  ,   "CD4+T_s33"  ,   "CD4+T_s34" ,    "CD4+T_s35", 
                                        "CD4+T_s36","CD4+T_s37"   ,  "CD4+T_s38"   ,  "CD4+T_s39"   ,  "CD4+T_s4"   ,   "CD4+T_s40" ,   
                                        "CD4+T_s41","CD4+T_s42"    , "CD4+T_s43"    , "CD4+T_s44" ,    "CD4+T_s45"   ,  "CD4+T_s5"   ,  
                                        "CD4+T_s6","CD4+T_s7",      "CD4+T_s8",      "CD4+T_s9"    ,  "CD8+T_s1"      ,"CD8+T_s10"    ,
                                        "CD8+T_s11","CD8+T_s12",     "CD8+T_s13",     "CD8+T_s14"   ,  "CD8+T_s15",     "CD8+T_s16"    ,
                                        "CD8+T_s17","CD8+T_s18" ,    "CD8+T_s19" ,    "CD8+T_s2"     , "CD8+T_s20" ,    "CD8+T_s21"    ,
                                        "CD8+T_s22","CD8+T_s23"  ,   "CD8+T_s24"  ,   "CD8+T_s25" ,    "CD8+T_s26"  ,   "CD8+T_s27"    ,
                                        "CD8+T_s28","CD8+T_s29"   ,  "CD8+T_s3"    ,  "CD8+T_s30"  ,   "CD8+T_s31"   ,  "CD8+T_s32"    ,
                                        "CD8+T_s33","CD8+T_s34"    , "CD8+T_s35"    , "CD8+T_s36"   ,  "CD8+T_s37"    , "CD8+T_s38"    ,
                                        "CD8+T_s39","CD8+T_s4"      ,"CD8+T_s40"     ,"CD8+T_s41"    , "CD8+T_s42"     ,"CD8+T_s43"    ,
                                        "CD8+T_s44","CD8+T_s45",     "CD8+T_s46",     "CD8+T_s47",     "CD8+T_s48",     "CD8+T_s49"    ,
                                        "CD8+T_s5","CD8+T_s50"  ,   "CD8+T_s51"  ,   "CD8+T_s52"  ,   "CD8+T_s53"  ,   "CD8+T_s54"    ,
                                        "CD8+T_s55","CD8+T_s56"   ,  "CD8+T_s57"   ,  "CD8+T_s58"   ,  "CD8+T_s59"   ,  "CD8+T_s6"     ,
                                        "CD8+T_s60" ,    "CD8+T_s61",     "CD8+T_s62",     "CD8+T_s63",     "CD8+T_s64",     "CD8+T_s65",    
                                        "CD8+T_s66"  ,   "CD8+T_s67" ,    "CD8+T_s68" ,    "CD8+T_s69" ,    "CD8+T_s7"  ,    "CD8+T_s70" ,   
                                        "CD8+T_s8"    ,  "CD8+T_s9"   ,   "Epithelia_s1",  "Epithelia_s10", "Epithelia_s11", "Epithelia_s12",
                                        "Epithelia_s13", "Epithelia_s14", "Epithelia_s15", "Epithelia_s16" ,"Epithelia_s17", "Epithelia_s18",
                                        "Epithelia_s19" ,"Epithelia_s2"  ,"Epithelia_s20" ,"Epithelia_s3",  "Epithelia_s4",  "Epithelia_s5" ,
                                        "Epithelia_s6",  "Epithelia_s7",  "Epithelia_s8",  "Epithelia_s9" , "Myeloid_s1"   , "Myeloid_s10"  ,
                                        "Myeloid_s11"  , "Myeloid_s12"  , "Myeloid_s13"  , "Myeloid_s14"   ,"Myeloid_s15"   ,"Myeloid_s16"  ,
                                        "Myeloid_s17" ,  "Myeloid_s18"  , "Myeloid_s19"   ,"Myeloid_s2",    "Myeloid_s20",   "Myeloid_s21"  ,
                                        "Myeloid_s22"  , "Myeloid_s23"  , "Myeloid_s24",   "Myeloid_s25",   "Myeloid_s26" ,  "Myeloid_s27"  ,
                                        "Myeloid_s28"   ,"Myeloid_s29"   ,"Myeloid_s3"  ,  "Myeloid_s30" ,  "Myeloid_s31"  , "Myeloid_s32"  ,
                                        "Myeloid_s33",   "Myeloid_s34",   "Myeloid_s35"  , "Myeloid_s4"   , "Myeloid_s5"    ,"Myeloid_s6"   ,
                                        "Myeloid_s7"  ,  "Myeloid_s8"  ,  "Myeloid_s9"    ,"NK_s1"         ,"NK_s10"    ,    "NK_s11" ,      
                                        "NK_s12"       , "NK_s13"       , "NK_s14"  ,      "NK_s15",        "NK_s16"     ,   "NK_s17"  ,     
                                        "NK_s18"        ,"NK_s19"        ,"NK_s2"    ,     "NK_s20" ,       "NK_s21"      ,  "NK_s22"   ,    
                                        "NK_s23"      ,  "NK_s24"       , "NK_s25"    ,    "NK_s26"  ,      "NK_s27"    ,    "NK_s28"    ,   
                                        "NK_s29"       , "NK_s3"         ,"NK_s30"     ,   "NK_s31"   ,     "NK_s32"     ,   "NK_s33"     ,  
                                        "NK_s4"         ,"NK_s5"         ,"NK_s6"       ,  "NK_s7"     ,    "NK_s8"       ,  "NK_s9"       , 
                                        "Tregs_s1"  ,    "Tregs_s10"  ,   "Tregs_s11"    , "Tregs_s12"  ,   "Tregs_s13"    , "Tregs_s14"    ,
                                        "Tregs_s15"  ,   "Tregs_s16"   ,  "Tregs_s17"     ,"Tregs_s18"   ,  "Tregs_s19"     ,"Tregs_s2"     ,
                                        "Tregs_s20"   ,  "Tregs_s21"    , "Tregs_s22" ,    "Tregs_s23"    , "Tregs_s24" ,    "Tregs_s25"    ,
                                        "Tregs_s26"    , "Tregs_s27"     ,"Tregs_s28"  ,   "Tregs_s29"     ,"Tregs_s3"   ,   "Tregs_s30"    ,
                                        "Tregs_s31"     ,"Tregs_s32"  ,   "Tregs_s33"   ,  "Tregs_s34"    , "Tregs_s4"    ,  "Tregs_s5"     ,
                                        "Tregs_s6",      "Tregs_s7"    ,  "Tregs_s8"     , "Tregs_s9" ))
row.names(grp)=row.names(cnv_score)
grp=grp[-1]
identical(row.names(cnv_score),row.names(grp))
cnv_score$cluster=grp$cluster
colnames(cnv_score)=c("Score","Cluster")
ggboxplot(cnv_score,"Cluster","Score",fill="Cluster")+
  coord_cartesian(ylim = c(1,1000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(cnv_score$Cluster)

##小提琴图

ggplot(cnv_score,aes(x=Cluster,y=Score,fill=Cluster))+
  geom_violin()+
  coord_cartesian(ylim = c(0,2000))+
  theme_minimal()+
  labs(x="Cluster",y="CNV Score")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

####cnvscore222
anno=read.table("cnv_anno_type.txt",row.names = 1)
data=read.table("infercnv/infercnv.observations.txt",header=T)
expr=data%>%as.matrix()
expr.scale=scale(t(expr))
tmp1=sweep(expr.scale,2,apply(expr.scale,2,min),'-')
tmp2=apply(expr.scale,2,max)-apply(expr.scale,2,min)
expr_1=t(2*sweep(tmp1,2,tmp2,"/")-1)
cnv_score=as.data.frame(colSums(expr_1*expr_1))
colnames(cnv_score)="cnv_score"
cnv_score=rownames_to_column(cnv_score,var='cell')
anno$cell=rownames(anno)
colnames(anno)[1]="cluster"

cnv_score$cell=gsub("\\.","-",cnv_score$cell)
head(cnv_score)
test=merge(cnv_score,anno,by="cell",all=F)
ggplot2::ggplot(test,aes(x=cluster,y=cnv_score))+
  geom_violin(aes(fill = cluster),ces=1.2)+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),legend.position = 'none')






#########################annotation 
infercnv_obj=CreateInfercnvObject(raw_counts_matrix = matrix_counts[common_genes,],
                                  annotations_file = "./cnv_anno_type.txt",
                                  delim = "\t",
                                  gene_order_file = "./hg38_gencode_v27.txt",
                                  ref_group_names = c("T cell","Monocytic","B_Plasma cell"))

infercnv_obj=infercnv::run(infercnv_obj,
                           cutoff = 0.1,
                           out_dir = "infercnv_annotation/",
                           no_prelim_plot = T,
                           cluster_by_groups = T,
                           denoise = T,
                           min_cells_per_gene = 10,
                           num_threads = 10,
                           write_expr_matrix = T)
infercnv::plot_cnv(infercnv_obj,
                   output_filename = "inferCNV_cluster_heatmap",
                   output_format = "pdf",
                   custom_color_pal = color.palette(c("#2067AE","white","#B11F2B")))



#####################################拷贝数变异2#########



#########################11.29
###cnv----以成纤维细胞和内皮细胞作为参考
#########################################################拷贝数变异
BiocManager::install("infercnv")
DimPlot(ovar2,reduction = "umap",label=T,split.by = "tissue_type")
DimPlot(ovar2,reduction = "tsne",label=T,group.by = "orig.ident")

table(ovar3$celltype.main.22)
table(ovar3$tissue_type)
Idents(ovar3)="celltype.main.22"
Idents(ovar3)###抽样
ovar3_part=ovar3[,sample(1:ncol(ovar3),round(ncol(ovar3)/5))]

#matrix_counts=as.matrix(GetAssayData(ovar2_part,assay = "RNA" , layer = "counts"))
matrix_counts=as.matrix(GetAssayData(ovar3,assay = "RNA" , layer = "counts"))
###创建新的分组
ovar3$cnv_type=paste0(ovar3$tissue_type,"_",ovar3$celltype.main.22)
unique(ovar3$celltype.main.22)
###根据cluster分析
write.table(ovar2$cnv_type,"cnv_type.txt",sep = "\t",quote = F , col.names = F)
###根据细胞注释分析
write.table(ovar2_part$celltype.main,"cnv_anno_type.txt",sep = "\t",quote = F , col.names = F)

write.table(ovar3$celltype.main.22,"cnv_anno_type.txt",sep = "\t",quote = F , col.names = F)

gencode=read_tsv("F:/毕设/data3/hg38_gencode_v27.txt",col_names = c("gene","chr","start","end"))
gencode=gencode[!duplicated(gencode$gene),]
###提取基因交集
common_genes=intersect(gencode$gene,rownames(matrix_counts))
sort(unique(ovar3$celltype.main.22))####$后改

###创建Infercnv的object####cluster
options(scipen = 100)
infercnv_obj=CreateInfercnvObject(raw_counts_matrix = matrix_counts[common_genes,],
                                  annotations_file = "./cnv_anno_type.txt",
                                  delim = "\t",
                                  gene_order_file = "./hg38_gencode_v27.txt",
                                  ref_group_names = c("Fibroblast","Endothelia"))



infercnv_obj=infercnv::run(infercnv_obj,
                           cutoff = 0.1,
                           out_dir = "infercnv/",
                           no_prelim_plot = T,
                           cluster_by_groups = T,
                           denoise = T,
                           min_cells_per_gene = 10,
                           num_threads = 10,
                           write_expr_matrix = T)

     save(infercnv_obj,file="infercnv_obj.Rdata",HMM = F)
infercnv::plot_cnv(infercnv_obj,
                   output_filename = "inferCNV_cluster_heatmap",
                   output_format = "png",
                   custom_color_pal = color.palette(c("#2067AE","white","#B11F2B")))

grp=read.table("F:/毕设/data3/finial/F3----上皮细胞/xincvv/infercnv/infercnv.observation_groupings.txt",sep = "",header = T,row.names = NULL)
obs=data.table::fread("F:/毕设/data3/finial/F3----上皮细胞/xincvv/infercnv/infercnv.observations.txt",data.table = F)%>%column_to_rownames(var = 'V1')
##F:\毕设\data3\finial\F3----上皮细胞\xincvv\infercnv

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
grp$cluster <- sub("_[^_]+$", "", grp$Dendrogram.Group)
#grp$cluster=sapply(strsplit(grp$Dendrogram.Group,"_"),function(x) paste(x[1],x[2],sep = "_"))
sort(unique(grp$cluster))
grp$cluster=factor(grp$cluster,levels=c("B_Plasma" ,"CD4+T" ,"CD8+T" ,"Epithelia" ,"Myeloid" ,"NK" ,"Tregs"))
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







