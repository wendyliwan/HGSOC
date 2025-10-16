######对每个亚型细胞进行Norma-Cancer差异分析------pseudobulks#############################3
#############基质细胞#########
bs=split(colnames(Stro_sce),Stro_sce$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(Stro_sce)%in%bs[[x]]
    data=GetAssayData(Stro_sce,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(Stro_sce@meta.data[,c('orig.ident',"tissue_type")])
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
#write.csv(DEG_DESeq2,"./差异分析/Stromal_DEG_DESeq2.csv")

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
p
ggsave("./差异分析/Stro_volcano_unlabel.png",p,width = 8,height = 7)

top_20 <-bind_rows(
  DEG_DESeq2 %>%
    filter(Type == 'up') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(13),
  DEG_DESeq2 %>%
    filter(Type == 'down') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(13)
)
#write.csv(top_20,"./差异分析/Stro_DEG_top10.csv")
rownames(top_20)=top_20$Gene_Sambol

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
volc_plot2 <- p +
  geom_label_repel(
    data = top_20,
    aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
    size = 3.5,  # 从 2 改为 3.5（建议值，可根据需求调整）
    segment.size = 0.4,  # 可选：同步调整连接线粗细
    max.overlaps = 30    # 可选：避免密集标签时的重叠
  )
volc_plot2
ggsave("./差异分析/Stro_volcano.png",volc_plot2,width = 8,height = 7)

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
###选取感兴趣的基因集
select_ego1_res=ego1_res%>%
  dplyr::filter(grepl("pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))

set.seed(123)

# 随机选择10个BP、MF、CC的Description
selected_bp <- sample(ego1_res$Description[ego1_res$ONTOLOGY == "BP"], 5)
selected_mf <- sample(ego1_res$Description[ego1_res$ONTOLOGY == "MF"], 5)
selected_cc <- sample(ego1_res$Description[ego1_res$ONTOLOGY == "CC"], 5)

selected_bp <- ego1_res$ONTOLOGY=="BP"
selected_mf <- ego1_res$ONTOLOGY == "MF"
selected_cc <- ego1_res$ONTOLOGY == "CC"


selected_bp <- subset(ego1_res, ONTOLOGY == "BP")
selected_bp=selected_bp[order(selected_bp$LogP, decreasing = TRUE), ]

selected_cc <- subset(ego1_res, ONTOLOGY == "CC")
selected_cc=selected_cc[order(selected_cc$LogP, decreasing = TRUE), ]

selected_mf <- subset(ego1_res, ONTOLOGY == "MF")
selected_mf=selected_mf[order(selected_mf$LogP, decreasing = TRUE), ]



# 合并选择的Description
selected_descriptions <- c(selected_bp[1:5,], selected_mf[1:5,], selected_cc[1:5,])
selected_descriptions <- c(
  selected_bp$Description[1:5],  # 提取前5个BP的Description
  selected_mf$Description[1:5],  # 提取前5个MF的Description
  selected_cc$Description[1:5]   # 提取前5个CC的Description
)
# 筛选数据
filtered_go_results <- ego1_res[ego1_res$Description %in% selected_descriptions, ]
ggplot(selected_descriptions, aes(x = reorder(Description, RichFactor), y = RichFactor, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +  # 按ONTOLOGY分面，纵向排列
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "GO Enrichment Analysis", x = "Description", y = "Rich Factor") +
  scale_fill_manual(values = c("BP" = "blue", "MF" = "red", "CC" = "green")) +
  coord_flip()  # 将条形图改为横向排列


# 绘制条形图，按logp设置渐变颜色
# ggplot(filtered_go_results, aes(x = reorder(Description,Count), y = Count, fill = LogP)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1)+  # 按ONTOLOGY分面，纵向排列
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "GO Enrichment Analysis", 
#         x = "Description", 
#         y = "Count", 
#         fill = "LogP") +  # 图例标题
#   # scale_fill_gradient(low = "grey", high = "purple") +  # 设置渐变颜色
#   coord_flip()  # 将条形图改为横向排列
#   # 
# 
# ggplot(ego1_res[1:40,],aes(Count,Description))+
#   geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
#   #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
#   ggtitle("GO-BP")+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))

ontology_colors <- list(
  BP = c("lightblue", "darkblue"),
  CC = c("pink", "red"),
  MF = c("lightgreen", "darkgreen") # 分子功能（MF）
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
  coord_flip()
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
  head(30)  # 取前20个

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

#################
# set.seed(123)
# selected_cp <- sample(kk1_res$Description[kk1_res$category == "Cellular Processes"], 5)
# selected_gp <- sample(kk1_res$Description[kk1_res$category == "Genetic Information Processing"], 5)
# selected_me <- sample(kk1_res$Description[kk1_res$category == "Metabolism"], 5)
# selected_ep <- sample(kk1_res$Description[kk1_res$category == "Environmental Information Processing"], 5)
# selected_hd <- sample(kk1_res$Description[kk1_res$category == "Human Diseases"], 5)
# selected_os <- sample(kk1_res$Description[kk1_res$category == "Organismal Systems"], 5)
# 
# selected_descriptions <- c(selected_cp, selected_gp, selected_me,selected_ep,selected_hd,selected_os)
# 
# filtered_go_results <- kk1_res[kk1_res$Description %in% selected_descriptions, ]
# ggplot(selected_descriptions, aes(x = reorder(Description, RichFactor), y = RichFactor, fill = category)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~ category, scales = "free_y", ncol = 1) +  # 按ONTOLOGY分面，纵向排列
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "GO Enrichment Analysis", x = "Description", y = "Rich Factor") +
#   scale_fill_manual(values = c("Cellular Processes" = "#FF0000FF", "Genetic Information Processing" = "#FF9900FF", 
#                                "Metabolism" = "#FFCC00FF" ,"Environmental Information Processing"="#00FF00FF",
#                                "Human Diseases"="#6699FFFF","Organismal Systems" ="#CC33FFFF")) +
#   coord_flip()  # 将条形图改为横向排列
#################
ggplot(kk1_res[1:30,],aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  ggtitle("KEGG Enrichment Analysis")+
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









#####内皮细胞##############################################################################
cell.ids <- which(Stro_sce@meta.data$subcelltype=="Endothelial")
Endo_sce=subset(Stro_sce,cells = cell.ids)
bs=split(colnames(Endo_sce),Endo_sce$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(Endo_sce)%in%bs[[x]]
    data=GetAssayData(Endo_sce,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(Endo_sce@meta.data[,c('orig.ident',"tissue_type")])
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
write.csv(DEG_DESeq2,"Endothelial_DEG_DESeq2.csv")
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
ggsave("./差异分析/Endothelial/Endo_volcano_unlabel.png",p,width = 8,height = 7)

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
write.csv(top_20,"./差异分析/Endothelial/Endo_DEG_top10.csv")
rownames(top_20)=top_20$Gene_Sambol

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("./差异分析/Endothelial/Endo_volcano.png",volc_plot2,width = 8,height = 7)

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
###选取感兴趣的基因集
select_ego1_res=ego1_res%>%
  dplyr::filter(grepl("pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))

ggplot(select_ego1_res,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  ggtitle("GO")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
#-----------------------------kegg
kk1=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk1_res=kk1@result
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










#############成纤维----平滑肌细胞#################################################################
Idents(Stro_sce)
cell.ids <- which(Stro_sce@meta.data$subcelltype=="myoFibroblast")
myFib_sce=subset(Stro_sce,cells = cell.ids)
bs=split(colnames(myFib_sce),myFib_sce$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(myFib_sce)%in%bs[[x]]
    data=GetAssayData(myFib_sce,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(myFib_sce@meta.data[,c('orig.ident',"tissue_type")])
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
write.csv(DEG_DESeq2,"./差异分析/myoFibroblast/myFib_DEG_DESeq2.csv")
table(DEG_DESeq2$Type)

p=ggplot(DEG_DESeq2,aes(log2FoldChange,-log10(pvalue)))+
  geom_point(size=2.5,
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
ggsave("./差异分析/myoFibroblast/myFib_volcano_unlabel.png",p,width = 8,height = 7)

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
rownames(top_20)=top_20$Gene_Sambol
write.csv(top_20,"./差异分析/myoFibroblast/myFib_DEG_top10.csv")


volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("./差异分析/myoFibroblast/myFib_volcano.png",volc_plot2,width = 8,height = 7)

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
###选取感兴趣的基因集
select_ego1_res=ego1_res%>%
  dplyr::filter(grepl("pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))

ggplot(select_ego1_res,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  ggtitle("GO")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
#-----------------------------kegg
kk2=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk2_res=kk2@result
kk2_res$Group="Cancer vs Normal"
kk2_res$LogP=-log(kk2_res$p.adjust)

ggplot(kk2_res[1:40,],aes(Count,Description))+
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








#############成纤维细胞##########################################################################
Idents(Stro_sce)
cell.ids <- which(Stro_sce@meta.data$subcelltype=="Fibroblast")
Fib_sce=subset(Stro_sce,cells = cell.ids)
bs=split(colnames(Fib_sce),Fib_sce$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(Fib_sce)%in%bs[[x]]
    data=GetAssayData(Fib_sce,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(Fib_sce@meta.data[,c('orig.ident',"tissue_type")])
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
write.csv(DEG_DESeq2,"./差异分析/Fibroblast/Fibroblast_DEG_DESeq2.csv")
table(DEG_DESeq2$Type)

p=ggplot(DEG_DESeq2,aes(log2FoldChange,-log10(pvalue)))+
  geom_point(size=2.5,
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
ggsave("./差异分析/Fibroblast/Fib_volcano_unlabel.png",p,width = 8,height = 7)

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
rownames(top_20)=top_20$Gene_Sambol
write.csv(top_20,"./差异分析/Fibroblast/Fib_DEG_top10.csv")


volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("./差异分析/Fibroblast/Fib_volcano.png",volc_plot2,width = 8,height = 7)


##---------------------------------------Go
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

ggplot(select_ego1_res,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  ggtitle("GO")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
#-----------------------------kegg
kk3=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk3_res=kk3@result
kk3_res$Group="Cancer vs Normal"
kk3_res$LogP=-log(kk3_res$p.adjust)

ggplot(kk3_res[1:40,],aes(Count,Description))+
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










#############间皮细胞#######################################################################
Idents(Stro_sce)
cell.ids <- which(Stro_sce@meta.data$subcelltype=="Mesothelial")
Meso_sce=subset(Stro_sce,cells = cell.ids)
bs=split(colnames(Meso_sce),Meso_sce$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(Meso_sce)%in%bs[[x]]
    data=GetAssayData(Meso_sce,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
colnames(exprset)=names(bs)
phe=unique(Meso_sce@meta.data[,c('orig.ident',"tissue_type")])
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
write.csv(DEG_DESeq2,"./差异分析/Mesothelial/Meso_DEG_DESeq2.csv")
table(DEG_DESeq2$Type)

p=ggplot(DEG_DESeq2,aes(log2FoldChange,-log10(pvalue)))+
  geom_point(size=2.5,
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
ggsave("./差异分析/Mesothelial/Meso_volcano_unlabel.png",p,width = 8,height = 7)

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
rownames(top_20)=top_20$Gene_Sambol
write.csv(top_20,"./差异分析/Mesothelial/Meso_DEG_top10.csv")

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("./差异分析/Mesothelial/Meso_volcano.png",volc_plot2,width = 8,height = 7)

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
###选取感兴趣的基因集
select_ego1_res=ego1_res%>%
  dplyr::filter(grepl("pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))

ggplot(ego1_res,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  #scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  ggtitle("GO")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
#-----------------------------kegg
kk4=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk4_res=kk4@result
kk4_res$Group="Cancer vs Normal"
kk4_res$LogP=-log(kk4_res$p.adjust)

ggplot(kk4_res[1:40,],aes(Count,Description))+
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




