install.packages("XML")
install.packages("methods")
library(XML)
library(rjson)
library(methods)
library(tidyverse)
setwd("F:/毕设/TCGA/")
json=jsonlite::fromJSON("metadata.cart.2024-12-04.json")
sample_id=sapply(json$associated_entities,function(x){x[,1]})
file_sample=data.frame(sample_id,file_name=json$file_name)
count_file=list.files('F:/毕设/TCGA/clinical/gdc_download_20241204_090135.560622/',
                      pattern = '*.xml',recursive = TRUE)
count_file_name=strsplit(count_file,split='/')
count_file_name=sapply(count_file_name,function(x){x[2]})
matrix=data.frame(matrix(nrow = 60660,ncol = 0))

for(i in 1:length(count_file)){
  path=paste0('F:/毕设/TCGA/clinical/gdc_download_20241204_090135.560622/',count_file[i])
  data=read.delim(path,fill=TRUE,header=FALSE,row.names = 1)
  colnames(data)=data[2,]
  data=data[-c(1:6),]
  data=data[6]
  colnames(data)=file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix=cbind(matrix,data)
}
matrix <- data.frame()  # 初始化

for(i in 1:length(count_file)) {
  path <- paste0('F:/毕设/TCGA/clinical/gdc_download_20241204_090135.560622/', count_file[i])
  
  # 读取数据并处理重复行名
  data <- read.delim(path, fill = TRUE, header = FALSE, row.names = NULL)
  rownames(data) <- make.unique(as.character(data[, 1]))  # 确保行名唯一
  data <- data[, -1]  # 移除原第一列
  
  # 设置列名并筛选数据
  colnames(data) <- data[2, ]
  data <- data[-c(1:6), ]
  data <- data[6]  # 保持为数据框
  
  # 匹配样本ID
  sample_id <- file_sample$sample_id[which(file_sample$file_name == count_file_name[i])]
  colnames(data) <- sample_id
  
  # 合并到矩阵
  matrix <- cbind(matrix, data)
}







dir="F:/毕设/TCGA/clinical/gdc_download_20241204_090135.560622"
all_files=list.files(path=dir,pattern = '*.xml$',recursive = T)
all_files[1:10]
cl=lapply(all_files,function(x){
  result=xmlParse(file = file.path(dir,x))
  rootnode=xmlRoot(result)
  xmldataframe=xmlToDataFrame(rootnode[2])
  return(t(xmldataframe))
})
clinical=t(do.call(cbind,cl))
write.table(clinical,file="clinical.txt",sep = "\t",quote = F,row.names = F)


setwd("F:/毕设/TCGA/rna/")
count_files=dir("gdc_download_20250403_100107.570226/",pattern="*isoforms.quantification.txt$",recursive=T)
exp=list()
for(i in 1:length(count_files)){
  exp[[i]]=read.table(paste0("gdc_download_20250403_100107.570226/",count_files[[i]]),sep = "\t",header = T)%>%
    dplyr::select(c(6,4))%>%
    group_by(miRNA_region)%>%
    summarise(reads_per_million_miRNA_mapped=sum(reads_per_million_miRNA_mapped))
}
m=Reduce(function(x,y) merge(x,y,by="miRNA_region",all=T),exp)
m[is.na(m)]=0
exp=column_to_rownames(m,var="miRNA_region")
exp=exp[-((nrow(exp)-2):nrow(exp)),]
rownames(exp)=str_remove(rownames(exp),"mature,")
library(miRBaseVersions.db)
mh=select(miRBaseVersions.db,
          keys=rownames(exp),
          keytype = "MIMAT",
          columns = c("ACCESSION","NAME","VERSION"))
mh=mh[mh$VERSION=="21",]
mh=mh[match(rownames(exp),mh$ACCESSION),]
rownames(exp)=mh$NAME
meta=jsonlite::fromJSON("metadata.cart.2025-04-03 (1).json")
ID=sapply(meta$associated_entities,
          function(x){x$entity_submitter_id})
file2id=data.frame(file_name=meta$file_id,ID=ID)

count_files2=stringr::str_split(count_files,"/",simplify = T)[,1]
file2id=file2id[match(count_files2,file2id$file_name),]
colnames(exp)=file2id$ID
k=apply(exp,1,function(x) sum(x>0)>100);table(k)
exp=exp[k,]
exp=as.matrix(exp)
exp1=data.frame(ID=rownames(exp),exp)
colnames(exp1)=gsub('[.]','-',colnames(exp1))
write.table(exp1,'miRNA.RPM.txt',sep = "\t",quote=F,row.names = F)



setwd("F:/毕设/TCGA/GE/")
json=jsonlite::fromJSON("metadata.cart.2025-04-04.json")
sample_id=sapply(json$associated_entities,function(x){x[,1]})
file_sample=data.frame(sample_id,file_name=json$file_name)
count_file=list.files('gdc_download_20250404_005508.897618/',
                      pattern = '*.tsv',recursive = TRUE)
count_file_name=strsplit(count_file,split='/')
count_file_name=sapply(count_file_name,function(x){x[2]})
matrix=data.frame(matrix(nrow = 60660,ncol = 0))

for(i in 1:length(count_file)){
  path=paste0('gdc_download_20250404_005508.897618/',count_file[i])
  data=read.delim(path,fill=TRUE,header=FALSE,row.names = 1)
  colnames(data)=data[2,]
  data=data[-c(1:6),]
  data=data[6]  ##TPM
  #=data[3] ##count
  colnames(data)=file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix=cbind(matrix,data)
}

path=paste0('gdc_download_20250404_005508.897618/',count_file[1])
data=as.matrix(read.delim(path,fill=TRUE,header=FALSE,row.names = 1))
gene_name=data[-c(1:6),1]
matrix0=cbind(gene_name,matrix)
gene_type=data[-c(1:6),2]
matrix0=cbind(gene_type,matrix0)
matrix0=aggregate(.~gene_name,data=matrix0,max)
matrix0=subset(x=matrix0,gene_type=="protein_coding")
rownames(matrix0)=matrix0[,1]
matrix0=matrix0[,-c(1,2)]
matrix1=data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1)=gsub('[.]','-',colnames(matrix1))
write.table(matrix1,'TCGA_OV_TPM.txt',sep = "\t",quote=F,row.names = F)


#############生存分析######
setwd("F:/毕设/TCGA/")
cli=read.csv("time_ov.csv",header=T,sep=",",check.names = F,row.names = 1)
cli$time=cli$time/365
data=read.table("TCGA_OV_TPM.txt",header=T,sep="\t",check.names = F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)
data=t(data)
rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]','-',rownames(data))
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli,data)


group=ifelse(rt[,"HAL-E"]> quantile(rt[,"HAL-E"]) [2],"High","Low")
rt[,"group"]=group
length=length(levels(factor(group)))

# diff=survdiff(Surv(time,state) ~group,data=rt)
# pValue=1-pchisq(diff$chisq,df=length-1)
# pValue=paste0("p=",sprintf("%.03f",pValue))

sur.cut <- surv_cutpoint(rt,   time= 'time',
                         event = 'state' ,
                         variables = 'HAL-E' )
sur.cut
sur.cat <- surv_categorize(sur.cut)
###这一步的主要原理是，放弃以前所用的中位值来定义高低组的方法，
# 采用不同的阈值来重新定义高低分组以达到最低的P值。我们看到经过这
# 样的步骤，高组别的患者明显多于低组别。
table(sur.cat$ITGAV)

fit <- survfit(Surv(time , state) ~HAL-E, data =sur.cat)

 # fit=survfit(Surv(time,state) ~group,data=rt)  

bioCol=c("Firebrick3","MediumSeaGreen","#6e568c","#223d6c")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(
  fit,   data = sur.cat,   
  pval = TRUE,           
  conf.int = TRUE, 
  # legend.title = "CUL3",
  xlim = c(0,5),  
  xlab = "Time in years",  
  break.time.by = 1,  
  pval.size = 8,
  risk.table = "absolute", 
  risk.table.y.text.col = T,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 5,
  # legend.labs =  c("high", "Low"),   
  palette = c("#E41A1C","#377EB8"),surv.median.line = "hv",
  font.y  = c(16, "bold"), 
  font.x  = c(16, "bold"),
  legend = "top",
  font.legend = c(16, "bold"))
# surPlot=ggsurvplot(fit,
#                    data=rt,
#                    conf.int=F,
#                    pval=TRUE,
#                    pval.size=6,
#                    legend.title="ITGAV experssion",
#                    legend.labs=levels(factor(rt[,"group"])),
#                    legend=c(0.8,0.8),
#                    font.legend=10,
#                    xlab="Time(years)",
#                    break.time.by=2,
#                    palette=bioCol,
#                    surv.median.line="hv",
#                    risk.table=T,
#                    cumevents=F,
#                    risk.table.height=.25
#                    )



surPlot

# 加载必要的包
library(survival)
library(survminer)
library(tidyverse)

# 1. 数据准备
cli <- read.csv("time_ov.csv", header=T, sep=",", check.names=F, row.names=1)
cli$time <- cli$time/365  # 将时间转换为年

data <- read.table("TCGA_OV_TPM.txt", header=T, sep="\t", check.names=F, row.names=1)
data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), 
               dimnames=list(rownames(data), colnames(data)))
data <- t(data)
rownames(data) <- substr(rownames(data), 1, 12)
rownames(data) <- gsub('[.]', '-', rownames(data))

# 2. 合并临床数据和表达数据
sameSample <- intersect(row.names(data), row.names(cli))
data <- data[sameSample, ]
cli <- cli[sameSample, ]
rt <- cbind(cli, data)

# 3. 定义要分析的基因列表（这里使用所有基因，也可以自定义）
genes <- colnames(rt)[3:ncol(rt)]  # 假设前两列是time和state
genes <- genes[grep("\\.", genes, invert = TRUE)]  # invert=TRUE表示反向匹配
# 4. 创建目录保存结果
if(!dir.exists("Survival_Plots")) dir.create("Survival_Plots")

# 5. 批量生存分析函数
batch_survival <- function(rt, genes) {
  results <- data.frame(Gene=character(), pValue=numeric(), stringsAsFactors=F)
  
  for(gene in genes) {
    tryCatch({
      # 5.1 分组（中位数）
      gene_exp <- rt[, gene]
      if(all(is.na(gene_exp))) next  # 跳过全NA的基因
      
      cutoff <- quantile(gene_exp, probs=0.5, na.rm=TRUE)
      rt$group <- ifelse(gene_exp > cutoff, "High", "Low")
      
      # 5.2 生存分析
      diff <- survdiff(Surv(time, state) ~ group, data=rt)
      pValue <- 1 - pchisq(diff$chisq, df=1)
      
      # 保存结果
      results <- rbind(results, data.frame(Gene=gene, pValue=pValue))
      
      # 5.3 仅对显著基因绘图
      if(pValue < 0.05) {
        fit <- survfit(Surv(time, state) ~ group, data=rt)
        
        # 绘图
        survPlot <- ggsurvplot(
          fit, data=rt,
          conf.int=FALSE,
          pval=paste0("p = ", format.pval(pValue, digits=3)),
          pval.size=6,
          legend.title=paste0(gene, " expression"),
          legend.labs=c("High", "Low"),
          legend=c(0.8, 0.8),
          font.legend=10,
          xlab="Time (years)",
          break.time.by=2,
          palette=c("#E7B800", "#2E9FDF"),
          surv.median.line="hv",
          risk.table=TRUE,
          cumevents=FALSE,
          risk.table.height=0.25
        )
        
        # 保存图片
        ggsave(
          filename=paste0("Survival_Plots/", gene, "_survival.pdf"),
          plot=print(survPlot),
          width=8, height=6
        )
      }
    }, error=function(e) {
      message(paste("Error in gene", gene, ":", e$message))
    })
  }
  
  # 6. 多重检验校正
  results$adj.p <- p.adjust(results$pValue, method="BH")
  results <- results[order(results$pValue), ]
  
  # 7. 保存结果表格
  write.csv(results, "Survival_Analysis_Results.csv", row.names=FALSE)
  
  return(results)
}

# 运行批量分析
survival_results <- batch_survival(rt, genes)

# 查看显著结果
significant_genes <- subset(survival_results, adj.p < 0.05)
print(paste("Found", nrow(significant_genes), "significant genes"))




  
  
  
  
  
  
  
  
  
  
  
  
  









#######################
DEG_DESeq2=read.csv("DEG_DESeq2.csv")
no_decimal_genes <- !grepl("\\.", DEG_DESeq2$Gene_Sambol)

# 过滤掉包含小数点的基因名
filtered_DEG <- DEG_DESeq2[no_decimal_genes, ]
table(DEG_DESeq2$Type)

p=ggplot(filtered_DEG,aes(log2FoldChange,-log10(pvalue)))+
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
ggsave("volcano_unlabel.png",p,width = 8,height = 7)

top_20 <-bind_rows(
  filtered_DEG %>%
    filter(Type == 'up') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(20),
  filtered_DEG %>%
    filter(Type == 'down') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(20)
)
rownames(top_20)=top_20$Gene_Sambol
write.csv(top_20,"filter_DEG_top10.csv")
top_20=na.omit(top_20)

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("volcano.png",volc_plot2,width = 8,height = 7)


####################################################
Idents(ovar3)=ovar3$tissue_type

bs=split(colnames(ovar3),ovar3$orig.ident)
exprset=do.call(
  cbind,lapply(names(bs),function(x){
    kp=colnames(ovar3)%in%bs[[x]]
    data=GetAssayData(ovar3,assay = "RNA",layer = "counts")###=exprmatrix
    data=as.matrix(data)
    rowSums(as.matrix(data[,kp]))
  })
)
no_decimal_genes <- !grepl("\\.", rownames(exprset))

# 过滤掉包含小数点的基因名
exprset <- exprset[no_decimal_genes, ]

# 找出不包含0的行的索引
no0 <- which(!apply(exprset, 1, function(x) any(x == 0)))

# 删除包含0的行
exprset.2 <- exprset[no0, , drop = FALSE]
write.csv(exprset,file="exprset.csv")
write.csv(exprset.2,file="exprsetno0.csv")
colnames(exprset)=names(bs)
phe=unique(ovar3@meta.data[,c('orig.ident',"tissue_type")])
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
table                                                                                                                                                                                                                                                                                                                                                                                    

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
ggsave("volcano_unlabel.png",p,width = 8,height = 7)

top_20 <-bind_rows(
  DEG_DESeq2 %>%
    filter(Type == 'up') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(25),
  DEG_DESeq2 %>%
    filter(Type == 'down') %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    head(60)
)
rownames(top_20)=top_20$Gene_Sambol
write.csv(top_20,"DEG_top10.csv")
top_20=na.omit(top_20)

volc_plot2 <- p +
  geom_label_repel(data = top_20,
                   aes(log2FoldChange, -log10(pvalue), label = rownames(top_20)),
                   size = 2)
ggsave("volcano.png",volc_plot2,width = 8,height = 7)

sig_DEG=subset(DEG_DESeq2,Type!='ns')
genelist=bitr(sig_DEG$Gene_Sambol,fromType = "SYMBOL",
              toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
ego.nvsc=enrichGO(gene=genelist$ENTREZID,
              OrgDb =org.Hs.eg.db,
              ont = "ALL",
              minGSSize = 1,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.2,
              readable = T)
ego.nvsc.res=ego.nvsc@result
ego.nvsc.res$Group="Cancer vs Normal"
ego.nvsc.res$LogP=-log(ego1_res$p.adjust)
###选取感兴趣的基因集
ego.nvsc.res.select=ego.nvsc.res%>%
  dplyr::filter(grepl("ovarian|pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))

top <-bind_rows(
  ego.nvsc.res %>%
    filter(ONTOLOGY == 'BP') %>%
    arrange(pvalue) %>%
    head(5),
  ego.nvsc.res %>%
    filter(ONTOLOGY == 'CC') %>%
    arrange(pvalue) %>%
    head(5),
  ego.nvsc.res %>%
    filter(ONTOLOGY == 'MF') %>%
    arrange(pvalue) %>%
    head(5)
)

top$term <- paste(top$ID, top$Description, sep = ': ') #将ID与Description合并成新的一列
top$term <- factor(top$term, levels = top$term,ordered = T)

p=ggplot(top,aes(y=term,x=zScore,fill=Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  ggtitle("GO")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))

p<- p+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')
Idents(ovar3)






kk4=enrichKEGG(gene = genelist$ENTREZID,
               keyType = 'kegg',
               organism = 'hsa',
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.1)
kk4=kk4@result
kk4$Group="Cancer vs Normal"
kk4$LogP=-log(kk4_res$p.adjust)

kk4=kk4%>%
  dplyr::filter(grepl("ovarian|pathway|signal|cancer",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))


convert_gene_ids_to_names <- function(gene_ids_str, OrgDb) {
  # 拆分基因ID字符串为单个ID
  gene_ids <- unlist(strsplit(gene_ids_str, "/"))
  
  # 使用bitr函数进行ID转换
  gene_info <- bitr(gene_ids, fromType="ENTREZID", toType="SYMBOL", OrgDb=OrgDb)
  
  # 提取转换后的基因名称（处理可能存在的NA值）
  gene_names <- ifelse(is.na(gene_info$SYMBOL), gene_ids, gene_info$SYMBOL)
  
  # 重新组合为以"/"分割的基因名称字符串
  gene_names_str <- paste(gene_names, collapse="/")
  
  return(gene_names_str)
}
kk4_res$gene <- apply(kk4_res, 1, function(row) convert_gene_ids_to_names(row["geneID"], OrgDb='org.Hs.eg.db'))

no_decimal_genes <- !grepl("\\.", rownames(exprset))

# 过滤掉包含小数点的基因名
filtered_exp <- exprset[no_decimal_genes, ]
kk4_res$LogP=-log(kk4_res$p.adjust)

kk4_res.f=kk4_res%>%
  dplyr::filter(grepl("pathway|immune|Ovarian",Description))%>%
  dplyr::arrange(dplyr::desc(LogP),dplyr::desc(Description))%>%
  mutate(Description=forcats::fct_inorder(Description))


ggplot(kk4_res.f,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP),stat = "identity")+
  ggtitle("KEGG")+
  facet_wrap(~category, ncol = 1, scales = "free_y", strip.position = "right")+
  scale_fill_gradient(low="#d5cabd",high = "#dc0000")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 0.5))
write.csv(kk4_res.f,file="kegg_pathway.csv")

kk4_res.f2=kk4_res.f%>%
  dplyr::filter(grepl("Organismal Systems|Environmental Information Processing",category))


ggplot(kk4_res.f2, aes(x = Count, y = reorder(Description, Count))) +  # 注意：x和y应该交换，因为我们要创建水平条形图
  geom_bar(aes(fill = LogP), stat = "identity") +  # 简化aes映射，因为x和y已经在全局aes中设置了
  ggtitle("KEGG") +
  facet_wrap(~category, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_fill_gradient(low="#d5cabd",high = "#dc0000") +  # 图例标题
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),  # 如果需要，可以隐藏x轴标题
        axis.title.y = element_text("Count"))


res = results(dds2, pAdjustMethod = "bonferroni")

## apply variance stabilizing transformation
v = vst(dds2, blind=FALSE)
vsted = assay(v)
## Plot PCA of VST values
DESeq2::plotPCA(v, intgroup="Condition")+
  theme_bw()

## Define the input genes, and use clusterProfiler::bitr to convert the ID.
sig = subset(res, padj<0.05)
cand.entrez = clusterProfiler::bitr(rownames(sig), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

## Perform enrichment analysis (ORA)
pway = ReactomePA::enrichPathway(gene = cand.entrez)
pwayGO = clusterProfiler::enrichGO(cand.entrez, ont = "BP", OrgDb = org.Hs.eg.db)

## Convert to SYMBOL
pway = setReadable(pway, OrgDb=org.Hs.eg.db)
pwayGO = setReadable(pwayGO, OrgDb=org.Hs.eg.db)

## Store the similarity
pway = enrichplot::pairwise_termsim(pway)

## Define including samples
incSample = rownames(subset(meta, Condition=="T"))

allEntrez = clusterProfiler::bitr(rownames(res), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
res$SYMBOL <- rownames(res)
lfc <- merge(data.frame(res), allEntrez, by="SYMBOL")
lfc <- lfc[order(lfc$log2FoldChange, decreasing=TRUE),]
geneList <- lfc$log2FoldChange
names(geneList) <- lfc$ENTREZID

pwayGSE <- ReactomePA::gsePathway(geneList)
sigpway <- subset(pway@result, p.adjust<0.06)
paste(mean(sigpway$Count), sd(sigpway$Count))


barplot(pway, showCategory = 15)
#使用bngeneplot函数绘图
bngeneplot(results = pway, exp = vsted, pathNum = 17)
#Change the label for the better readability.
bngeneplot(results = pway, exp = vsted, pathNum = 17, labelSize=7, shadowText=TRUE)

# Show the confidence of direction
bngeneplot(results = pway,
           exp = vsted,
           expSample = incSample,
           pathNum = 13, R = 50, showDir = T,
           convertSymbol = T,
           expRow = "ENSEMBL",
           strThresh = 0.7)

heat=read.csv("HEATMAP.csv")
rownames(heat)=heat[,1]
heat=heat[,-1]

# 加载必要的库（如果有的话，这里不需要特定库）

# 定义稳健标准化函数
robust_standardize <- function(heat, cols = NULL) {
  # 如果cols为NULL，则选择所有数值列
  if (is.null(cols)) {
    cols <- sapply(heat, is.numeric)
    cols <- names(heat)[cols]
  }
  
  # 初始化一个新的数据框来存储标准化后的数据
  standardized_heat <- heat[, !(names(heat) %in% cols), drop = FALSE]  # 保留非数值列
  
  # 对每一列进行稳健标准化
  for (col in cols) {
    median_val <- median(heat[[col]])
    iqr_val <- IQR(heat[[col]])
    standardized_heat[[col]] <- (heat[[col]] - median_val) / iqr_val
  }
  
  # 返回标准化后的数据框
  return(standardized_heat)
}

# 创建一个示例数据框
set.seed(123)
example_heat <- data.frame(
  ID = 1:110,
  Value1 = c(rnorm(100, mean = 500, sd = 200), runif(10, min = 7500, max = 8000)),
  Value2 = c(rnorm(110, mean = 100, sd = 50))
)

# 对数据框进行稳健标准化
standardized_heat <- robust_standardize(heat)
# 假设您的数据框名为df，且已经正确加载到R环境中
# df <- read.csv("your_data.csv")  # 如果数据是从CSV文件中读取的，可以使用这行代码

# 选择要标准化的列（假设是除了第一列以外的所有列）
columns_to_standardize <- names(heat)[-1]  # 排除第一列（GSM编号）

# 对选定的列执行Z-score标准化
df_standardized <- as.data.frame(lapply(heat[columns_to_standardize], scale))

# 将标准化的数据框与GSM编号列合并
df_standardized$gene <- heat$gene

# 查看标准化后的数据框
head(df_standardized)
df_standardized <- as.data.frame(lapply(heat, scale))
rownames(df_standardized)=df_standardized[,13]
df_standardized=df_standardized[,1:12]

pheatmap(df_standardized)


# 加载必要包
library(survival)
library(survminer)
library(tidyverse)

# 1. 数据准备（与之前相同）
cli <- read.csv("time_ov.csv", header=T, sep=",", check.names=F, row.names=1)
cli$time <- cli$time/365

data <- read.table("TCGA_OV_TPM.txt", header=T, sep="\t", check.names=F, row.names=1)
data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), 
               dimnames=list(rownames(data), colnames(data)))
data <- t(data)
rownames(data) <- substr(rownames(data), 1, 12)
rownames(data) <- gsub('[.]', '-', rownames(data))

sameSample <- intersect(row.names(data), row.names(cli))
data <- data[sameSample, ]
cli <- cli[sameSample, ]
rt <- cbind(cli, data)

# 2. 获取基因列表（排除带点号的基因）
genes <- colnames(rt)[3:ncol(rt)]  
genes <- genes[!grepl("\\.", genes)]  # 移除含点号的基因名

# 3. 批量生存分析并记录结果
results <- data.frame(
  Gene = character(),
  pValue = numeric(),
  HR = numeric(),
  HR_lower = numeric(),
  HR_upper = numeric(),
  stringsAsFactors = FALSE
)

for(gene in genes) {
  tryCatch({
    # 分组（中位数）
    gene_exp <- rt[, gene]
    if(all(is.na(gene_exp))) next
    
    cutoff <- quantile(gene_exp, probs=0.5, na.rm=TRUE)
    rt$group <- ifelse(gene_exp > cutoff, "High", "Low")
    
    # 生存分析
    diff <- survdiff(Surv(time, state) ~ group, data=rt)
    pValue <- 1 - pchisq(diff$chisq, df=1)
    
    # Cox回归计算HR值
    cox_fit <- coxph(Surv(time, state) ~ group, data=rt)
    hr_summary <- summary(cox_fit)
    hr <- hr_summary$coefficients[2]  # HR值
    hr_ci <- hr_summary$conf.int[1, c(3,4)]  # 95% CI
    
    # 记录结果
    results <- rbind(results, data.frame(
      Gene = gene,
      pValue = pValue,
      HR = hr,
      HR_lower = hr_ci[1],
      HR_upper = hr_ci[2]
    ))
    
  }, error=function(e) {
    message(paste("Error in gene", gene, ":", e$message))
  })
}

# 4. 多重检验校正
results$adj.p <- p.adjust(results$pValue, method="BH")

# 5. 筛选显著基因（p < 0.05）并排序
significant_genes <- results %>%
  filter(pValue < 0.05) %>%
  arrange(pValue)

# 6. 输出CSV文件
write.csv(significant_genes, "significant_genes_p0.05.csv", row.names=FALSE)

# 7. 打印汇总信息
cat(paste("\nFound", nrow(significant_genes), "significant genes (p < 0.05)\n"))
cat(paste("Top 5 most significant genes:\n"))
print(head(significant_genes, 5))



# 定义目标基因
target_gene <- "MEK1"
group=ifelse(rt[,target_gene]> quantile(rt[,target_gene],seq(0,1,1/2)) [2],"High","Low")

  # 生存分析
  fit <- survfit(Surv(time, state) ~ group, data = rt)
  
  # 绘图
  surv_plot <- ggsurvplot(
    fit, data = rt,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = c("red", "blue"),
    title = paste("Survival by", target_gene, "Expression"),
    legend.title = "Expression",
    legend.labs = c("High", "Low")
  )

  group=ifelse(rt[,"PEX3"]> quantile(rt[,"PEX3"],seq(0,1,1/2)) [2],"High","Low")
  rt[,"group"]=group
  length=length(levels(factor(group)))
  
  diff=survdiff(Surv(time,state) ~group,data=rt)
  pValue=1-pchisq(diff$chisq,df=length-1)

  pValue=paste0("p=",sprintf("%.03f",pValue))
  fit=survfit(Surv(time,state) ~group,data=rt)  
    
surPlot=ggsurvplot(fit,
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="FGF7 experssion",
                     legend.labs=levels(factor(rt[,"group"])),
                     legend=c(0.8,0.8),
                     font.legend=10,
                     xlab="Time(years)",
                     break.time.by=2,
                     palette=bioCol,
                     surv.median.line="hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25
  )


# 检查该基因是否在显著列表中
if(target_gene %in% significant_genes$Gene){
  
  # 分组（中位数）
  rt$group <- ifelse(rt[, target_gene] > median(rt[, target_gene], "High", "Low")
                     
                     # 生存分析
                     fit <- survfit(Surv(time, state) ~ group, data = rt)
                     
                     # 绘图
                     surv_plot <- ggsurvplot(
                       fit, data = rt,
                       pval = TRUE,
                       pval.method = TRUE,
                       conf.int = FALSE,
                       risk.table = TRUE,
                       palette = c("red", "blue"),
                       title = paste("Survival by", target_gene, "Expression"),
                       legend.title = "Expression",
                       legend.labs = c("High", "Low")
                     )
                     
                     # 保存图片
                     ggsave(paste0(target_gene, "_survival_curve.png"), 
                            plot = print(surv_plot),
                            width = 8, height = 6)
                     
                     print(paste("成功生成", target_gene, "的生存曲线"))
} else {
  print(paste(target_gene, "不在显著基因列表中（p ≥ 0.05）"))
}
