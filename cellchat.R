###############################allcellchat--normal#########
data1=GetAssayData(ovar_Normal,assay = "RNA",layer = "data")
table(ovar_Normal@meta.data$celltype.main)
table(ovar_Normal@meta.data$orig.ident)
cellchat=createCellChat(data1,meta = ovar_Normal@meta.data,
                        group.by = "celltype.main")
CellChatDB=CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB=CellChatDB
cellchat=subsetData(cellchat)
future::plan("multisession",workers=10)
#options(future.globals.maxSize = 1000 * 1024^2)  # 1000 MiB
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat@idents <- droplevels(cellchat@idents)
cellchat=computeCommunProb(cellchat,raw.use = F)

cellchat_normal=filterCommunication(cellchat,min.cells=10)
cellchat_normal=computeCommunProbPathway(cellchat_normal)
df.net.cancer=subsetCommunication(cellchat_normal)
write.csv(df.net.cancer,"Gene.csv",row.names = F)
df.netP.cancer=subsetCommunication(cellchat_normal,slot.name = "netP")
write.csv(df.netP.cancer,"Pathway.csv",row.names = F)
cellchat_normal=aggregateNet(cellchat_normal)
save(cellchat_normal,file="ovar_Normal_cellchat.Rdata")
groupSize=as.numeric(table(cellchat_normal@idents))

pdf("Normal_NetVisual_overview_all.pdf",width = 8,height = 6)
par(xpd=TRUE)
netVisual_circle(cellchat_normal@net$count,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Number of interactions")
netVisual_circle(cellchat_normal@net$weight,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Interaction weights/strength")
dev.off()

pdf("Normal_NetVisual_overview_split.pdf",width = 6,height = 5)
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

pdf("Normal_NetVisual_pathways_circle.pdf",width=6,height=5)
for(pathways.show in mypathways1){
  par(xpd=T)
  netVisual_aggregate(cellchat_normal,signaling = pathways.show,layout="circle")
}
dev.off()

pdf("Normal_NetVisual_pathways_chord.pdf",width = 10,height = 8)
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
ggsave("Normal_CCI_all.pdf",p,width = 15,height=20,limitsize = F)

p=netVisual_bubble(cellchat_normal,sources.use = c("CD4+_T","CD8_TEM","CD8_TNaive","CD8_TEMRA"),
                   targets.use =c("CD4+_T","CD8_TEM","CD8_TNaive","CD8_TEMRA"),remove.isolate = F)
ggsave("CCI_part.pdf",p,width = 15,height=20,limitsize = F)

# p=netVisual_bubble(cellchat_normal,sources.use ="T",targets.use = c("CD8+ T","CD4+ T","Tregs"),remove.isolate = F)
# ggsave("CCI_T.pdf",p,width = 8,height=20,limitsize = F)
# 
# p=netVisual_bubble(cellchat_normal,sources.use =c("T","B","NK"),targets.use = c("CD8+ T","CD4+ T","Tregs","Plasma"),signaling=c("CCL","TNF"),remove.isolate = F)
# ggsave("CCI_T_pathway.pdf",p,width = 8,height=6,limitsize = F)

pairLR.use=c("CCL3_CCR1","CCL4_CCR5","CCL5_CCR3","TNF_TNFRSF1A","TNF_TNFRSF1B","GZMA_F2R")
pairLR.use=data.frame(interaction_name=pairLR.use)
p=netVisual_bubble(cellchat_normal,sources.use = c("CD4+_T","CD8_TEM","CD8_TNaive","CD8_TEMRA"),
                   targets.use =c("CD4+_T","CD8_TEM","CD8_TNaive","CD8_TEMRA"),
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
netAnalysis_signalingRole_heatmap(cellchat_normal2,pattern = "outgoing",signaling = c("CCL","LCK","CD6","MHC-I","CD80"))
netAnalysis_signalingRole_heatmap(cellchat_normal2,pattern = "incoming",signaling = c("CCL","LCK","CD6","MHC-I","CD80"))

