data <- read.table('different_gene.txt',sep=',',header = T)

suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))

#描述上调、下调关系
data$logFC[which(data$logFC < 0)] <- 'down'
data$logFC[which(data$logFC != 'down' & data$logFC != 0)] <- 'up'

#检查基因symbol
sample_unique <- unique(as.character(data[,1]))
mat <- match(sample_unique,data[,1])
data <- data[mat,]
data <- na.omit(data)

#提取gene symbol做GO富集分析
de<-data[,1]
de <- as.character(de)
rm(mat,sample_unique)

#symbol转换为ENTREZID
symbol <- select(org.Hs.eg.db,  #数据库
                 keys=de, #我们提供的ID
                 columns="ENTREZID", #需要在数据库里取出哪些ID信息
                 keytype="SYMBOL" #我们提供的ID属于哪种类型
                 )

de <- symbol$ENTREZID
#GO enrichment reuqires  clusterProfiler package
BP <- enrichGO(de, organism="human", ont="BP", pvalueCutoff=0.01)
CC <- enrichGO(de, organism="human", ont="CC", pvalueCutoff=0.01)
MF <- enrichGO(de, organism="human", ont="MF", pvalueCutoff=0.01)

BP_analysis <- as.data.frame(summary(BP))
CC_analysis <- as.data.frame(summary(CC))
MF_analysis <- as.data.frame(summary(MF))
rm(BP,CC,MF,DOSEEnv,GOSemSimEnv,ICEnv,SemSimCache,clusterProfilesEnv,de)


#构建up/down向量
up_type <- symbol$ENTREZID[which(data$logFC == 'up')]
down_type <- symbol$ENTREZID[which(data$logFC == 'down')]


#取出排名前十的GO富集结果
BP_analysis <- BP_analysis[order(BP_analysis$pvalue, decreasing = TRUE),]
if(length(BP_analysis$pvalue) > 10){
  BP_analysis <- BP_analysis[1:10,c(2,8)]
} else {
  BP_analysis <- BP_analysis[,c(2,8)]
}

#构建循环对BP_analysis中每一行的geneID进行分割
bplist <- alist()
for(i in 1:nrow(BP_analysis)){
  bplist[[i]] <- strsplit(as.character(BP_analysis$geneID)[i],'/')[[1]]
}

#分别统计上下调基因个数
up_cha <- as.character()
down_cha <- as.character()
for(i in 1:length(bplist)){
  up_cha[i] <- as.numeric(length(intersect(bplist[[i]],up_type)))
  down_cha[i] <- as.numeric(length(intersect(bplist[[i]],down_type)))
}
BP_analysis_up <- cbind(BP_analysis,as.data.frame(up_cha))
BP_analysis_down <- cbind(BP_analysis,as.data.frame(down_cha))

#上下调分别命名
colnames(BP_analysis_up)[3] <- 'up_down'
BP_analysis_up$geneID <- rep('up',each=nrow(BP_analysis_up))
colnames(BP_analysis_down)[3] <- 'up_down'
BP_analysis_down$geneID <- rep('down',each=nrow(BP_analysis_down))

#合并表格（上下调合在一起）
BP_analysis <- rbind(BP_analysis_up,BP_analysis_down)
BP_analysis <- cbind(BP_analysis,as.matrix(rep('BP',each=nrow(BP_analysis))))
colnames(BP_analysis)[4] <- 'type'
rm(BP_analysis_up,BP_analysis_down)



CC_analysis <- CC_analysis[order(CC_analysis$pvalue, decreasing = TRUE),]
if(length(CC_analysis$pvalue) > 10){
  CC_analysis <- CC_analysis[1:10,c(2,8)]
} else {
  CC_analysis <- CC_analysis[,c(2,8)]
}
bplist <- alist()
for(i in 1:nrow(CC_analysis)){
  bplist[[i]] <- strsplit(as.character(CC_analysis$geneID)[i],'/')[[1]]
}
up_cha <- as.character()
down_cha <- as.character()
for(i in 1:length(bplist)){
  up_cha[i] <- as.numeric(length(intersect(bplist[[i]],up_type)))
  down_cha[i] <- as.numeric(length(intersect(bplist[[i]],down_type)))
}
CC_analysis_up <- cbind(CC_analysis,as.data.frame(up_cha))
CC_analysis_down <- cbind(CC_analysis,as.data.frame(down_cha))
colnames(CC_analysis_up)[3] <- 'up_down'
CC_analysis_up$geneID <- rep('up',each=nrow(CC_analysis_up))
colnames(CC_analysis_down)[3] <- 'up_down'
CC_analysis_down$geneID <- rep('down',each=nrow(CC_analysis_down))
CC_analysis <- rbind(CC_analysis_up,CC_analysis_down)
CC_analysis <- cbind(CC_analysis,as.matrix(rep('CC',each=nrow(CC_analysis))))
colnames(CC_analysis)[4] <- 'type'
rm(CC_analysis_up,CC_analysis_down)


MF_analysis <- MF_analysis[order(MF_analysis$pvalue, decreasing = TRUE),]
if(length(MF_analysis$pvalue) > 10){
  MF_analysis <- MF_analysis[1:10,c(2,8)]
} else {
  MF_analysis <- MF_analysis[,c(2,8)]
}
bplist <- alist()
for(i in 1:nrow(MF_analysis)){
  bplist[[i]] <- strsplit(as.character(MF_analysis$geneID)[i],'/')[[1]]
}
up_cha <- as.character()
down_cha <- as.character()
for(i in 1:length(bplist)){
  up_cha[i] <- as.numeric(length(intersect(bplist[[i]],up_type)))
  down_cha[i] <- as.numeric(length(intersect(bplist[[i]],down_type)))
}
MF_analysis_up <- cbind(MF_analysis,as.data.frame(up_cha))
MF_analysis_down <- cbind(MF_analysis,as.data.frame(down_cha))
colnames(MF_analysis_up)[3] <- 'up_down'
MF_analysis_up$geneID <- rep('up',each=nrow(MF_analysis_up))
colnames(MF_analysis_down)[3] <- 'up_down'
MF_analysis_down$geneID <- rep('down',each=nrow(MF_analysis_down))
MF_analysis <- rbind(MF_analysis_up,MF_analysis_down)
MF_analysis <- cbind(MF_analysis,as.matrix(rep('MF',each=nrow(MF_analysis))))
colnames(MF_analysis)[4] <- 'type'
rm(MF_analysis_up,MF_analysis_down,up_cha,down_cha,bplist,down_type,up_type,i,data,symbol)

# CC_analysis <- CC_analysis[1:10,c(1,2,9)]
# MF_analysis <- MF_analysis[1:10,c(1,2,9)]
# 
# BP_analysis <- cbind(BP_analysis,as.matrix(rep('BP',each=10)))
# BP_analysis <- cbind(BP_analysis,as.matrix(rep(c('up','down'),each=5)))
# colnames(BP_analysis)[[4]][1] <- 'type'
# colnames(BP_analysis)[[5]][1] <- 'regulate'
# 
# CC_analysis <- cbind(CC_analysis,as.matrix(rep('CC',each=10)))
# CC_analysis <- cbind(CC_analysis,as.matrix(rep(c('up','down'),each=5)))
# colnames(CC_analysis)[[4]][1] <- 'type'
# colnames(CC_analysis)[[5]][1] <- 'regulate'
# 
# MF_analysis <- cbind(MF_analysis,as.matrix(rep('MF',each=10)))
# MF_analysis <- cbind(MF_analysis,as.matrix(rep(c('up','down'),each=5)))
# colnames(MF_analysis)[[4]][1] <- 'type'
# colnames(MF_analysis)[[5]][1] <- 'regulate'
# 
# analysis <- rbind(BP_analysis,CC_analysis,MF_analysis)

GO_analysis <- rbind(BP_analysis,CC_analysis,MF_analysis)
GO_analysis$up_down <- as.numeric(as.character(GO_analysis$up_down))
rm(BP_analysis,CC_analysis,MF_analysis)

p <- ggplot(GO_analysis,aes(Description, up_down, fill=geneID)) + theme_bw() + ylab('GO enrichment number')

#p+geom_bar(stat='identity') +coord_flip()+ labs(y="-log(P-value)",x="KEGG enrichment")
#P + geom_bar(position='dodge')

p + facet_wrap(~type,scales='free_x') + theme(axis.text.x=element_text(angle=80,hjust=1),axis.title.x=NULL) + 
  geom_bar(position='dodge',stat="identity") + theme(panel.grid =element_blank())

