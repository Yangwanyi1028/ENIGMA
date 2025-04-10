# 清空环境变量
rm(list=ls())

# 加载包
library(ggplot2)
library(ggExtra)
library(vegan)
library(ggthemes)

# 设置工作路径
setwd('/Users/yangwanyi/Downloads/ENIGMA_STUDY_DATA/2pcoa/')

# 加载数据
# 加载种水平物种丰度表，并设置列名和分隔符
spe <- read.csv("species_abun_df.csv",row.names = 1)
spe <- t(spe)
# 计算所有物种相对丰度
data1 <- spe / apply(spe, 2, sum)
data2 <- t(data1)# 翻转

# 加载分组表group
group <- read.csv('group_df.csv',stringsAsFactors = FALSE)

# 进行PCoA分析
# vegdist函数，计算距离；method，选择距离类型
# 可选euclidean、manhattan、jaccard
bray <- vegdist(data2, method = 'bray',na.rm = T)
# euclidean <- vegdist(data, method = 'euclidean')
# manhattan <- vegdist(data, method = 'manhattan')
# jaccard <- vegdist(data, method = 'jaccard')
# 计算加权bray-curtis距离
# dune_bray <- vegdist(data, method="bray", binary=F)
bray <- as.matrix(bray)
write.table(bray,"bray-crutis.txt",sep = "\t") 

# 对非加权距离进行PCoA分析,k等于3,前三个维度。
# 可选k = (nrow(data) - 1)
pcoa <- cmdscale(bray, k = 3, eig = T)

# 提取样本点坐标
pcoa_data <- data.frame({pcoa$point})

# 提取列名
pcoa_data$Sample_ID <- rownames(pcoa_data)
names(pcoa_data)[1:3] <- paste0("PCoA", 1:3)

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig = pcoa$eig
eig_percent <- round(pcoa$eig/sum(eig)*100,1)

# 提取样本点坐标
poi = pcoa$points
poi = as.data.frame(poi)

#为样本点坐标添加分组信息
pcoa_result <- cbind(pcoa_data,group ) # PCoA值和tax表/OTU表合并
head(pcoa_result)

# 进行置换多元（因素）方差分析（PERMANOVA）
# 基于bray-curtis距离进行
dune.div <- adonis2(data2 ~ Region, data = group, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis
# 绘制简单的pcoa图
p = ggplot(pcoa_result, aes(x=PCoA1, y=PCoA2, color=Region)) +
  geom_point(aes(color=Region),size=5)+ # ,shape=Region
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       caption = dune_adonis)  +# 也可用title
  scale_colour_manual(values = c("#85a894","#3f6b7e","#fce4a8"))+
  theme(legend.position = c(0.95,0.1),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.tictle=element_text(hjust=0.5),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text = element_text(color = "black",size=10))+
  geom_hline(aes(yintercept=0), colour="#BEBEBE", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="#BEBEBE", linetype="dashed")

p = p + stat_ellipse(data=pcoa_result, 
                     geom = "polygon", 
                     level=0.9, 
                     linetype = 2, 
                     linewidth=0.5, 
                     aes(fill=Region), 
                     alpha=0.3, 
                     show.legend = T) +
  scale_fill_manual(values = c("#85a894","#3f6b7e","#fce4a8")) # 这里添加或修改颜色


pdf(file="bray_pcoa_region.pdf", width=10, height=10)
ggMarginal(
  p,
  type =c('density'),
  margins = 'both',
  size = 3.5,
  groupColour = F,
  groupFill = T
) 
dev.off()

# 进行置换多元（因素）方差分析（PERMANOVA）
# 基于bray-curtis距离进行



noncdai_sample <- group[group$Group=='CD',]$study_id
group <- group[group$Group!='CD',]
data2 <- data2[!rownames(data2) %in% noncdai_sample, ]
pcoa_data <- pcoa_data[!rownames(pcoa_data) %in% noncdai_sample,]

#为样本点坐标添加分组信息
pcoa_result <- cbind(pcoa_data,group ) # PCoA值和tax表/OTU表合并
head(pcoa_result)

dune.div <- adonis2(data2 ~ Group, data = group, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("adonis R2: ", round(dune.div$R2, 2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis

# 绘制简单的PCoA图
p <- ggplot(pcoa_result, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(aes(color = Group), size = 5) + #, shape = Group
  labs(
    x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
    y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
    caption = dune_adonis  # 也可用title
  ) +
  scale_colour_manual(values = c("firebrick","#f96b0a","#90a5a6")) +
  theme(
    legend.position = c(0.9, 0.1),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    axis.text = element_text(color = "black", size = 10)
  ) +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")

# 添加椭圆
p <- p + stat_ellipse(
  data = pcoa_result,
  geom = "polygon",
  level = 0.9,
  linetype = 2,
  linewidth = 0.5,
  aes(fill = Group),
  alpha = 0.3,
  show.legend = TRUE
) +
  scale_fill_manual(values = c("firebrick","#f96b0a","#90a5a6"))  # 这里添加或修改颜色


pdf(file="bray_pcoa_group.pdf", width=10, height=10) 
ggMarginal(
  p,
  type =c('density'),
  margins = 'both',
  size = 3.5,
  groupColour = F,
  groupFill = T
) 
dev.off()

