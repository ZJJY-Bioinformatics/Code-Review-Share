# 普氏分析
require(tidyverse)
library(phyloseq)
require(microbiomeutilities)
require(ggpubr)
require(microbiome)
require(vegan)

rm(list = ls())

#----------
data1 = read.delim('data/16s_table.xls', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
data2 = read.delim('data/metabo_table.xls', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#----------
# 需要先对数据计算其样本间的距离
# 一般物种群落使用bray-curtis距离，而环境因子使用欧式距离
amp_mat = data1
gen_mat = data2

nid = intersect(colnames(amp_mat),colnames(gen_mat))

amp_mat = amp_mat[,nid]
gen_mat = gen_mat[,nid]

amp.dist <- vegdist(t(amp_mat)) # 默认Bray-Curtis
gen.dist <- vegdist(t(gen_mat))

# 可以换其他降维方式-------
mds.amp <- monoMDS(amp.dist)
mds.gen <- monoMDS(gen.dist)

pro.s.e <- procrustes(mds.amp,mds.gen, symmetric = TRUE)
summary(pro.s.e)

pdf("the_sample_residual.pdf")
plot(pro.s.e, kind = 2)
dev.off()

write.csv(residuals(pro.s.e),"the_sample_residual.csv")

set.seed(1)
pro.s.e_t <- protest(mds.amp,mds.gen, permutations = 999)

# 偏差平方和（M2统计量）
pro.s.e_t$ss
# 对应p值结果
pro.s.e_t$signif

library(ggplot2)
# 获得x和y轴的坐标及旋转过的坐标
Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
Pro_X <- data.frame(pro.s.e$rotation)

# 绘图
ggplot(data=Pro_Y) +
  geom_segment(aes(x = X1, y = X2,
                   xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2),
               # geom_segment 绘制两点间的直线
               arrow = arrow(length = unit(0, 'cm')),
               color = "grey50", size = 0.8) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2,
                   xend = MDS1, yend = MDS2),
               arrow = arrow(length = unit(0.3, 'cm')),
               color = "grey50", size = 0.8) +
  geom_point(aes(X1, X2, color = "data2"), size = 3, shape = 16) +
  geom_point(aes(MDS1, MDS2, color = "data1"), size = 3, shape = 16) +
  theme(panel.grid = element_blank(), # 绘制背景
        panel.background = element_rect(color = 'black',
                                        fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between 16S and Metagenome") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  scale_color_manual(values = c(
                                "data2" = "#ff6fab", 
                                "data1" = "#442dd4"))+
  annotate('text', label = paste0('   Procrustes analysis:\n   M2 = ',round(pro.s.e_t$ss,3),', p-value = ',pro.s.e_t$signif,'\n\n'),
           x = -Inf, y = -Inf, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=14,colour = "black",
                                  hjust = 0.5,face = "bold"))
ggsave("Procrustes.pdf",width = 8,height = 6)
#---------------



