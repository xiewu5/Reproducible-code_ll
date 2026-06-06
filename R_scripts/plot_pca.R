setwd("")
# 读取 PCA 结果
pca <- read.table("../input/PCA_plink.out.eigenvec", header = FALSE, stringsAsFactors = FALSE)
colnames(pca) <- c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

# 读取样本性别分组文件
pop <- read.table("../input/sample.pop", header = FALSE, stringsAsFactors = FALSE)
colnames(pop) <- c("IID","Sex")

# 合并数据
dat <- merge(pca, pop, by = "IID")

# 读取特征值，计算 PC1 和 PC2 的解释率
eig <- scan("../input/PCA_plink.out.eigenval")
pc1_var <- round(eig[1] / sum(eig) * 100, 2)
pc2_var <- round(eig[2] / sum(eig) * 100, 2)

# 加载作图包
library(ggplot2)

# 设定分组顺序（可选）
dat$Sex <- factor(dat$Sex, levels = c("XX", "XY"))

# 画图
p <- ggplot(dat, aes(x = PC1, y = PC2, color = Sex, shape = Sex)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("XX" = 16, "XY" = 17)) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

# 显示图
print(p)


