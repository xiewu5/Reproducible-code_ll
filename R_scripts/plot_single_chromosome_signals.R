## ===============================
## 0. 加载包
## ===============================
library(CMplot)

setwd('')
## ===============================
## 1. 读入 Chr06 per-site Fst
## ===============================
data <- read.table(
  "../Chr06/Chr06.weir.fst",
  header = TRUE,
  stringsAsFactors = FALSE
)

## ===============================
## 2. 整理为 CMplot 标准格式
## ===============================
Fst <- data.frame(
  SNP = paste(data$CHROM, data$POS, sep = ":"),
  Chromosome = data$CHROM,
  Position = data$POS,
  Fst = data$WEIR_AND_COCKERHAM_FST
)

## 去掉 NA / Inf
Fst <- Fst[is.finite(Fst$Fst), ]

## Weir & Cockerham Fst 允许负值
## 论文中通常直接置 0
Fst$Fst[Fst$Fst < 0] <- 0

## ===============================
## 3. 阈值 & y 轴
## ===============================
cutoff <- quantile(Fst$Fst, 0.95, na.rm = TRUE)
ymax <- max(Fst$Fst, na.rm = TRUE)

## ===============================
## 4. 画 Chr06 Fst 曼哈顿
## ===============================
CMplot(
  Fst,
  type = "p",
  plot.type = "m",
  LOG10 = FALSE,
  col = "steelblue4",   # 单染色体 = 单色，论文更干净
  cex = 0.35,
  band = 0.6,
  ylab = "Fst",
  ylim = c(0, ymax),
  threshold = cutoff,
  #threshold.col = "red",
  threshold.lty = 2,
  threshold.lwd = 1,
  amplify = FALSE,
  file.output = TRUE,
  file = "pdf",
  dpi = 600,
  file.name = "Chr06_Fst_per_site"
)

