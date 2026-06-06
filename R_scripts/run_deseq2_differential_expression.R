## ============================================================
##  性别差异分析（DESeq2）
##  输入文件：
##    counts.final.int.tsv  <- 整数 counts 矩阵
##    metadata.tsv          <- 样本信息，第二列 Sex 为分组
## ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## 1. 读入数据 =====================================================
counts <- read.delim("", check.names = FALSE, row.names = 1)
meta   <- read.delim("", stringsAsFactors = FALSE)

## 保证 counts 与 meta 样本顺序一致
meta <- meta[match(colnames(counts), meta$Sample), ]
stopifnot(all(colnames(counts) == meta$Sample))

## 设置分组：female 为对照组
meta$Group <- factor(meta$Sex, levels = c("female", "male"))

## 2. 过滤低表达基因 ================================================
# 至少 3 个样本 counts >= 10 才保留
keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, ]

## 3. DESeq2 差异分析 ==============================================
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ Group
)

dds <- DESeq(dds)

# 获取性别差异结果
res <- results(dds, contrast = c("Group", "male", "female"))

# log2FC 收缩
res <- lfcShrink(dds, coef = "Group_male_vs_female", type = "ashr")

## 4. 结果整理 ======================================================
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
res_df <- res_df[order(res_df$padj), ]

## 5. 显著基因筛选 ================================================
padj_threshold  <- 0.05
log2FC_threshold <- 1

diff_genes <- subset(res_df, padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold)

## 6. 输出结果 ======================================================
write.csv(res_df, "AllGenes_male_vs_female.csv", quote = FALSE, row.names = TRUE)
write.csv(diff_genes, "DEGs_male_vs_female.csv", quote = FALSE, row.names = TRUE)

## 上调 / 下调基因分别输出
up_genes   <- subset(diff_genes, log2FoldChange > 0)
down_genes <- subset(diff_genes, log2FoldChange < 0)

write.csv(up_genes, "UP_male_vs_female.csv", quote = FALSE, row.names = TRUE)
write.csv(down_genes, "DOWN_male_vs_female.csv", quote = FALSE, row.names = TRUE)

## 7. 可选：快速检查差异基因数量 ==================================
cat("总基因数:", nrow(res_df), "\n")
cat("显著上调基因数:", nrow(up_genes), "\n")
cat("显著下调基因数:", nrow(down_genes), "\n")
