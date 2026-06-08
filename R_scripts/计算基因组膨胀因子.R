# 读取 GEMMA 结果文件


gwas <- read.table(
  "C:/Users/13799/Desktop/R/谢武-r脚本/计数基因组膨胀因子/长吻鮠/no_chr8_indel.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)
# 提取 P 值
p <- gwas$p_wald

# 基本清洗
p <- p[!is.na(p) & p > 0 & p <= 1]

# 转换为卡方统计量
chisq <- qchisq(1 - p, df = 1)

# 计算 λ
lambda <- median(chisq) / qchisq(0.5, df = 1)

lambda

