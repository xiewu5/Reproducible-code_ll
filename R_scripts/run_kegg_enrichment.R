############################################################
## 0. 加载包
############################################################
library(tidyverse)
library(clusterProfiler)

############################################################
## 1. 设置路径与读入数据
############################################################
setwd("")
bg_file <- "../go_input_data/background_withkegg.tsv"
deg_file <- "../go_input_data/DEGs_male_vs_female13.noNewGene.tsv"

# 读入背景并整理为长表 (TERM2GENE)
bg_raw <- read.table(bg_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, sep = "\t")
term2gene <- bg_raw %>%
  pivot_longer(cols = -V1, names_to = "tmp", values_to = "ko") %>%
  filter(ko != "" & !is.na(ko)) %>%
  filter(str_detect(ko, "ko")) %>%
  select(ko, V1) %>%
  distinct()

# 读入差异基因，保留 logFC 用于计算 zScore
deg_df <- read.delim(deg_file, header = TRUE, stringsAsFactors = FALSE)
deg_genes <- unique(deg_df$id)
gene_logfc <- setNames(deg_df$logFC, deg_df$id)

############################################################
## 2. 获取描述信息 (TERM2NAME)
############################################################
message("正在联网获取最新 KEGG 描述信息...")
kegg_list_raw <- clusterProfiler:::kegg_list("pathway")
term2name <- data.frame(
  ko = gsub("path:", "", kegg_list_raw$from) %>% gsub("map", "ko", .),
  Description = kegg_list_raw$to,
  stringsAsFactors = FALSE
)

############################################################
## 3. 运行富集分析
############################################################
kegg_res <- enricher(
  gene = deg_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

############################################################
## 4. 生成 13 列标准结果
############################################################
if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
  cat("⚠️ 警告：没有筛选到显著富集的 KEGG 通路。\n")
} else {
  # 生成完整结果表
  kegg_output <- as.data.frame(kegg_res) %>%
    mutate(ONTOLOGY = "KEGG Pathway") %>%
    mutate(
      k = as.numeric(gsub("/.*", "", GeneRatio)), 
      n = as.numeric(gsub(".*/", "", GeneRatio)), 
      M = as.numeric(gsub("/.*", "", BgRatio)),   
      N = as.numeric(gsub(".*/", "", BgRatio)),
      RichFactor = k / M,
      FoldEnrichment = (k / n) / (M / N)
    ) %>%
    rowwise() %>%
    mutate(
      up = sum(gene_logfc[unlist(strsplit(geneID, "/"))] > 0, na.rm = TRUE),
      down = sum(gene_logfc[unlist(strsplit(geneID, "/"))] < 0, na.rm = TRUE),
      zScore = (up - down) / sqrt(Count)
    ) %>%
    ungroup() %>%
    select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, RichFactor, 
           FoldEnrichment, zScore, pvalue, p.adjust, qvalue, geneID, Count)
  
  # 保存全集
  write.table(kegg_output, file = "KEGG_enrichment_results_13columns.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  ############################################################
  ## 5. 提取 KEGG 最显著的前 5 个条目
  ############################################################
  kegg_top5 <- kegg_output %>%
    arrange(p.adjust) %>%   # 按显著性排序
    slice_head(n = 5)       # 取前 5 个
  
  # 保存 Top5 结果
  write.table(kegg_top5, file = "KEGG_enrichment_top5_significant.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("✅ KEGG 任务完成！\n")
  cat("1. 完整结果已保存至：KEGG_enrichment_results_13columns.tsv\n")
  cat("2. 最显著 Top5 已保存至：KEGG_enrichment_top5_significant.tsv\n")
}

