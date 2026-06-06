############################################################
## 0. 加载包
############################################################
library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)

############################################################
## 1. 设置路径与读入数据
############################################################
setwd("")
bg_file <- "../go_input_data/background_withGO.tsv"
deg_file <- "../go_input_data/DEGs_male_vs_female13.noNewGene.tsv"

# 读取背景并整理为长表
bg_raw <- read.table(bg_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, sep = "\t")
term2gene <- bg_raw %>%
  pivot_longer(cols = -V1, names_to = "tmp", values_to = "GO") %>%
  filter(GO != "" & !is.na(GO)) %>%
  select(GO, V1) %>%
  distinct()

# 读取差异基因
deg_df <- read.delim(deg_file, header = TRUE, stringsAsFactors = FALSE)
deg_genes <- unique(deg_df$id)
gene_logfc <- setNames(deg_df$logFC, deg_df$id)

############################################################
## 2. 准备描述信息
############################################################
term2name <- AnnotationDbi::select(
  GO.db, keys = unique(term2gene$GO), 
  columns = c("TERM", "ONTOLOGY"), keytype = "GOID"
) %>% dplyr::rename(GO = GOID, Description = TERM)

############################################################
## 3. 运行富集分析
############################################################
go_res <- enricher(
  gene = deg_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name[, c("GO", "Description")],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

############################################################
## 4. 生成 13 列标准结果
############################################################
if (is.null(go_res) || nrow(as.data.frame(go_res)) == 0) {
  cat("⚠️ 警告：没有显著富集结果。\n")
} else {
  # 生成完整结果表
  go_output <- as.data.frame(go_res) %>%
    left_join(distinct(term2name[, c("GO", "ONTOLOGY")]), by = c("ID" = "GO")) %>%
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
  write.table(go_output, file = "GO_enrichment_results_13columns.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  ############################################################
  ## 5. 提取 BP、CC、MF 各前 5 个显著条目
  ############################################################
  go_top5_each <- go_output %>%
    group_by(ONTOLOGY) %>%            # 按分类分组
    arrange(p.adjust) %>%             # 组内按显著性排序
    slice_head(n = 5) %>%             # 每组取前 5 行
    ungroup()                         # 取消分组
  
  # 保存精简集
  write.table(go_top5_each, file = "GO_enrichment_top5_per_Ontology.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("✅ 任务完成！\n")
  cat("1. 完整结果已保存至：GO_enrichment_results_13columns.tsv\n")
  cat("2. 分组 Top5 已保存至：GO_enrichment_top5_per_Ontology.tsv\n")
}

