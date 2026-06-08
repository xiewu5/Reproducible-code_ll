#!/usr/bin/bash
#SBATCH --job-name=fastp_clean      # 任务名称
#SBATCH --partition=sommi          # 计算节点分区
#SBATCH --output=fastp_output.log  # 日志文件
#SBATCH --error=fastp_error.log    # 错误文件
#SBATCH --nodelist=compare-0-0     # 指定节点（按需修改）
#SBATCH -n 16                      # 使用16核（根据需求调整）

# 解除系统限制
ulimit -s unlimited
ulimit -l unlimited

# 创建输出目录
mkdir -p cleandata

# 运行fastp（直接指定你的文件名）
fastp \
  -i 5-L_raw_1.fq.gz -I 5-L_raw_2.fq.gz \
  -o cleandata/5-L_clean_1.fq.gz -O cleandata/5-L_clean_2.fq.gz \
  -l 150 -w 16 \
  --html cleandata/5-L_fastp_report.html \
  --json cleandata/5-L_fastp_report.json
