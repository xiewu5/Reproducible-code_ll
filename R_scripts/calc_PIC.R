awk 'NR>1' AUS_NG_diversity.hwe > tmp.hwe
R --no-save -q << 'EOF'
hwe <- read.table('tmp.hwe', header=FALSE,
                  col.names=c('CHR','POS','OBS','E','Chi','P','P_def','P_exc'))
obs <- do.call(rbind, strsplit(as.character(hwe$OBS), '/'))
obs <- as.data.frame(apply(obs, 2, as.numeric))
colnames(obs) <- c('nAA','nAB','nBB')
n   <- rowSums(obs)
p   <- (2*obs$nAA + obs$nAB) / (2*n)
q   <- 1 - p

## 过滤：观测>0 且 非固定位点
valid <- n > 0 & p > 0 & q > 0
pic <- 1 - (p[valid]^2 + q[valid]^2) - 2 * p[valid]^2 * q[valid]^2
mean_pic <- mean(pic, na.rm = TRUE)

cat(sprintf('AUS_NG\tPIC=%.3f\n', mean_pic))
writeLines(sprintf('AUS_NG\t%.3f', mean_pic), 'AUS_NG_PIC.txt')
EOF
rm tmp.hwe
