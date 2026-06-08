# 保存为 calc_HwF.R
awk 'NR>1' IND_diversity.hwe > tmp.hwe
R --no-save -q << 'EOF'
hwe <- read.table('tmp.hwe', header=FALSE,
                  col.names=c('CHR','POS','OBS','E','Chi','P','P_def','P_exc'))
obs <- do.call(rbind, strsplit(as.character(hwe$OBS), '/'))
exp <- do.call(rbind, strsplit(as.character(hwe$E),   '/'))
obs <- as.data.frame(apply(obs, 2, as.numeric))
exp <- as.data.frame(apply(exp, 2, as.numeric))
colnames(obs) <- colnames(exp) <- c('HOM1','HET','HOM2')

valid <- rowSums(obs) > 0 & rowSums(exp) > 0 & exp$HET > 0
obs <- obs[valid, ]; exp <- exp[valid, ]

Ho  <- mean(obs$HET / rowSums(obs))
He  <- mean(exp$HET / rowSums(exp))
Fis <- mean(1 - Ho/He, na.rm = TRUE)

cat(sprintf('IND\tHo=%.3f\tHe=%.3f\tFis=%.3f\n', Ho, He, Fis))
writeLines(sprintf('IND\t%.3f\t%.3f\t%.3f', Ho, He, Fis), 'IND_HwF.txt')
EOF
rm tmp.hwe
