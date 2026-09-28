#这个项目是为了提取featurecounts 比对数量的TPM值
#原始代码来自于featurecounts.R
expr_df <- read.table( "/xtdisk/xueyb_group/zhangby/lyl/data2024/RNAseq/20240128yahs_canu/featurecount",header=T, row.names=1, check.names=F, sep="\t") 
dim(expr_df); names(expr_df) 
head(expr_df[,1:7])
#提取基因信息,featureCounts前几列
featureCounts_meta <- expr_df[,1:5] ;head(featureCounts_meta)
expr_df <- expr_df[,6:ncol(expr_df)]
## 保存counts矩阵
write.table(expr_df, "merged.Counts.txt",quote=F, sep="\t", row.names=T, col.names=T )

prefix <-"SV"   #设置输出文件前缀名

# ----- TPM计算 ------

# 基因长度，目标基因的外显子长度之和除以1000，单位是Kb，不是bp
kb <- featureCounts_meta$Length / 1000 
rpk <- expr_df / kb   #每千碱基reads (“per million” scaling factor) 长度标准化
tpm <- t(t(rpk)/colSums(rpk) * 1000000)  # 每百万缩放因子 (“per million” scaling factor ) 深度标准化
avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))
colnames(tpm) <- c("ID","Anther","Pistil","Petal","Leaf")
# 保存
write.table(avg_tpm, paste0(prefix,"_avg_tpm.xls"),quote=F, sep="\t", row.names=T, col.names=T )
write.table(tpm, paste0(prefix,"_tpm"), quote=F, sep="\t", row.names=T, col.names=T )