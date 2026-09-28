library(pheatmap)
library(RColorBrewer)
library(dplyr)
setwd("D:/Rscript/02_Sc03_rnaseq")
countData <- read.csv('Slocus_tpm03', header = TRUE, sep = "\t")
# 对所有数值列取自然对数，并创建新列（或覆盖原列）

df_scaled <- countData  # 复制原数据框
numeric_cols <- sapply(df_scaled, is.numeric)  # 找出数值列
df_scaled[numeric_cols] <- scale(countData[numeric_cols])  # 替换为标准化后的值

expr_matrix <- countData %>%
  mutate(across(2:5, log1p))
rownames(expr_matrix) <- expr_matrix$Geneid
expr_matrix$Geneid <- NULL 
rownames(countData) <- countData$Geneid
countData$Geneid  <- NULL
rownames(df_scaled) <- df_scaled$Geneid
df_scaled$Geneid  <- NULL
# 查看数据前几行
head(expr_matrix)
# 转置矩阵
expr_t <- t(df_scaled)

# 使用pheatmap绘制
pheatmap(expr_matrix,
         main = "Gene Expression Heatmap",
         fontsize_row = 6,       # 行标签字体大小
         fontsize_col = 8,      # 列标签字体大小
         cluster_rows = FALSE,   # 基因顺序不变
         cluster_cols = FALSE,    # 列（样本）仍可聚类（可选）
         show_rownames = TRUE,   # 显示基因名（基因多时可设为FALSE）
         show_colnames = TRUE)
# 定义颜色渐变
my_colors <- colorRampPalette(c("white", "#d77f3d", "#721c39"))(10)

# 保存为PDF
pdf("heatmap01.pdf", width = 15, height = 4)
pheatmap(expr_t,
         main = "Gene Expression Heatmap",
         fontsize_row = 6,       # 行标签字体大小
         fontsize_col = 8,      # 列标签字体大小
         color = my_colors,
         cluster_rows = FALSE,   # 基因顺序不变
         cluster_cols = FALSE,    # 列（样本）仍可聚类（可选）
         show_rownames = TRUE,   # 显示基因名（基因多时可设为FALSE）
         show_colnames = TRUE)
dev.off()

# 保存为PNG（高分辨率）
png("heatmap01.png", width = 2500, height = 1000, res = 300)
pheatmap(expr_t,
         main = "Gene Expression Heatmap",
         fontsize_row = 6,       # 行标签字体大小
         fontsize_col = 8,      # 列标签字体大小
         color = my_colors,
         cluster_rows = FALSE,   # 基因顺序不变
         cluster_cols = FALSE,    # 列（样本）仍可聚类（可选）
         show_rownames = TRUE,   # 显示基因名（基因多时可设为FALSE）
         show_colnames = TRUE)
dev.off()

