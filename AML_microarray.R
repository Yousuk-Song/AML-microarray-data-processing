# 필요한 패키지 로드
library(affy)         # 마이크로어레이 데이터 로드
library(oligo)        # Affymetrix 마이크로어레이 데이터 처리
library(limma)        # 정규화 및 분석
library(pheatmap)     # 히트맵 시각화
library(ggplot2)      # PCA 및 클러스터링 시각화
library(biomaRt)      # 유전자 annotation
library(cluster)      # Silhouette Score 계산

# 데이터 경로 설정
data_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578"

# .CEL.gz 파일 리스트 가져오기
cel_files <- list.files(path = data_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)

# Affymetrix 데이터 로드 및 정규화
raw_data <- read.celfiles(filenames = cel_files)       # .CEL 파일 로드
norm_data <- rma(raw_data)                             # RMA 정규화 수행 (log2 스케일)
exprs_matrix <- exprs(norm_data)

# Non-AMKL 샘플 필터링
amkl_samples <- c("GSM2601197_SM01.CEL.gz", "GSM2601209_SM13.CEL.gz", "GSM2601226_SM30.CEL.gz", 
                  "GSM2601227_SM31.CEL.gz", "GSM2601228_SM32.CEL.gz", "GSM2601234_SM38.CEL.gz", 
                  "GSM2601236_SM40.CEL.gz", "GSM2601240_SM44.CEL.gz", "GSM2601243_SM47.CEL.gz")
non_amkl_indices <- !(colnames(exprs_matrix) %in% amkl_samples)
combat_data <- exprs_matrix[, non_amkl_indices]

# 유전자 annotation 설정
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Probe ID와 Gene Symbol 매핑
probe_annotation <- getBM(
  attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
  filters = "affy_hg_u133_plus_2",
  values = rownames(combat_data),
  mart = mart
)

# Probe ID와 발현값 매핑
annotated_exprs <- merge(
  data.frame(probe_id = rownames(combat_data), combat_data),
  probe_annotation,
  by.x = "probe_id",
  by.y = "affy_hg_u133_plus_2"
)

# Gene별 평균 계산
gene_exprs <- aggregate(
  . ~ hgnc_symbol,
  data = annotated_exprs[, -1], # probe_id 제거
  FUN = mean
)

# Gene 이름을 행 이름으로 설정
rownames(gene_exprs) <- gene_exprs$hgnc_symbol
gene_exprs <- gene_exprs[, -1] # hgnc_symbol 열 제거

# 상위 변이 유전자 선택
n_top_genes <- 5000  # 사용할 유전자 수 설정
gene_variances <- apply(gene_exprs, 1, var)  # 각 유전자의 분산 계산
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:n_top_genes]  # 상위 n개 유전자 선택

# 선택된 유전자만 사용
filtered_gene_exprs <- gene_exprs[top_genes, ]

mhc_class2_genes <- c(
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
  "HLA-DPA1", "HLA-DPA2", "HLA-DPB1", "HLA-DPB2",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "CD74"
)

# MHC Class II 관련 유전자 데이터 필터링
mhc_genes_filtered_exprs <- filtered_gene_exprs[rownames(filtered_gene_exprs) %in% mhc_class2_genes, ]
# 필터링된 유전자 데이터가 비어있는지 확인
if (nrow(mhc_genes_filtered_exprs) == 0) {
  stop("No matching MHC Class II genes found in the expression data.")
}

# PCA 수행
mhc_pca <- prcomp(t(mhc_genes_filtered_exprs), scale. = TRUE)

# 설명된 분산 비율 계산
mhc_explained_variance <- mhc_pca$sdev^2 / sum(mhc_pca$sdev^2)
mhc_cumulative_variance <- cumsum(mhc_explained_variance)

# 누적 분산이 80% 이상을 차지하는 첫 번째 PC 찾기
pc_80 <- which(mhc_cumulative_variance >= 0.8)[1]

cat("Cumulative Variance Explained:\n", mhc_cumulative_variance, "\n")
cat("First PC capturing >=80% variance:", pc_80, "\n")

# 설명된 분산 비율 시각화
variance_df <- data.frame(
  PC = 1:length(mhc_explained_variance),
  ExplainedVariance = mhc_explained_variance,
  CumulativeVariance = mhc_cumulative_variance
)

variance_plot <- ggplot(variance_df, aes(x = PC)) +
  geom_bar(aes(y = ExplainedVariance), stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = CumulativeVariance, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = CumulativeVariance), color = "red", size = 2) +
  geom_vline(xintercept = pc_80, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  annotate("text", x = pc_80, y = 0.8, label = paste("PC", pc_80, "\n80%"), vjust = -1, color = "darkgreen") +
  labs(
    title = "Explained Variance and Cumulative Variance by PCs",
    x = "Principal Component (PC)",
    y = "Variance Explained",
    caption = "Bar: Explained Variance, Line: Cumulative Variance"
  ) +
  theme_minimal()

print(variance_plot)

# Elbow Method를 이용해 최적 클러스터 수 탐색
max_k <- 15  # 최대 클러스터 수 설정
wcss <- sapply(1:max_k, function(k) {
  kmeans_result <- kmeans(mhc_pca$x[, 1:pc_80], centers = k, nstart = 25)
  kmeans_result$tot.withinss
})

# Elbow Plot 생성
elbow_plot <- ggplot(data.frame(K = 1:max_k, WCSS = wcss), aes(x = K, y = WCSS)) +
  geom_line(color = "blue") +
  geom_point(size = 3, color = "red") +
  labs(
    title = "Elbow Method for Optimal Clusters",
    x = "Number of Clusters (K)",
    y = "WCSS (Within-Cluster Sum of Squares)"
  ) +
  theme_minimal()
print(elbow_plot)

# 사용자가 원하는 기준으로 클러스터 수 설정
num_clusters <- 2  # 기본값 설정, Elbow Plot 및 PCA 결과에 따라 수정 가능

# 클러스터링 수행
mhc_kmeans <- kmeans(mhc_pca$x[, 1:pc_80], centers = num_clusters, nstart = 25)

# 클러스터 정렬 기준 설정 (예: PC1의 평균값으로 클러스터 정렬)
cluster_means <- tapply(mhc_pca$x[, 1], mhc_kmeans$cluster, mean)
cluster_order <- order(cluster_means)

# 클러스터 레이블 재정렬
mhc_kmeans$cluster <- factor(mhc_kmeans$cluster, levels = cluster_order, labels = paste0("Cluster", 1:num_clusters))

# PCA 결과와 클러스터 추가
mhc_pca_df <- data.frame(
  PC1 = mhc_pca$x[, 1],
  PC2 = mhc_pca$x[, 2],
  Samples = colnames(mhc_genes_filtered_exprs),
  Cluster = mhc_kmeans$cluster
)

# 'Cluster2' 샘플 필터링
cluster2_samples <- mhc_pca_df$Samples[mhc_pca_df$Cluster == "Cluster2"]
cluster2_exprs <- mhc_genes_filtered_exprs[, cluster2_samples, drop = FALSE]

# 이후 동일한 과정 반복
cluster2_pca <- prcomp(t(cluster2_exprs), scale. = TRUE)

cluster2_kmeans <- kmeans(cluster2_pca$x[, 1:cluster2_pc_80], centers = 2, nstart = 25)

# 클러스터 정렬 및 재정렬
cluster2_means <- tapply(cluster2_pca$x[, 1], cluster2_kmeans$cluster, mean)
cluster2_order <- order(cluster2_means)

cluster2_kmeans$cluster <- factor(cluster2_kmeans$cluster, levels = cluster2_order, labels = paste0("Cluster", 1:2))

# PCA 데이터프레임 업데이트
cluster2_pca_df <- data.frame(
  PC1 = cluster2_pca$x[, 1],
  PC2 = cluster2_pca$x[, 2],
  Samples = cluster2_samples,
  Cluster = cluster2_kmeans$cluster
)

# PCA 시각화
mhc_pca_plot <- ggplot(mhc_pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Samples)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 0.5, size = 3) +
  labs(
    title = "PCA Plot of MHC Class II Genes",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  ) +
  theme_minimal()
print(mhc_pca_plot)


# 히트맵 생성
mhc_cluster_annotation <- data.frame(Sample = colnames(mhc_genes_filtered_exprs), Cluster = mhc_pca_df$Cluster)
rownames(mhc_cluster_annotation) <- mhc_cluster_annotation$Sample

pheatmap(
  mhc_genes_filtered_exprs,
  scale = "row",
  annotation_col = mhc_cluster_annotation,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "MHC Class II Gene Expression Clustering"
)



# Sample과 CellLine 매핑
cell_line_data <- data.frame(
  Sample = c(
    "GSM2601198_SM02.CEL.gz", "GSM2601199_SM03.CEL.gz", "GSM2601200_SM04.CEL.gz",
    "GSM2601201_SM05.CEL.gz", "GSM2601202_SM06.CEL.gz", "GSM2601203_SM07.CEL.gz",
    "GSM2601204_SM08.CEL.gz", "GSM2601205_SM09.CEL.gz", "GSM2601206_SM10.CEL.gz",
    "GSM2601207_SM11.CEL.gz", "GSM2601208_SM12.CEL.gz", "GSM2601210_SM14.CEL.gz",
    "GSM2601211_SM15.CEL.gz", "GSM2601212_SM16.CEL.gz", "GSM2601213_SM17.CEL.gz",
    "GSM2601214_SM18.CEL.gz", "GSM2601215_SM19.CEL.gz", "GSM2601216_SM20.CEL.gz",
    "GSM2601217_SM21.CEL.gz", "GSM2601218_SM22.CEL.gz", "GSM2601219_SM23.CEL.gz",
    "GSM2601220_SM24.CEL.gz", "GSM2601221_SM25.CEL.gz", "GSM2601222_SM26.CEL.gz",
    "GSM2601223_SM27.CEL.gz", "GSM2601224_SM28.CEL.gz", "GSM2601225_SM29.CEL.gz",
    "GSM2601229_SM33.CEL.gz", "GSM2601230_SM34.CEL.gz", "GSM2601231_SM35.CEL.gz",
    "GSM2601232_SM36.CEL.gz", "GSM2601233_SM37.CEL.gz", "GSM2601235_SM39.CEL.gz",
    "GSM2601237_SM41.CEL.gz", "GSM2601238_SM42.CEL.gz", "GSM2601239_SM43.CEL.gz",
    "GSM2601241_SM45.CEL.gz", "GSM2601242_SM46.CEL.gz", "GSM2601244_SM48.CEL.gz"
  ),
  CellLine = c(
    "ML-2", "OCI-AML2", "OCI-AML3", "OCI-M1", "OCI-M2", "SKM-1", "SIG-M5", "PLB-985",
    "MOLM-13", "EOL-1", "HNT-34", "U937", "THP-1", "KG-1", "HL60/MX1", "MOLM14",
    "MV4;11", "GDM-1", "KU812", "TUR", "K562", "TF-1a", "MM1", "MEG-A2", "Kasumi-1",
    "NOMO-1", "HL60/MX2", "TF-1", "HEL", "KG1a", "Kasumi-6", "HL60", "SKNO-1",
    "MM6", "AML-193", "Kasumi-3", "NB-4", "OCI-AML5", "AP-1060"
  )
)

cell_line_map <- setNames(cell_line_data$CellLine, cell_line_data$Sample)

# Cluster annotation에서 Sample 이름을 CellLine으로 변경
mhc_cluster_annotation <- data.frame(
  CellLine = cell_line_map[colnames(mhc_genes_filtered_exprs)],
  Cluster = mhc_pca_df$Cluster
)
rownames(mhc_cluster_annotation) <- mhc_cluster_annotation$CellLine

# Column 이름을 CellLine으로 변경
colnames(mhc_genes_filtered_exprs) <- cell_line_map[colnames(mhc_genes_filtered_exprs)]

# 히트맵 생성
pheatmap(
  mhc_genes_filtered_exprs,
  scale = "row",
  annotation_col = mhc_cluster_annotation,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "MHC Class II Gene Expression Clustering (Cell Line)"
)

# PCA 시각화에서도 CellLine 반영
mhc_pca_df$CellLine <- cell_line_map[mhc_pca_df$Samples]

mhc_pca_plot <- ggplot(mhc_pca_df, aes(x = PC1, y = PC2, color = Cluster, label = CellLine)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 0.5, size = 3) +
  labs(
    title = "PCA Plot of MHC Class II Genes (Cell Line)",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  ) +
  theme_minimal()
print(mhc_pca_plot)




# Cluster2만 가지고 한 번 더 clustering

# Cluster2 샘플 필터링
cluster2_samples <- mhc_pca_df$Samples[mhc_pca_df$Cluster == "Cluster2"]
# cluster2_samples를 Cell Line 이름으로 변환
cluster2_cell_lines <- cell_line_data$CellLine[match(cluster2_samples, cell_line_data$Sample)]

cluster2_exprs <- mhc_genes_filtered_exprs[,cluster2_cell_lines, drop = FALSE]

# PCA 수행
cluster2_pca <- prcomp(t(cluster2_exprs), scale. = TRUE)

# 설명된 분산 비율 계산
cluster2_explained_variance <- cluster2_pca$sdev^2 / sum(cluster2_pca$sdev^2)
cluster2_cumulative_variance <- cumsum(cluster2_explained_variance)

# 누적 분산이 80% 이상을 차지하는 첫 번째 PC 찾기
cluster2_pc_80 <- which(cluster2_cumulative_variance >= 0.8)[1]

cat("Cumulative Variance Explained:\n", cluster2_cumulative_variance, "\n")
cat("First PC capturing >=80% variance:", cluster2_pc_80, "\n")


# 누적 분산과 설명된 분산 비율 시각화
cluster2_variance_df <- data.frame(
  PC = 1:length(cluster2_explained_variance),
  ExplainedVariance = cluster2_explained_variance,
  CumulativeVariance = cluster2_cumulative_variance
)

variance_plot <- ggplot(cluster2_variance_df, aes(x = PC)) +
  geom_bar(aes(y = ExplainedVariance), stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = CumulativeVariance, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = CumulativeVariance), color = "red", size = 2) +
  geom_vline(xintercept = cluster2_pc_80, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  annotate("text", x = cluster2_pc_80, y = 0.8, label = paste("PC", cluster2_pc_80, "\n80%"), 
           vjust = -1, color = "darkgreen") +
  labs(
    title = "Explained and Cumulative Variance by PCs (Cluster2)",
    x = "Principal Component (PC)",
    y = "Variance Explained",
    caption = "Bar: Explained Variance, Line: Cumulative Variance"
  ) +
  theme_minimal()

print(variance_plot)

# Elbow Method로 최적 클러스터 수 탐색
max_k <- 15  # 최대 클러스터 수 설정
cluster2_wcss <- sapply(1:max_k, function(k) {
  kmeans_result <- kmeans(cluster2_pca$x[, 1:cluster2_pc_80], centers = k, nstart = 25)
  kmeans_result$tot.withinss
})

# Elbow Plot 생성
elbow_plot <- ggplot(data.frame(K = 1:max_k, WCSS = cluster2_wcss), aes(x = K, y = WCSS)) +
  geom_line(color = "blue") +
  geom_point(size = 3, color = "red") +
  labs(
    title = "Elbow Method for Optimal Clusters (Cluster2)",
    x = "Number of Clusters (K)",
    y = "WCSS (Within-Cluster Sum of Squares)"
  ) +
  theme_minimal()
print(elbow_plot)

# 최적 클러스터 수 설정 (사용자 입력 또는 Elbow Plot 참조)
optimal_clusters <- 4  # Elbow Plot 보고 수정 가능

# K-means 클러스터링 수행
cluster2_kmeans <- kmeans(cluster2_pca$x[, 1:cluster2_pc_80], centers = optimal_clusters, nstart = 25)

# PCA 결과에 클러스터 추가
cluster2_pca_df <- data.frame(
  PC1 = cluster2_pca$x[, 1],
  PC2 = cluster2_pca$x[, 2],
  Samples = cluster2_samples,
  Cluster = factor(cluster2_kmeans$cluster, labels = paste0("Cluster", 1:optimal_clusters))
)

# PCA 시각화
cluster2_pca_plot <- ggplot(cluster2_pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Samples)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 0.5, size = 3) +
  labs(
    title = "PCA Plot of MHC Class II Genes - Cluster2",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  ) +
  theme_minimal()
print(cluster2_pca_plot)

# 히트맵 생성
cluster2_annotation <- data.frame(Sample = cluster2_samples, Cluster = cluster2_pca_df$Cluster)
rownames(cluster2_annotation) <- cluster2_annotation$Sample

#cluster2_samples를 Cell Line 이름으로 변환
cluster2_cell_lines <- cell_line_data$CellLine[match(cluster2_samples, cell_line_data$Sample)]

# cluster2_annotation 업데이트 (Cell Line 이름으로 변경)
cluster2_annotation <- data.frame(
  CellLine = cluster2_cell_lines,
  Cluster = cluster2_kmeans$cluster
)
rownames(cluster2_annotation) <- cluster2_annotation$CellLine

# 클러스터별 색상 정의 (리스트로 설정)
cluster2_annotation_colors <- list(
  Cluster = c("Cluster1" = "red", "Cluster2" = "darkgreen", "Cluster3" = "skyblue", "Cluster4" = "purple")
)
pheatmap(
  cluster2_exprs,
  scale = "row",
  annotation_col = cluster2_annotation,
  annotation_colors = cluster2_annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "MHC Class II Gene Expression Clustering - Cluster2"
)

library(cluster)
silhouette_scores <- silhouette(cluster2_kmeans$cluster, dist(t(cluster2_exprs)))
plot(silhouette_scores, main = "Silhouette Plot")

# Silhouette 값이 0.1보다 작은 샘플의 인덱스 확인
low_silhouette_indices <- which(silhouette_scores[, "sil_width"] < 0.1)

# 해당 샘플의 이름 추출 (colnames(cluster2_exprs) 사용)
low_silhouette_samples <- colnames(cluster2_exprs)[low_silhouette_indices]

# 결과 출력
print(low_silhouette_samples)

# 문제가 되는 샘플 제외
filtered_samples <- setdiff(colnames(cluster2_exprs), low_silhouette_samples)
cluster2_exprs_filtered <- cluster2_exprs[, filtered_samples, drop = FALSE]

# 클러스터링 다시 수행
cluster2_pca_filtered <- prcomp(t(cluster2_exprs_filtered), scale. = TRUE)
cluster2_kmeans_filtered <- kmeans(cluster2_pca_filtered$x[, 1:cluster2_pc_80], centers = 4, nstart = 25)

# Silhouette Score 재확인
filtered_silhouette_scores <- silhouette(cluster2_kmeans_filtered$cluster, dist(t(cluster2_exprs_filtered)))
plot(filtered_silhouette_scores, main = "Silhouette Plot (Filtered Samples)")



# AML 전체 top 5000 gene으로 만든 cluster
aml_5000_cluster_annotation <- data.frame(
  Sample = c(
    "GSM2601198_SM02.CEL.gz", "GSM2601199_SM03.CEL.gz", "GSM2601200_SM04.CEL.gz",
    "GSM2601201_SM05.CEL.gz", "GSM2601202_SM06.CEL.gz", "GSM2601203_SM07.CEL.gz",
    "GSM2601204_SM08.CEL.gz", "GSM2601205_SM09.CEL.gz", "GSM2601206_SM10.CEL.gz",
    "GSM2601207_SM11.CEL.gz", "GSM2601208_SM12.CEL.gz", "GSM2601210_SM14.CEL.gz",
    "GSM2601211_SM15.CEL.gz", "GSM2601212_SM16.CEL.gz", "GSM2601213_SM17.CEL.gz",
    "GSM2601214_SM18.CEL.gz", "GSM2601215_SM19.CEL.gz", "GSM2601216_SM20.CEL.gz",
    "GSM2601217_SM21.CEL.gz", "GSM2601218_SM22.CEL.gz", "GSM2601219_SM23.CEL.gz",
    "GSM2601220_SM24.CEL.gz", "GSM2601221_SM25.CEL.gz", "GSM2601222_SM26.CEL.gz",
    "GSM2601223_SM27.CEL.gz", "GSM2601224_SM28.CEL.gz", "GSM2601225_SM29.CEL.gz",
    "GSM2601229_SM33.CEL.gz", "GSM2601230_SM34.CEL.gz", "GSM2601231_SM35.CEL.gz",
    "GSM2601232_SM36.CEL.gz", "GSM2601233_SM37.CEL.gz", "GSM2601235_SM39.CEL.gz",
    "GSM2601237_SM41.CEL.gz", "GSM2601238_SM42.CEL.gz", "GSM2601239_SM43.CEL.gz",
    "GSM2601241_SM45.CEL.gz", "GSM2601242_SM46.CEL.gz", "GSM2601244_SM48.CEL.gz"
  ),
  AML_5000_cluster = c(
    "Cluster1", "Cluster1", "Cluster1", "Cluster2", "Cluster2", "Cluster1",
    "Cluster1", "Cluster1", "Cluster1", "Cluster1", "Cluster3", "Cluster1",
    "Cluster1", "Cluster3", "Cluster1", "Cluster1", "Cluster1", "Cluster1",
    "Cluster2", "Cluster1", "Cluster2", "Cluster2", "Cluster1", "Cluster2",
    "Cluster1", "Cluster1", "Cluster1", "Cluster2", "Cluster2", "Cluster3",
    "Cluster2", "Cluster2", "Cluster1", "Cluster1", "Cluster1", "Cluster3",
    "Cluster1", "Cluster1", "Cluster1"
  )
)

# MHC 클러스터 및 Cell Line 정보 추가
aml_5000_cluster_annotation$MHC_cluster <- mhc_cluster_annotation$Cluster[match(aml_5000_cluster_annotation$Sample, mhc_cluster_annotation$Sample)]
aml_5000_cluster_annotation$CellLine <- cell_line_data$CellLine[match(aml_5000_cluster_annotation$Sample, cell_line_data$Sample)]

print(aml_5000_cluster_annotation)

# CSV 파일로 저장
output_path <- file.path(data_dir, "DEG_results", "aml_5000_vs_mhc_cluster_annotation.csv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE) # 디렉토리가 없는 경우 생성
write.csv(aml_5000_cluster_annotation, file = output_path, row.names = FALSE)

cat("AML 5000 vs MHC cluster annotation saved to:", output_path, "\n")




# 모든 샘플의 유전자별 평균 발현값 계산
overall_mean_expr <- rowMeans(mhc_genes_filtered_exprs, na.rm = TRUE)

# DEG 분석을 위한 데이터프레임 준비
deg_results <- list()

# MHC 클러스터별 DEG 분석
for (cluster in unique(aml_5000_cluster_annotation$MHC_cluster)) {
  # 해당 클러스터의 샘플 선택
  cluster_samples <- aml_5000_cluster_annotation$Sample[aml_5000_cluster_annotation$MHC_cluster == cluster]
  
  # 클러스터 내 발현값 추출
  cluster_exprs <- mhc_genes_filtered_exprs[, colnames(mhc_genes_filtered_exprs) %in% cluster_samples, drop = FALSE]
  
  # 클러스터 평균 계산 (샘플 부족 체크)
  if (ncol(cluster_exprs) < 2) {
    cat("Skipping cluster", cluster, "due to insufficient samples.\n")
    next
  }
  
  cluster_mean_expr <- rowMeans(cluster_exprs, na.rm = TRUE)
  
  # DEG 분석 (t-검정)
  deg_df <- data.frame(
    Gene = rownames(mhc_genes_filtered_exprs),
    ClusterMean = cluster_mean_expr,
    OverallMean = overall_mean_expr,
    FoldChange = log2(cluster_mean_expr / overall_mean_expr), # Fold Change 계산
    PValue = apply(mhc_genes_filtered_exprs, 1, function(expr_values) {
      t.test(
        expr_values[colnames(mhc_genes_filtered_exprs) %in% cluster_samples],
        expr_values[!colnames(mhc_genes_filtered_exprs) %in% cluster_samples]
      )$p.value
    })
  )
  
  # FDR 보정
  deg_df$AdjustedPValue <- p.adjust(deg_df$PValue, method = "fdr")
  
  # Upregulation/Downregulation 추가
  deg_df$Regulation <- ifelse(deg_df$FoldChange > 0, "Upregulated", "Downregulated")
  
  # 유의미한 DEG 필터링 (p < 0.05)
  significant_deg <- subset(deg_df, AdjustedPValue < 0.05)
  
  # 결과 저장
  deg_results[[cluster]] <- significant_deg
}

# 결과 병합
deg_combined <- do.call(rbind, lapply(names(deg_results), function(cluster) {
  deg_results[[cluster]]$Cluster <- cluster
  deg_results[[cluster]]
}))

# 결과 출력
if (nrow(deg_combined) > 0) {
  print(head(deg_combined))
} else {
  cat("No significant DEGs found.\n")
}

# DEG 결과 저장
output_deg_path <- file.path(data_dir, "DEG_results", "mhc_cluster_deg_results_with_regulation.csv")
dir.create(dirname(output_deg_path), recursive = TRUE, showWarnings = FALSE) # 디렉토리가 없는 경우 생성
write.csv(deg_combined, file = output_deg_path, row.names = FALSE)

cat("DEG analysis results saved to:", output_deg_path, "\n")

# 클러스터별로 유전자 발현량 정렬 및 CSV 저장
output_cluster_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results"
#dir.create(output_cluster_dir, recursive = TRUE, showWarnings = FALSE) # 디렉토리가 없는 경우 생성

# GSM ID 추출 함수 정의
extract_gsm <- function(sample_name) {
  sub("_SM.*", "", sample_name)  # "_" 이후 부분 제거
}

cell_line_data$GSM <- sapply(cell_line_data$Sample, extract_gsm)

# 클러스터별 데이터 저장
for (cluster in unique(filtered_pca_df$Cluster)) {
  cluster_dir <- file.path(output_cluster_dir, paste0("Cluster_", cluster))
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)  # 디렉토리 생성
  
  # 해당 클러스터의 샘플 선택
  cluster_samples <- filtered_pca_df$Samples[filtered_pca_df$Cluster == cluster]
  
  for (gene in rownames(cd36_genes_filtered_exprs)) {
    # 해당 유전자의 발현량 데이터 추출
    gene_expr <- cd36_genes_filtered_exprs[gene, cluster_samples, drop = FALSE]
    gene_sample_exprs <- data.frame(
      Sample = names(gene_expr),
      Gene = gene,
      Expression = as.numeric(gene_expr)
    )
    
    # CellName 추가 (샘플 이름 그대로 매핑)
    gene_sample_exprs$CellName <- sample_names$CellName[match(gene_sample_exprs$Sample, sample_names$GSM)]
    orubt
    
    # 발현량 기준 내림차순 정렬
    gene_sample_exprs_sorted <- gene_sample_exprs[order(-gene_sample_exprs$Expression), ]
    
    # 클러스터별 유전자 CSV 저장
    gene_file <- file.path(cluster_dir, paste0(gene, "_Expression_in_Cluster_", cluster, ".csv"))
    write.csv(gene_sample_exprs_sorted, file = gene_file, row.names = FALSE)
    
    cat("Saved gene expression data for", gene, "in Cluster", cluster, "to:", gene_file, "\n")
  }
}
