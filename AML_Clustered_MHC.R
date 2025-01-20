 # 필요한 패키지 로드
library(affy)
library(oligo)
library(limma)
library(pheatmap)
library(ggplot2)
library(biomaRt)
library(cluster)
library(clusterProfiler)

# 데이터 경로 및 결과 저장 경로 설정
data_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578"
output_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/DEG_result/Array_DEG"

# 결과 디렉토리 생성
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# .CEL.gz 파일 리스트 가져오기
cel_files <- list.files(path = data_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)

# Affymetrix 데이터 로드 및 정규화
raw_data <- read.celfiles(filenames = cel_files)
norm_data <- rma(raw_data)
exprs_matrix <- exprs(norm_data)

# Non-AMKL 샘플 필터링
amkl_samples <- c("GSM2601197_SM01.CEL.gz", "GSM2601209_SM13.CEL.gz", "GSM2601226_SM30.CEL.gz", 
                  "GSM2601227_SM31.CEL.gz", "GSM2601228_SM32.CEL.gz", "GSM2601234_SM38.CEL.gz", 
                  "GSM2601236_SM40.CEL.gz", "GSM2601240_SM44.CEL.gz", "GSM2601243_SM47.CEL.gz")
non_amkl_indices <- !(colnames(exprs_matrix) %in% amkl_samples)
combat_data <- exprs_matrix[, non_amkl_indices]

# 유전자 annotation 설정
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
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
  data = annotated_exprs[, -1],
  FUN = mean
)

# Gene 이름을 행 이름으로 설정
rownames(gene_exprs) <- gene_exprs$hgnc_symbol
gene_exprs <- gene_exprs[, -1]

# 상위 변이 유전자 선택
n_top_genes <- 5000
gene_variances <- apply(gene_exprs, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:n_top_genes]
filtered_gene_exprs <- gene_exprs[top_genes, ]

# PCA 수행
filtered_pca <- prcomp(t(filtered_gene_exprs), scale. = TRUE)



# 필요한 패키지 로드
library(cluster)
library(ggplot2)
library(factoextra)

# 클러스터 수 결정 함수
determine_optimal_clusters <- function(data, max_clusters = 10, seed = 123) {
  set.seed(seed)
  
  # Elbow Method
  wss <- sapply(1:max_clusters, function(k) {
    kmeans_result <- kmeans(t(data), centers = k, nstart = 25)
    kmeans_result$tot.withinss
  })
  
  elbow_plot <- data.frame(Clusters = 1:max_clusters, WSS = wss)
  
  elbow_plot_graph <- ggplot(elbow_plot, aes(x = Clusters, y = WSS)) +
    geom_line(size = 1, color = "blue") +
    geom_point(size = 3, color = "red") +
    theme_minimal() +
    labs(
      title = "Elbow Method for Optimal Clusters",
      x = "Number of Clusters",
      y = "Within-Cluster Sum of Squares"
    )
  
  print(elbow_plot_graph)
  
  # Silhouette Method
  silhouette_scores <- sapply(2:max_clusters, function(k) {
    kmeans_result <- kmeans(t(data), centers = k, nstart = 25)
    silhouette_values <- silhouette(kmeans_result$cluster, dist(t(data)))
    mean(silhouette_values[, 3]) # Average silhouette score
  })
  
  silhouette_df <- data.frame(Clusters = 2:max_clusters, Silhouette = silhouette_scores)
  
  silhouette_plot_graph <- ggplot(silhouette_df, aes(x = Clusters, y = Silhouette)) +
    geom_line(size = 1, color = "blue") +
    geom_point(size = 3, color = "red") +
    theme_minimal() +
    labs(
      title = "Silhouette Method for Optimal Clusters",
      x = "Number of Clusters",
      y = "Average Silhouette Score"
    )
  
  print(silhouette_plot_graph)
  
  # Gap Statistic
  gap_stat <- clusGap(
    t(data),
    FUN = kmeans,
    K.max = max_clusters,
    B = 50
  )
  
  gap_stat_plot <- fviz_gap_stat(gap_stat)
  print(gap_stat_plot)
  
  # Optimal Clusters from Each Method
  optimal_clusters_elbow <- which.min(diff(diff(wss))) + 1
  optimal_clusters_silhouette <- which.max(silhouette_scores) + 1
  optimal_clusters_gap <- which.max(gap_stat$Tab[, "gap"])
  
  results <- list(
    Elbow = optimal_clusters_elbow,
    Silhouette = optimal_clusters_silhouette,
    GapStatistic = optimal_clusters_gap
  )
  
  return(results)
}

# 데이터 입력 및 클러스터 수 결정
optimal_clusters <- determine_optimal_clusters(filtered_gene_exprs, max_clusters = 10)
print(optimal_clusters)



# 클러스터링 수행
num_clusters <- 3
dist_matrix <- dist(filtered_pca$x[, 1:18])
clustering <- hclust(dist_matrix, method = "ward.D2")
cluster_groups <- cutree(clustering, k = num_clusters)


# 클러스터 결과를 데이터프레임에 추가
filtered_pca_df <- data.frame(
  PC1 = filtered_pca$x[, 1],
  PC2 = filtered_pca$x[, 2],
  Samples = colnames(filtered_gene_exprs),
  Cluster = factor(cluster_groups, labels = paste0("Cluster", 1:num_clusters))
)
cluster_annotation <- data.frame(Sample = colnames(filtered_gene_exprs), Cluster = filtered_pca_df$Cluster)
rownames(cluster_annotation) <- cluster_annotation$Sample

write.csv(cluster_annotation, "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/DEG_result/Array_DEG/Sample_Cluster_Info.csv", row.names = FALSE)

## 다른 PCA 방법들과 비교!


# 필요한 패키지 로드
library(ConsensusClusterPlus)
library(NMF)
library(cluster)
library(ggplot2)
library(pheatmap)

# 1. PCA 수행 및 기본 클러스터링
filtered_pca <- prcomp(t(filtered_gene_exprs), scale. = TRUE)

# Hierarchical Clustering
dist_matrix <- dist(filtered_pca$x[, 1:18])
hclust_clustering <- hclust(dist_matrix, method = "ward.D2")
hclust_clusters <- cutree(hclust_clustering, k = num_clusters)

# PCA 결과 데이터프레임 생성
filtered_pca_df <- data.frame(
  PC1 = filtered_pca$x[, 1],
  PC2 = filtered_pca$x[, 2],
  Samples = colnames(filtered_gene_exprs),
  HierarchicalCluster = factor(hclust_clusters, labels = paste0("Cluster", 1:num_clusters))
)

# 2. Consensus Clustering
consensus_results <- ConsensusClusterPlus(
  as.matrix(filtered_gene_exprs), 
  maxK = 6,  # 최대 클러스터 수 설정
  reps = 50,  # 반복 횟수
  pItem = 0.8,  # 샘플 비율
  pFeature = 1,  # 특징 비율
  clusterAlg = "hc",  # 클러스터링 알고리즘
  distance = "euclidean",  # 거리 계산 방식
  seed = 123,  # 랜덤 시드
  plot = "pdf",  # 결과를 PDF로 저장
  title = file.path(output_dir, "Consensus_Clustering")  # 결과 저장 디렉토리
)

# Consensus Clustering 결과에서 클러스터 번호 추출
consensus_clusters_raw <- consensus_results[[num_clusters]]$consensusClass

# "Cluster" 접두사를 추가하여 수정
consensus_clusters <- paste0("Cluster", consensus_clusters_raw)


# 3. NMF Clustering
nmf_results <- nmf(filtered_gene_exprs, rank = num_clusters, method = "brunet", seed = 123)
nmf_clusters <- predict(nmf_results)

# 4. K-Means Clustering
set.seed(123)
kmeans_results <- kmeans(t(filtered_gene_exprs), centers = num_clusters, nstart = 25)
kmeans_clusters <- kmeans_results$cluster

# 5. 클러스터링 결과 통합
filtered_pca_df$ConsensusCluster <- factor(consensus_clusters)
filtered_pca_df$NMFCluster <- factor(nmf_clusters)
filtered_pca_df$KMeansCluster <- factor(kmeans_clusters)

# 6. PCA 플롯 생성
# Hierarchical Clustering PCA Plot
pca_plot_hierarchical <- ggplot(filtered_pca_df, aes(x = PC1, y = PC2, color = HierarchicalCluster)) +
  geom_point(size = 4, alpha = 0.8) +   # 점 크기와 투명도 설정
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +  # 샘플 이름 추가
  theme_minimal() +
  labs(
    title = "PCA Plot with Hierarchical Clustering",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  )

# Consensus Clustering PCA Plot
pca_plot_consensus <- ggplot(filtered_pca_df, aes(x = PC1, y = PC2, color = ConsensusCluster)) +
  geom_point(size = 4, alpha = 0.8) +   # 점 크기와 투명도 설정
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +  # 샘플 이름 추가
  theme_minimal() +
  labs(
    title = "PCA Plot with Consensus Clustering",
    x = "PC1",
    y = "PC2",
    color = "Consensus Cluster"
  )

# NMF Clustering PCA Plot
pca_plot_nmf <- ggplot(filtered_pca_df, aes(x = PC1, y = PC2, color = NMFCluster)) +
  geom_point(size = 4, alpha = 0.8) +   # 점 크기와 투명도 설정
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +  # 샘플 이름 추가
  theme_minimal() +
  labs(
    title = "PCA Plot with NMF Clustering",
    x = "PC1",
    y = "PC2",
    color = "NMF Cluster"
  )

# K-Means Clustering PCA Plot
pca_plot_kmeans <- ggplot(filtered_pca_df, aes(x = PC1, y = PC2, color = KMeansCluster)) +
  geom_point(size = 4, alpha = 0.8) +   # 점 크기와 투명도 설정
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +  # 샘플 이름 추가
  theme_minimal() +
  labs(
    title = "PCA Plot with K-Means Clustering",
    x = "PC1",
    y = "PC2",
    color = "K-Means Cluster"
  )


# 7. PCA 플롯 출력
print(pca_plot_hierarchical)
print(pca_plot_consensus)
print(pca_plot_nmf)
print(pca_plot_kmeans)   

# 8. PCA 플롯 저장
ggsave(file.path(output_dir, "PCA_Hierarchical_Clustering.png"), pca_plot_hierarchical, width = 8, height = 6)
ggsave(file.path(output_dir, "PCA_Consensus_Clustering.png"), pca_plot_consensus, width = 8, height = 6)
ggsave(file.path(output_dir, "PCA_NMF_Clustering.png"), pca_plot_nmf, width = 8, height = 6)
ggsave(file.path(output_dir, "PCA_KMeans_Clustering.png"), pca_plot_kmeans, width = 8, height = 6)

cat("PCA plots for all clustering methods saved to:", output_dir, "\n")



library(plotly)
# 3D PCA 플롯 생성
pca_3d_plot <- plot_ly(
  data = filtered_pca_df,
  x = ~PC1, y = ~PC2, z = ~filtered_pca$x[, 3],  # PC1, PC2, PC3 사용
  color = ~ConsensusCluster,                    # 클러스터 정보로 색상 지정
  text = ~Samples,                              # 샘플 이름 추가 (툴팁에 표시)
  type = "scatter3d",                           # 3D 산점도
  mode = "markers",                             # 점으로 표시
  marker = list(size = 6, opacity = 0.8)        # 마커 크기와 투명도
) %>%
  layout(
    title = "3D PCA Plot with Consensus Clustering",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

# 플롯 출력
pca_3d_plot




# 9. Silhouette 분석을 통한 클러스터링 품질 평가
silhouette_scores <- list(
  Hierarchical = silhouette(hclust_clusters, dist(t(filtered_gene_exprs))),
  Consensus = silhouette(consensus_clusters, dist(t(filtered_gene_exprs))),
  NMF = silhouette(as.numeric(nmf_clusters), dist(t(filtered_gene_exprs))),
  KMeans = silhouette(as.numeric(kmeans_clusters), dist(t(filtered_gene_exprs)))
)

# Silhouette 평균값 출력
for (method in names(silhouette_scores)) {
  avg_silhouette <- mean(silhouette_scores[[method]][, 3])
  cat(sprintf("Average silhouette width for %s: %.3f\n", method, avg_silhouette))
}



# 클러스터 결과를 데이터프레임에 추가 (Consensus Clustering 기반)
filtered_pca_df <- data.frame(
  PC1 = filtered_pca$x[, 1],
  PC2 = filtered_pca$x[, 2],
  Samples = colnames(filtered_gene_exprs),
  ConsensusCluster = factor(consensus_clusters, labels = paste0("Cluster", 1:length(unique(consensus_clusters))))
)

# 클러스터 annotation 데이터프레임 생성
cluster_annotation <- data.frame(
  Sample = colnames(filtered_gene_exprs),
  Cluster = filtered_pca_df$ConsensusCluster
)
rownames(cluster_annotation) <- cluster_annotation$Sample

# 저장된 클러스터 정보 확인
print(filtered_pca_df)
print(cluster_annotation)

# DEG 분석 및 저장 (upregulated, downregulated)
pairwise_degs_list <- list()

for (comp in names(pairwise_comparisons)) {
  # case와 control 정의
  case <- pairwise_comparisons[[comp]][[1]]
  control <- pairwise_comparisons[[comp]][[2]]
  
  # 그룹 변수 생성
  group_assignment <- ifelse(cluster_annotation$Cluster == case, "Case",
                             ifelse(cluster_annotation$Cluster %in% control, "Control", NA))
  selected_samples <- !is.na(group_assignment)
  group <- factor(group_assignment[selected_samples])
  
  # 디자인 매트릭스 생성
  design <- model.matrix(~ 0 + group)  # Intercept 없이 그룹만 포함
  colnames(design) <- levels(group)   # "Case", "Control"로 이름 설정
  
  # 대조군 설정
  contrast <- makeContrasts(Case_vs_Control = Case - Control, levels = design)
  
  # limma 분석
  exprs_data <- filtered_gene_exprs[, selected_samples, drop = FALSE]
  fit <- lmFit(exprs_data, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)
  deg_results <- topTable(fit, coef = "Case_vs_Control", number = Inf, sort.by = "p")
  
  # Upregulated 및 Downregulated 유전자 필터링
  upregulated_genes <- deg_results[deg_results$adj.P.Val < 0.01 & deg_results$logFC > 3, ]
  downregulated_genes <- deg_results[deg_results$adj.P.Val < 0.01 & deg_results$logFC < -3, ]
  
  # DEG 결과 저장
  write.csv(deg_results, file = file.path(output_dir, paste0(comp, "_All_DEGs.csv")), row.names = TRUE)
  
  # Upregulated 저장
  if (nrow(upregulated_genes) > 0) {
    write.csv(upregulated_genes, file = file.path(output_dir, paste0(comp, "_Upregulated_DEGs.csv")), row.names = TRUE)
  } else {
    cat(paste0("No upregulated genes for ", comp, ". Skipping save.\n"))
  }
  
  # Downregulated 저장
  if (nrow(downregulated_genes) > 0) {
    write.csv(downregulated_genes, file = file.path(output_dir, paste0(comp, "_Downregulated_DEGs.csv")), row.names = TRUE)
  } else {
    cat(paste0("No downregulated genes for ", comp, ". Skipping save.\n"))
  }
  
  # 결과 저장
  pairwise_degs_list[[comp]] <- list(
    "deg_results" = deg_results,
    "upregulated_genes" = upregulated_genes,
    "downregulated_genes" = downregulated_genes
  )
}

cat("DEG analysis complete. Results saved to output directory.\n")


# 통합 히트맵 생성 및 저장 (upregulated + downregulated)
all_top_genes <- unique(unlist(lapply(pairwise_degs_list, function(x) {
  c(rownames(x$upregulated_genes), rownames(x$downregulated_genes))
})))

combined_heatmap_data <- filtered_gene_exprs[all_top_genes, , drop = FALSE]
write.csv(combined_heatmap_data, file = file.path(output_dir, "Combined_Heatmap_Data.csv"), row.names = TRUE)

pheatmap(
  combined_heatmap_data,
  scale = "row",
  annotation_col = cluster_annotation,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Combined Heatmap of Top DEG Genes"
)

cat("All results (including upregulated and downregulated) have been saved to:", output_dir, "\n")

# 결과 확인
print("Upregulated Genes List:")
print(upregulated_genes_list)

print("Downregulated Genes List:")
print(downregulated_genes_list)


library("ggvenn")
# Upregulated Venn Diagram
upregulated_venn_plot <- ggvenn(
  data = upregulated_genes_list,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),  # 색상 추가
  text_size = 5,
  set_name_size = 4
)
print(upregulated_venn_plot)
# 저장
ggsave(file.path(output_dir, "Upregulated_Venn_Diagram.png"), upregulated_venn_plot, width = 8, height = 6)

# Downregulated Venn Diagram
downregulated_venn_plot <- ggvenn(
  data = downregulated_genes_list,
  fill_color = c("#984EA3", "#FFFF33", "#A65628", "#F781BF"),  # 색상 추가
  text_size = 5,
  set_name_size = 4
)
print(downregulated_venn_plot)
# 저장
ggsave(file.path(output_dir, "Downregulated_Venn_Diagram.png"), downregulated_venn_plot, width = 8, height = 6)

cat("Enhanced Venn Diagrams saved to:", output_dir, "\n")




# 교집합 및 각 그룹 간 겹치는 유전자 확인 함수 (logFC와 adj.P.Val 포함)
extract_overlaps_with_details <- function(genes_list, deg_data_list) {
  # 모든 비교군 쌍 확인
  combn_results <- combn(names(genes_list), 2, function(groups) {
    overlap <- intersect(genes_list[[groups[1]]], genes_list[[groups[2]]])
    if (length(overlap) > 0) {
      overlap_details <- lapply(overlap, function(gene) {
        gene_info <- lapply(groups, function(group) {
          cluster_data <- deg_data_list[[group]]$deg_results
          cluster_row <- cluster_data[rownames(cluster_data) == gene, c("logFC", "adj.P.Val")]
          data.frame(Group = group, Gene = gene, logFC = cluster_row$logFC, adj.P.Val = cluster_row$adj.P.Val)
        })
        do.call(rbind, gene_info)
      })
      do.call(rbind, overlap_details)
    } else {
      return(NULL)  # 교집합이 없으면 NULL 반환
    }
  }, simplify = FALSE)
  
  # 모든 비교군 간 완전한 교집합
  common_genes <- Reduce(intersect, genes_list)
  if (length(common_genes) > 0) {
    common_genes_details <- lapply(common_genes, function(gene) {
      gene_info <- lapply(names(deg_data_list), function(group) {
        cluster_data <- deg_data_list[[group]]$deg_results
        cluster_row <- cluster_data[rownames(cluster_data) == gene, c("logFC", "adj.P.Val")]
        data.frame(Group = group, Gene = gene, logFC = cluster_row$logFC, adj.P.Val = cluster_row$adj.P.Val)
      })
      do.call(rbind, gene_info)
    })
    common_genes_df <- do.call(rbind, common_genes_details)
    common_genes_df$Group1 <- "All"
    common_genes_df$Group2 <- "All"
    combn_results <- c(combn_results, list(common_genes_df))
  }
  
  # NULL 값 제거 및 결과 합치기
  overlap_df <- do.call(rbind, combn_results[!sapply(combn_results, is.null)])
  return(overlap_df)
}

# Upregulated 유전자 간 겹치는 유전자 확인 및 상세 정보 추가
upregulated_overlap_details <- extract_overlaps_with_details(upregulated_genes_list, pairwise_degs_list)
upregulated_overlap_details
write.csv(upregulated_overlap_details, file = file.path(output_dir, "Upregulated_Gene_Overlap_Details.csv"), row.names = FALSE)

# Downregulated 유전자 간 겹치는 유전자 확인 및 상세 정보 추가
downregulated_overlap_details <- extract_overlaps_with_details(downregulated_genes_list, pairwise_degs_list)
downregulated_overlap_details
write.csv(downregulated_overlap_details, file = file.path(output_dir, "Downregulated_Gene_Overlap_Details.csv"), row.names = FALSE)

# 출력
cat("Upregulated Gene Overlap Details:\n")
print(upregulated_overlap_details)

cat("\nDownregulated Gene Overlap Details:\n")
print(downregulated_overlap_details)

cat("\nOverlap details saved to:", output_dir, "\n")





# Upregulated 및 Downregulated DEG 데이터 추출
total_upregulated_genes <- unique(unlist(lapply(pairwise_degs_list, function(x) {
  rownames(x$upregulated_genes)
})))

total_downregulated_genes <- unique(unlist(lapply(pairwise_degs_list, function(x) {
  rownames(x$downregulated_genes)
})))

# Upregulated 유전자에 대한 히트맵 데이터 생성
upregulated_heatmap_data <- filtered_gene_exprs[total_upregulated_genes, , drop = FALSE]
write.csv(upregulated_heatmap_data, file = file.path(output_dir, "Upregulated_DEG_Heatmap_Data.csv"), row.names = TRUE)

# Downregulated 유전자에 대한 히트맵 데이터 생성
downregulated_heatmap_data <- filtered_gene_exprs[total_downregulated_genes, , drop = FALSE]
write.csv(downregulated_heatmap_data, file = file.path(output_dir, "Downregulated_DEG_Heatmap_Data.csv"), row.names = TRUE)

# 히트맵 생성 및 저장
library(pheatmap)

# Upregulated DEG 히트맵 생성
pheatmap(
  upregulated_heatmap_data,
  scale = "row",
  annotation_col = cluster_annotation,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Upregulated DEG Heatmap"
)
ggsave(file.path(output_dir, "Upregulated_DEG_Heatmap.png"), width = 10, height = 8)

# # Downregulated DEG 히트맵 생성
# pheatmap(
#   downregulated_heatmap_data,
#   scale = "row",
#   annotation_col = cluster_annotation,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   main = "Downregulated DEG Heatmap"
# )
# ggsave(file.path(output_dir, "Downregulated_DEG_Heatmap.png"), width = 10, height = 8)
# 
# cat("Heatmaps for Upregulated and Downregulated DEGs saved to:", output_dir, "\n")
# 






library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# GO term 분석 함수 정의
perform_go_analysis <- function(gene_list, output_dir, comp_name, ontology = "BP") {
  if (length(gene_list) == 0) {
    cat(sprintf("No genes available for GO analysis: %s\n", comp_name))
    return(NULL)
  }
  
  # GO term 분석 수행
  go_results <- enrichGO(
    gene          = gene_list,        # 유전자 리스트
    OrgDb         = org.Hs.eg.db,     # 인간 유전자 데이터베이스
    keyType       = "SYMBOL",         # 유전자 ID 형식
    ont           = ontology,         # "BP", "MF", "CC" 선택 가능
    pAdjustMethod = "BH",             # 다중 검정 방법 (Benjamini-Hochberg)
    pvalueCutoff  = 0.05,             # p-value 컷오프
    qvalueCutoff  = 0.2               # q-value 컷오프
  )
  
  if (is.null(go_results) || nrow(as.data.frame(go_results)) == 0) {
    cat(sprintf("No significant GO terms found for %s\n", comp_name))
    return(NULL)
  }
  
  # GO term 결과 저장
  go_results_df <- as.data.frame(go_results)
  output_file <- file.path(output_dir, paste0(comp_name, "_GO_", ontology, ".csv"))
  write.csv(go_results_df, file = output_file, row.names = FALSE)
  
  cat(sprintf("GO term results saved to: %s\n", output_file))
  return(go_results)
}

# GO term 분석 수행 및 저장
go_results_list_up <- list()
go_results_list_down <- list()

# Upregulated 분석
for (comp in names(upregulated_genes_list)) {
  gene_list <- upregulated_genes_list[[comp]]
  if (length(gene_list) > 0) {
    cat(sprintf("Performing GO term analysis for Upregulated: %s...\n", comp))
    
    # Biological Process (BP)
    go_results_bp <- perform_go_analysis(
      gene_list    = gene_list,
      output_dir   = output_dir,
      comp_name    = paste0(comp, "_Upregulated"),
      ontology     = "BP"
    )
    
    # Molecular Function (MF)
    go_results_mf <- perform_go_analysis(
      gene_list    = gene_list,
      output_dir   = output_dir,
      comp_name    = paste0(comp, "_Upregulated"),
      ontology     = "MF"
    )
    
    # Cellular Component (CC)
    go_results_cc <- perform_go_analysis(
      gene_list    = gene_list,
      output_dir   = output_dir,
      comp_name    = paste0(comp, "_Upregulated"),
      ontology     = "CC"
    )
    
    # 결과 저장
    go_results_list_up[[comp]] <- list(
      BP = go_results_bp,
      MF = go_results_mf,
      CC = go_results_cc
    )
  } else {
    cat(sprintf("No genes available for Upregulated GO term analysis for %s.\n", comp))
  }
}

# Downregulated 분석
for (comp in names(downregulated_genes_list)) {
  gene_list <- downregulated_genes_list[[comp]]
  if (length(gene_list) > 0) {
    cat(sprintf("Performing GO term analysis for Downregulated: %s...\n", comp))
    
    # Biological Process (BP)
    go_results_bp <- perform_go_analysis(
      gene_list    = gene_list,
      output_dir   = output_dir,
      comp_name    = paste0(comp, "_Downregulated"),
      ontology     = "BP"
    )
    
    # Molecular Function (MF)
    go_results_mf <- perform_go_analysis(
      gene_list    = gene_list,
      output_dir   = output_dir,
      comp_name    = paste0(comp, "_Downregulated"),
      ontology     = "MF"
    )
    
    # Cellular Component (CC)
    go_results_cc <- perform_go_analysis(
      gene_list    = gene_list,
      output_dir   = output_dir,
      comp_name    = paste0(comp, "_Downregulated"),
      ontology     = "CC"
    )
    
    # 결과 저장
    go_results_list_down[[comp]] <- list(
      BP = go_results_bp,
      MF = go_results_mf,
      CC = go_results_cc
    )
  } else {
    cat(sprintf("No genes available for Downregulated GO term analysis for %s.\n", comp))
  }
}

# 시각화
# Upregulated
for (comp in names(go_results_list_up)) {
  if (!is.null(go_results_list_up[[comp]]$BP)) {
    print(dotplot(go_results_list_up[[comp]]$BP, showCategory = 10) +
            ggtitle(paste("Top 10 GO Terms (BP) for Upregulated:", comp)))
  }
}

# Downregulated
for (comp in names(go_results_list_down)) {
  if (!is.null(go_results_list_down[[comp]]$BP)) {
    print(dotplot(go_results_list_down[[comp]]$BP, showCategory = 10) +
            ggtitle(paste("Top 10 GO Terms (BP) for Downregulated:", comp)))
  }
}

#HLA ONLY!!
hla_genes <- grep("HLA", rownames(gene_exprs), value = TRUE)
hla_gene_exprs <- gene_exprs[hla_genes, ]

# Highly Variable Genes (HVG) 추출
# 유전자별 분산 계산
hla_gene_variances <- apply(hla_gene_exprs, 1, var)

# 상위 변이를 가진 유전자 필터링 (상위 500개 기준)
num_top_hvg <- 16
hla_top_genes <- names(sort(hla_gene_variances, decreasing = TRUE))[1:num_top_hvg]
hla_gene_exprs <- hla_gene_exprs[hla_top_genes, ]

cat("Number of HLA genes selected as highly variable:", length(hla_top_genes), "\n")

# 1. PCA 수행
hla_filtered_pca <- prcomp(t(hla_gene_exprs), scale. = TRUE)

# 분산 비율 계산
explained_variance <- summary(hla_filtered_pca)$importance[2, ]  # 각 PC의 분산 비율
cumulative_variance <- cumsum(explained_variance)               # 누적 분산 비율

# PC 몇 개가 80% 이상의 분산을 설명하는지 확인
num_pcs_80 <- which(cumulative_variance >= 0.8)[1]
cat("Number of PCs explaining >= 80% variance:", num_pcs_80, "\n")

# 분산 비율 시각화
library(ggplot2)
variance_df <- data.frame(
  PC = 1:length(explained_variance),
  Variance = explained_variance,
  CumulativeVariance = cumulative_variance
)

variance_plot <- ggplot(variance_df, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = CumulativeVariance), color = "red", size = 1) +
  geom_point(aes(y = CumulativeVariance), color = "red", size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen", size = 0.8) +
  theme_minimal() +
  labs(
    title = "Explained Variance by Principal Components",
    x = "Principal Component (PC)",
    y = "Variance Explained"
  )

print(variance_plot)

# 2. 클러스터 수 결정 방법
library(cluster)

# 2.1 Elbow Method
wss <- sapply(1:10, function(k) {
  kmeans_result <- kmeans(t(hla_gene_exprs), centers = k, nstart = 25)
  kmeans_result$tot.withinss
})

elbow_df <- data.frame(
  Clusters = 1:10,
  WSS = wss
)

elbow_plot <- ggplot(elbow_df, aes(x = Clusters, y = WSS)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2, color = "red") +
  theme_minimal() +
  labs(
    title = "Elbow Method for Optimal Cluster Selection",
    x = "Number of Clusters",
    y = "Within-Cluster Sum of Squares (WSS)"
  )

print(elbow_plot)

# 2.2 Gap Statistic
set.seed(123)
gap_stat <- clusGap(
  t(hla_gene_exprs), 
  FUN = kmeans, 
  K.max = 10, 
  B = 50 # Number of bootstrap samples
)

# Gap Plot
library(factoextra)
gap_plot <- fviz_gap_stat(gap_stat)
print(gap_plot)

# Optimal cluster based on Gap Statistic
optimal_clusters_gap <- which.max(gap_stat$Tab[, "gap"])
cat("Optimal number of clusters (based on Gap Statistic):", optimal_clusters_gap, "\n")

# 2.3 Silhouette Score
silhouette_scores <- sapply(2:10, function(k) {
  kmeans_result <- kmeans(t(hla_gene_exprs), centers = k, nstart = 25)
  silhouette_values <- silhouette(kmeans_result$cluster, dist(t(hla_gene_exprs)))
  mean(silhouette_values[, 3])
})

silhouette_df <- data.frame(
  Clusters = 2:10,
  SilhouetteScore = silhouette_scores
)

silhouette_plot <- ggplot(silhouette_df, aes(x = Clusters, y = SilhouetteScore)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2, color = "red") +
  theme_minimal() +
  labs(
    title = "Silhouette Scores for Different Cluster Numbers",
    x = "Number of Clusters",
    y = "Average Silhouette Score"
  )

print(silhouette_plot)

# 3. 클러스터링 수행 (최적 클러스터 수 기반)
hla_num_clusters <- max(optimal_clusters_gap, which.max(silhouette_scores) + 1)
cat("Selected optimal number of clusters:", hla_num_clusters, "\n")

#hla_num_clusters <- 5

# PCA 기반 거리 행렬 생성
dist_matrix <- dist(hla_filtered_pca$x[, 1:num_pcs_80])
hla_hclust_clustering <- hclust(dist_matrix, method = "ward.D2")
hla_hclust_clusters <- cutree(hclust(dist_matrix, method = "ward.D2"), k = hla_num_clusters)
hla_hclust_clusters

# Consensus Clustering
library(ConsensusClusterPlus)
hla_consensus_results <- ConsensusClusterPlus(
  as.matrix(hla_gene_exprs), 
  maxK = 10,
  reps = 50,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 123,
  plot = "pdf",
  title = file.path(output_dir, "HLA_Consensus_Clustering")
)

hla_consensus_clusters <- hla_consensus_results[[hla_num_clusters]]$consensusClass

# NMF Clustering
library(NMF)
hla_nmf_results <- nmf(hla_gene_exprs, rank = hla_num_clusters, method = "brunet", seed = 123)
hla_nmf_clusters <- predict(hla_nmf_results)

# K-Means Clustering
set.seed(123)
hla_kmeans_results <- kmeans(t(hla_gene_exprs), centers = hla_num_clusters, nstart = 25)
hla_kmeans_clusters <- hla_kmeans_results$cluster

# PCA 결과 데이터프레임 생성
hla_filtered_pca_df <- data.frame(
  PC1 = hla_filtered_pca$x[, 1],
  PC2 = hla_filtered_pca$x[, 2],
  Samples = colnames(hla_gene_exprs),
  HierarchicalCluster = factor(hla_hclust_clusters, labels = paste0("Cluster", 1:hla_num_clusters)),
  ConsensusCluster = factor(hla_consensus_clusters),
  NMFCluster = factor(hla_nmf_clusters),
  KMeansCluster = factor(hla_kmeans_clusters)
)

# 4. PCA 플롯 생성
library(ggplot2)

# Hierarchical Clustering PCA Plot
hla_pca_plot_hierarchical <- ggplot(hla_filtered_pca_df, aes(x = PC1, y = PC2, color = HierarchicalCluster)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +
  theme_minimal() +
  labs(
    title = "PCA Plot with Hierarchical Clustering (HLA)",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  )

# Consensus Clustering PCA Plot
hla_pca_plot_consensus <- ggplot(hla_filtered_pca_df, aes(x = PC1, y = PC2, color = ConsensusCluster)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +
  theme_minimal() +
  labs(
    title = "PCA Plot with Consensus Clustering (HLA)",
    x = "PC1",
    y = "PC2",
    color = "Consensus Cluster"
  )

# NMF Clustering PCA Plot
hla_pca_plot_nmf <- ggplot(hla_filtered_pca_df, aes(x = PC1, y = PC2, color = NMFCluster)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +
  theme_minimal() +
  labs(
    title = "PCA Plot with NMF Clustering (HLA)",
    x = "PC1",
    y = "PC2",
    color = "NMF Cluster"
  )

# K-Means Clustering PCA Plot
hla_pca_plot_kmeans <- ggplot(hla_filtered_pca_df, aes(x = PC1, y = PC2, color = KMeansCluster)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Samples), hjust = 1.2, vjust = 1.2, size = 2, color = "black") +
  theme_minimal() +
  labs(
    title = "PCA Plot with K-Means Clustering (HLA)",
    x = "PC1",
    y = "PC2",
    color = "K-Means Cluster"
  )

# Display all PCA Plots
print(hla_pca_plot_hierarchical)
print(hla_pca_plot_consensus)
print(hla_pca_plot_nmf)
print(hla_pca_plot_kmeans)

# Save PCA plots
ggsave(
  filename = file.path(output_dir, "PCA_HLA_Hierarchical.png"),
  plot = hla_pca_plot_hierarchical,
  width = 10,
  height = 8
)
ggsave(
  filename = file.path(output_dir, "PCA_HLA_Consensus.png"),
  plot = hla_pca_plot_consensus,
  width = 10,
  height = 8
)
ggsave(
  filename = file.path(output_dir, "PCA_HLA_NMF.png"),
  plot = hla_pca_plot_nmf,
  width = 10,
  height = 8
)
ggsave(
  filename = file.path(output_dir, "PCA_HLA_KMeans.png"),
  plot = hla_pca_plot_kmeans,
  width = 10,
  height = 8
)

# Silhouette Score 계산
library(cluster)

# Silhouette Score를 계산할 거리 행렬
distance_matrix <- dist(hla_filtered_pca$x[, 1:num_pcs_80])

# Silhouette Score 계산 함수
calculate_silhouette <- function(clusters, distance_matrix) {
  silhouette_values <- silhouette(clusters, distance_matrix)
  mean(silhouette_values[, 3]) # 평균 Silhouette Score 반환
}

# Silhouette Scores 계산
silhouette_scores <- list(
  Hierarchical = calculate_silhouette(hla_hclust_clusters, distance_matrix),
  Consensus = calculate_silhouette(hla_consensus_clusters, distance_matrix),
  NMF = calculate_silhouette(as.numeric(hla_nmf_clusters), distance_matrix),
  KMeans = calculate_silhouette(hla_kmeans_clusters, distance_matrix)
)

# Silhouette Score 출력
print(silhouette_scores)






# 1. 클러스터링용 데이터 생성 (HLA 데이터 제외)
clustering_data <- as.matrix(upregulated_heatmap_data)  # Upregulated 데이터만 사용

# 2. 클러스터링 결과 계산 (사전 클러스터링)
row_dist <- dist(clustering_data)                  # 행 거리 계산
row_clustering <- hclust(row_dist, method = "ward.D2")  # 행 클러스터링

col_dist <- dist(t(clustering_data))               # 열 거리 계산
col_clustering <- hclust(col_dist, method = "ward.D2")  # 열 클러스터링

# 3. HLA 데이터 반복으로 시각적 추가를 위한 데이터 생성
hla_repeated <- hla_gene_exprs[rep(rownames(hla_gene_exprs), each = 30), ]
rownames(hla_repeated) <- paste0(rownames(hla_repeated), "_rep")  # 이름 변경하여 구분

# 4. 최종 히트맵 데이터 병합
final_heatmap_data <- rbind(clustering_data, hla_repeated)

# 5. 행 이름 설정 (HLA 유전자만 표시)
row_labels <- rownames(final_heatmap_data)
row_labels[!(rownames(final_heatmap_data) %in% rownames(hla_gene_exprs))] <- ""


# 5. 행 이름 설정 (HLA 유전자만 표시)
#rep_hla_gene <- paste0(rownames(hla_gene_exprs), "_rep")
row_labels <- rownames(final_heatmap_data)
row_labels[!(rownames(final_heatmap_data) %in% rownames(hla_repeated))] <- ""  # HLA 유전자만 이름 표시


# 전체 데이터에서 30개 간격으로 그룹화하여 처리
for (i in seq(1, length(row_labels), by = 30)) {
  group_indices <- i:(i + 29)  # 30개씩 그룹화
  group_indices <- group_indices[group_indices <= length(row_labels)]  # 범위를 초과하지 않도록 조정
  
  if (length(group_indices) >= 15) {  # 15번째 요소가 존재하는 경우
    row_labels[group_indices[-15]] <- ""  # 15번째를 제외한 나머지를 빈 문자열로 설정
  } else {
    row_labels[group_indices] <- ""  # 30개 미만인 그룹은 모두 빈 문자열로 처리
  }
}

# 6. . 기준으로 split하고 0번 인덱스만 남기기
row_labels <- sapply(row_labels, function(label) {
  if (label == "") {
    return("")  # 빈 문자열은 그대로 유지
  } else {
    return(strsplit(label, "\\.")[[1]][1])  # . 기준으로 split 후 첫 번째 인덱스만 반환
  }
})

# 결과 확인
tail(row_labels, 50)  # 결과 확인



final_heatmap_data  <- as.matrix(final_heatmap_data)

# 6. 클러스터링 결과 적용 및 HLA 데이터 제외 히트맵 생성
pheatmap(
  final_heatmap_data,
  cluster_rows = FALSE,                # HLA 데이터가 클러스터링에 영향을 미치지 않도록 설정
  cluster_cols = col_clustering,       # 열 클러스터링 결과 적용
  annotation_col = cluster_annotation, # 클러스터 주석 추가
#  annotation_colors = annotation_colors, # 주석 색상 추가
  labels_row = row_labels,             # HLA 유전자만 이름 표시
  show_colnames = TRUE,                # 샘플 이름 표시
  #  gaps_row = nrow(clustering_data),    # HLA와 Upregulated 데이터 사이에 간격 추가
  scale = "row",                       # 행별 표준화
  color = colorRampPalette(c("blue", "white", "red"))(50), # 색상
  main = "Heatmap with HLA Genes Always at Bottom"
)


# cell_line_data 정의
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

# 히트맵 데이터 열 이름 변경
# colnames(final_heatmap_data)를 cell_line_data의 CellLine으로 매핑
colnames(final_heatmap_data) <- cell_line_data$CellLine[
  match(colnames(final_heatmap_data), cell_line_data$Sample)
]

# 히트맵 생성
pheatmap(
  final_heatmap_data,
  cluster_rows = FALSE,
  cluster_cols = col_clustering,
  annotation_col = cluster_annotation,  # (Optional) 다른 주석 추가
#  annotation_colors = annotation_colors,
  labels_row = row_labels,
  show_colnames = TRUE,  # 이제 Cell Line 이름이 열 이름으로 표시됨
  scale = "row",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap with CellLine Names"
)


# Step 2: 히트맵 생성 시 Annotation 업데이트
# 클러스터 annotation 데이터프레임 생성
cluster_annotation <- data.frame(
  Sample = colnames(filtered_gene_exprs),
  Cluster = filtered_pca_df$ConsensusCluster
)
rownames(cluster_annotation) <- cluster_annotation$Sample




# # PDF 출력 디바이스 설정
# pdf("/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/DEG_result/Array_DEG/AML_heatmap_with_cellline.pdf", width = 10, height = 8)  # PDF 파일 이름과 크기 지정
# 
# # 히트맵 생성
# pheatmap(
#   final_heatmap_data,
#   cluster_rows = FALSE,
#   cluster_cols = col_clustering,
#   annotation_col = cluster_annotation,  # (Optional) 다른 주석 추가
#   annotation_colors = annotation_colors,
#   labels_row = row_labels,
#   show_colnames = TRUE,  # 열 이름 표시
#   scale = "row",
#   color = colorRampPalette(c("blue", "white", "red"))(50),
#   main = "Heatmap with CellLine Names"
# )
# 
# # PDF 저장 종료
# dev.off()





library(pheatmap)

# 1. 데이터 클러스터링 및 정렬 함수 정의
perform_clustering <- function(data, cluster_indices) {
  row_dist <- dist(data)  # 거리 계산
  row_clustering <- hclust(row_dist, method = "ward.D2")  # 클러스터링 수행
  cluster_order <- order.dendrogram(as.dendrogram(row_clustering))  # 순서 추출
  return(cluster_indices[cluster_order])
}

# 2. Gap 생성 함수 정의
generate_gap <- function(ncol, gap_length = 20) {
  gap_row <- matrix(NA, nrow = gap_length, ncol = ncol)
  colnames(gap_row) <- colnames(final_heatmap_data)
  rownames(gap_row) <- paste0("gap_", seq_len(gap_length))
  return(gap_row)
}

# 3. 라벨 생성 함수 정의
create_row_labels <- function(non_hla_labels, hla_labels, gap_length = 20) {
  # Non-HLA 라벨 모두 빈 문자열로 설정
  non_hla_labels <- rep("", length(non_hla_labels))
  
  # HLA 라벨 조건 처리
  hla_labels_processed <- rep("", length(hla_labels))
  for (i in seq(1, length(hla_labels), by = 30)) {
    group_indices <- i:(i + 29)
    group_indices <- group_indices[group_indices <= length(hla_labels)]
    
    if (length(group_indices) >= 15) {
      hla_labels_processed[group_indices[15]] <- hla_labels[group_indices[15]]
    }
    if (length(group_indices) >= 45) {
      hla_labels_processed[group_indices[45]] <- hla_labels[group_indices[45]]
    }
  }
  
  # Gap 라벨 생성
  gap_labels <- rep("", gap_length)
  
  # 라벨 결합
  return(c(non_hla_labels, gap_labels, hla_labels_processed))
}

# 4. Non-HLA와 HLA 데이터 분리
non_hla_data <- final_heatmap_data[!(rownames(final_heatmap_data) %in% rownames(hla_repeated)), ]
hla_data <- final_heatmap_data[rownames(final_heatmap_data) %in% rownames(hla_repeated), ]

# 5. 클러스터링 수행
non_hla_row_order <- perform_clustering(non_hla_data, rownames(non_hla_data))
hla_row_order <- perform_clustering(hla_data, rownames(hla_data))

# 6. 데이터 정렬 및 Gap 결합
non_hla_data <- final_heatmap_data[non_hla_row_order, ]
hla_data <- final_heatmap_data[hla_row_order, ]
gap_row <- generate_gap(ncol(final_heatmap_data))
final_ordered_with_gap <- rbind(non_hla_data, gap_row, hla_data)

# 7. 행 라벨 생성
non_hla_labels <- rownames(non_hla_data)
hla_labels <- sapply(rownames(hla_data), function(label) strsplit(label, "\\.")[[1]][1])
row_labels_with_gap <- create_row_labels(non_hla_labels, hla_labels, nrow(gap_row))

# 8. Annotation 정렬
final_annotation <- cluster_annotation[colnames(final_ordered_with_gap), , drop = FALSE]

# 9. CellLine 데이터 매핑
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

colnames(final_ordered_with_gap) <- cell_line_data$CellLine[
  match(colnames(final_ordered_with_gap), cell_line_data$Sample)
]

# 10. 히트맵 생성
pheatmap(
  final_ordered_with_gap,
  cluster_rows = FALSE,                     # 행 클러스터링 결과 고정
  cluster_cols = col_clustering,            # 열 클러스터링 추가
  annotation_col = final_annotation,        # 수정된 주석 데이터
  labels_row = row_labels_with_gap,         # Gap 포함 라벨
  show_colnames = TRUE,                     # Cell Line 이름 표시
  scale = "row",                            # 행 기준 표준화
  color = colorRampPalette(c("blue", "white", "red"))(50),  # 색상 조합
  na_col = "white",                         # Gap을 하얀색으로 표시
  main = "Heatmap with CellLine Names and Conditional HLA Labels"
)


# PDF 출력 디바이스 설정
pdf("/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/DEG_result/Array_DEG/AML_heatmap_with_cellline.pdf", width = 10, height = 8)  # PDF 파일 이름과 크기 지정

# 히트맵 생성
pheatmap(
  final_ordered_with_gap,
  cluster_rows = FALSE,                     # 행 클러스터링 결과 고정
  cluster_cols = col_clustering,            # 열 클러스터링 추가
  annotation_col = final_annotation,        # 수정된 주석 데이터
  labels_row = row_labels_with_gap,         # Gap 포함 라벨
  show_colnames = TRUE,                     # 열 이름 표시 (Cell Line 이름)
  scale = "row",                            # 행 기준 표준화
  color = colorRampPalette(c("blue", "white", "red"))(50),  # 색상 조합
  main = "Heatmap with CellLine Names"
)

# PDF 저장 종료
dev.off()
    
