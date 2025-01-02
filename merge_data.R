# 필요한 패키지 로드
library(dplyr)

# 파일 경로 설정
cluster1_input_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Cluster_Cluster1/"
cluster2_input_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Cluster_Cluster2/"
output_file <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Combined_Gene_Expression_with_Cluster.csv"

# 디렉토리 내 Cluster1 및 Cluster2 CSV 파일 가져오기
cluster1_files <- list.files(cluster1_input_dir, pattern = "*_Cluster_Cluster1.csv", full.names = TRUE)
cluster2_files <- list.files(cluster2_input_dir, pattern = "*_Cluster_Cluster2.csv", full.names = TRUE)

# Cluster1과 Cluster2 파일 처리 및 병합
process_files <- function(files, cluster_label) {
  lapply(files, function(file) {
    # CSV 파일 읽기
    data <- read.csv(file)
    
    # Sample과 CellLine 결합 및 Cluster1_2 열 추가
    data <- data %>%
      mutate(Sample = paste(Sample, CellLine, sep = "_"),
             MHC_gene_Cluster1_2 = cluster_label) %>%
      select(MHC_gene_Cluster1_2, Sample, Gene, Expression)
    
    return(data)
  }) %>%
    bind_rows()
}

# Cluster1 및 Cluster2 데이터 병합
cluster1_data <- process_files(cluster1_files, "Cluster1")
cluster2_data <- process_files(cluster2_files, "Cluster2")
combined_data <- bind_rows(cluster1_data, cluster2_data)

# Gene을 열로 변환
reshaped_data <- combined_data %>%
  pivot_wider(names_from = Gene, values_from = Expression)

# 결과 저장
write.csv(reshaped_data, output_file, row.names = FALSE)

cat("Combined gene expression data with Cluster1_2 column saved to:", output_file, "\n")













# 필요한 패키지 로드
library(dplyr)

# 파일 경로 설정
# 첫 번째 그룹 디렉토리
cluster1_input_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Cluster_Cluster1/"
cluster2_input_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Cluster_Cluster2/"

# 두 번째 그룹 디렉토리
clustered_cd36_dir1 <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Clustered_CD36_Expression/Cluster_Cluster1/"
clustered_cd36_dir2 <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Clustered_CD36_Expression/Cluster_Cluster2/"
clustered_cd36_dir3 <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Clustered_CD36_Expression/Cluster_Cluster3/"

output_file <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE98578/DEG_results/Combined_Gene_Expression_with_All_Clusters.csv"

# 디렉토리 내 CSV 파일 가져오기
cluster1_files <- list.files(cluster1_input_dir, pattern = "*_Cluster_Cluster1.csv", full.names = TRUE)
cluster2_files <- list.files(cluster2_input_dir, pattern = "*_Cluster_Cluster2.csv", full.names = TRUE)
clustered_cd36_files1 <- list.files(clustered_cd36_dir1, pattern = "*_Cluster_Cluster1.csv", full.names = TRUE)
clustered_cd36_files2 <- list.files(clustered_cd36_dir2, pattern = "*_Cluster_Cluster2.csv", full.names = TRUE)
clustered_cd36_files3 <- list.files(clustered_cd36_dir3, pattern = "*_Cluster_Cluster3.csv", full.names = TRUE)

# 파일 처리 및 병합 함수 정의
process_files <- function(files, cluster_label, column_name) {
  lapply(files, function(file) {
    # CSV 파일 읽기
    data <- read.csv(file)
    
    # CellLine 열이 없는 경우 기본값 추가
    if (!"CellLine" %in% colnames(data)) {
      data <- data %>% mutate(CellLine = "NA")
    }
    
    # Sample과 CellLine 결합 및 Cluster 열 추가
    data <- data %>%
      mutate(Sample = paste(Sample, CellLine, sep = "_"),
             {{ column_name }} := cluster_label) %>%
      select({{ column_name }}, Sample, Gene, Expression)
    
    return(data)
  }) %>%
    bind_rows()
}

# 첫 번째 그룹 데이터 병합
cluster1_data <- process_files(cluster1_files, "Cluster1", MHC_Cluster_12)
cluster2_data <- process_files(cluster2_files, "Cluster2", MHC_Cluster_12)

# 두 번째 그룹 데이터 병합
clustered_cd36_data1 <- process_files(clustered_cd36_files1, "Cluster1", Top5000_gene_cluster_123_For_CD36)
clustered_cd36_data2 <- process_files(clustered_cd36_files2, "Cluster2", Top5000_gene_cluster_123_For_CD36)
clustered_cd36_data3 <- process_files(clustered_cd36_files3, "Cluster3", Top5000_gene_cluster_123_For_CD36)

# 모든 데이터 병합
combined_data <- bind_rows(cluster1_data, cluster2_data, clustered_cd36_data1, clustered_cd36_data2, clustered_cd36_data3)

# Gene을 열로 변환, NA로 채우기
reshaped_data <- combined_data %>%
  pivot_wider(names_from = Gene, values_from = Expression, values_fill = NA) %>%
  relocate(MHC_Cluster_12, Top5000_gene_cluster_123_For_CD36, Sample)

# 결과 저장
write.csv(reshaped_data, output_file, row.names = FALSE)

cat("Combined gene expression data with MHC_Cluster_12 and Top5000_gene_cluster_123_For_CD36 columns saved to:", output_file, "\n")

print(reshaped_data)

