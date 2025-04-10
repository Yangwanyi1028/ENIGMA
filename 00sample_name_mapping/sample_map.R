setwd('/home/wanyiyang/ENIGMA_STUDY_DATA/00sample_name_mapping/')
# 读取数据
mapping_df <- readxl::read_excel('../source_data/ID_match_20241018HK.xlsx')
mapping_df <- na.omit(mapping_df)  # 去除NA值
path <- read.csv('../source_data/humann_pathabundance_relab_AUS_HK_KM_479_unstratified.csv')
sp <- readxl::read_excel('../source_data/merged_species_abundance_mp4_20240614_605samples.xlsx')
gene <- read.csv('../source_data/ARG68_gene_abundance_20250224.csv',row.names = 1)
metadata <- readxl::read_excel('../source_data/ENIGMA_metadata.xlsx')

# 获取path的列名（排除第一列）
path_id <- colnames(path)[-1]
gene_id <- colnames(gene)[-c(1,2)]
meta_id <- metadata$study_id 
# 根据mapping_df映射列名
new_colnames_path <- ifelse(
  path_id %in% mapping_df$another_name,  # 检查是否在映射规则中
  mapping_df$study_id[match(path_id, mapping_df$another_name)],  # 映射
  path_id  # 如果未匹配，保留原始列名
)

# 根据mapping_df映射列名
new_colnames_gene <- ifelse(
  gene_id %in% mapping_df$another_name,  # 检查是否在映射规则中
  mapping_df$study_id[match(gene_id, mapping_df$another_name)],  # 映射
  gene_id  # 如果未匹配，保留原始列名
)

new_colnames_meta <- ifelse(
  meta_id %in% mapping_df$another_name,  # 检查是否在映射规则中
  mapping_df$study_id[match(meta_id, mapping_df$another_name)],  # 映射
  meta_id  # 如果未匹配，保留原始列名
)

# 更新path的列名
colnames(path)[-1] <- new_colnames_path
colnames(gene)[-c(1,2)] <- new_colnames_gene
# metadata$study_id <- new_colnames_meta
metadata <- metadata[metadata$study_id %in% new_colnames_path,]
# 检查sp$clade_name与映射后的列名的交集
mapped_list_path <- colnames(path)[-1]
intersect_result <- intersect(sp$clade_name, mapped_list_path)
print(intersect_result)

mapped_list_gene <- colnames(gene)[-1]
intersect_result <- intersect(sp$clade_name, mapped_list_gene)
print(intersect_result)

write.csv(data.frame(clade_name=intersect_result),'sample_479.csv')
write.csv(path,'../source_data/new_humann_pathabundance_relab_AUS_HK_KM_479_unstratified.csv',
          row.names = F)
write.csv(gene,'../source_data/new_ARG68_gene_abundance_20250224.csv',
          row.names = F)
write.csv(metadata,'../source_data/new_meta_data.csv',
          row.names = F)
