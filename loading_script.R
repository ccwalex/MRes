#script for scRNA lung autopsy, university of columbia dataset
#https://www.nature.com/articles/s41586-021-03569-1#data-availability

packages = c('dplyr', 'parallel', 'Seurat', 'SeuratObject', 'ggplot2', 'patchwork', 'ggrepel', 'SingleR', 'celldex', 'speckle', 'tidyverse')
lapply(packages, require, character.only = TRUE)


lung = ReadMtx(
  cells = 'lung_cellNames.csv',
  mtx = 'gene_sorted-lung_expression_data.mtx.gz',
  features = 'lung_geneNames_upload.csv',
  mtx.transpose = FALSE,
  feature.column = 1
  
)
lungmeta = read.delim2('lung_metaData.txt', header = TRUE)
lungmetaraw = lungmeta[-1,]

rownames(lungmetaraw) = lungmetaraw$NAME

lung_try = CreateSeuratObject(counts = lung)

lungobject = CreateSeuratObject(counts = lung
                                , meta.data = lungmetaraw
)
lungobject[['percent.mt']] = PercentageFeatureSet(lungobject, pattern = '^MT-')
lungobject = PercentageFeatureSet(lungobject, pattern = '^RPS', col.name =  'percent.rps')
lungobject = PercentageFeatureSet(lungobject, pattern = '^RPL', col.name =  'percent.rpl')
lungobject$percent.rp = lungobject$percent.rps + lungobject$percent.rpl
#all filtered, doublets also removed in data


lungobject = NormalizeData(lungobject)
lungobject = FindVariableFeatures(lungobject)
lungobject = ScaleData(lungobject)
lungobject = RunPCA(lungobject)
lungobject = JackStraw(lungobject, num.replicate = 100)
lungobject = ScoreJackStraw(lungobject, reduction = 'pca', dims = 1:20)
lungobject = FindNeighbors(lungobject, dims = 1:30)
lungobject = FindClusters(lungobject)
lungobject = RunUMAP(lungobject, dims = 1:30)

#removing broken/corrupted cells data
lungobject = subset(lungobject, subset = NAME != 'NA')

lungobject[['RNA']] = split(lungobject[['RNA']], f = lungobject@meta.data$donor_id)

#use this
lung_try = CreateSeuratObject(counts = lung)
lung_try[['RNA']] = split(lung_try[['RNA']], f = lung_try$orig.ident)
lung_try[['RNA']] = split(lung_try[['RNA']], f = lung_try$orig.ident)
lung_try[['RNA']] = split(lung_try[['RNA']], f = lung_try$orig.ident)
lung_try = AddMetaData(lung_try, metadata = lungmetaraw)

lung_try = NormalizeData(lung_try)
lung_try = FindVariableFeatures(lung_try)
lung_try = ScaleData(lung_try)
lung_try = RunPCA(lung_try)
lung_try = JackStraw(lung_try, num.replicate = 100)
lung_try = ScoreJackStraw(lung_try, reduction = 'pca', dims = 1:20)
lung_try = FindNeighbors(lung_try, dims = 1:30)
lung_try = FindClusters(lung_try)
lung_try = RunUMAP(lung_try, dims = 1:30)

options(future.globals.maxSize= 2017897210)

lung_integrate = IntegrateLayers(object = lungobject
                                 , method = RPCAIntegration
                                 , orig.reduction = 'pca'
                                 , new.reduction = 'integrated.rpca'
                                 #, normalization.method = "SCT"
                                 , verbose = TRUE)

set.seed(388388)
lung_integrate[["RNA"]] <- JoinLayers(lung_integrate[["RNA"]])
lung_integrate <- FindNeighbors(lung_integrate, reduction = "integrated.rpca", dims = 1:30)
lung_integrate <- FindClusters(lung_integrate, resolution = 1)
lung_integrate = RunUMAP(lung_integrate, reduction = "integrated.rpca", dims = 1:30)
DimPlot(lung_integrate, reduction = "umap", group.by = c("disease__ontology_label", "cell_type_main"))


EC = subset(lungobject, cells = WhichCells(lungobject, expression = cell_type_main == 'Endothelial cells'))
#EC = NormalizeData(EC)
EC = FindVariableFeatures(EC)
ECgenes = rownames(EC)
EC = ScaleData(EC, features = ECgenes)
EC = RunPCA(EC, features = VariableFeatures(object = EC))
EC = FindNeighbors(EC, dims = 1:25)
EC = FindClusters(EC, resolution = 2, cluster.name='EC.clusters')
EC = RunUMAP(EC, dims = 1:25, reduction.name = "umap.EC")

DimPlot(EC, reduction = 'umap.EC')
DimPlot(EC, group.by = 'disease__ontology_label', reduction = 'umap.EC')
DimPlot(EC, group.by = 'donor_id', reduction = 'umap.EC')
VariableFeaturePlot(EC)
RidgePlot(EC, features = c('CD34', 'VWF'), group.by = 'disease__ontology_label')
FeaturePlot(EC, features = c('CD34', 'VWF'), reduction = 'umap.EC', split.by = 'disease__ontology_label')
FeatureScatter(EC, 'CD34', 'VWF', split.by = 'disease__ontology_label')
RidgePlot(EC, features = c('CD34', 'VWF'), group.by = 'EC.clusters')

EC = JoinLayers(EC)
ECvar = FindAllMarkers(EC)

prop = propeller(clusters=EC$EC.clusters, sample=EC$donor_id, group=EC$disease__ontology_label)
propeller(clusters=lung_integrate$seurat_clusters, sample=lung_integrate$donor_id, group=lung_integrate$disease__ontology_label)

normal = subset(EC, cells = WhichCells(EC, expression = disease__ontology_label =='normal'))
disease = subset(EC, cells = WhichCells(EC, expression = disease__ontology_label =='COVID-19'))

FeaturePlot(lung_integrate, features = c('CD34', 'VWF'))
FeaturePlot(EC, features = c('CD34', 'VWF'))

RidgePlot(lung_integrate, features = c('CD34', 'VWF'))

#s.genes = cc.genes.updated.2019$s.genes
#g2m.genes = cc.genes.updated.2019$g2m.genes
#EC = CellCycleScoring(EC, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

DimPlot(EC, reduction = 'umap.EC', group.by = 'Phase', split.by = 'disease__ontology_label')


DimPlot(EC, reduction = 'umap.EC', split.by = 'disease__ontology_label')


Idents(EC_testdiff) = 'disease__ontology_label'
ECdiff = FindMarkers(EC_testdiff, ident.1 = 'COVID-19', ident.2 = 'normal')
write.csv(ECdiff, 'EC_conditionDE.csv')

gene_list = c('ICAM1', 'IL1B', 'IL6', 'VCAM1')
FeaturePlot(EC, features = gene_list, split.by = 'disease__ontology_label')
FeaturePlot(lung_integrate, features = gene_list, split.by = 'disease__ontology_label')
normal_id = c('data.C51ctr', 'data.C52ctr', 'data.C53ctr', 'data.C54ctr', 'data.C55ctr', 'data.C56ctr', 'data.C57ctr')
normalexp = list()
#define functions for gene t test
getcount = function(mtx, outputlist){
  expression = as.data.frame(rowMeans(mtx))
  outputlist = append(outputlist, list(expression[gene_list, ]))
  return(outputlist)
  
}
for (k in normal_id){
  normalexp = getcount(LayerData(normal, assay = 'RNA', layer = k), normalexp)
}

normal_genes = as.data.frame(do.call(rbind,normalexp))


#disease groups
dis_id = c('data.L01cov', 'data.L03cov', 'data.L05cov', 'data.L06cov', 'data.L07cov', 'data.L08cov', 'data.L09cov',
           'data.L10cov', 'data.L11cov', 'data.L12cov', 'data.L13cov', 'data.L15cov', 'data.L16cov', 'data.L17cov', 'data.L18cov', 'data.L19cov',
           'data.L21cov', 'data.L22cov')

disexp = list()
for (k in dis_id){
  disexp = getcount(LayerData(disease, assay = 'RNA', layer = k), disexp)
}

disease_genes = as.data.frame(do.call(rbind,disexp))

colnames(normal_genes) = gene_list
colnames(disease_genes) = gene_list

normal_genes
disease_genes

gene_list = list(gene_list)
#all genes log-transformed and scaled
#scaled expression from seurat object

disease_genes['label'] = 'disease'
normal_genes['label'] = 'normal'

merged_gene = bind_rows(normal_genes, disease_genes)
library(rstatix)
library(ggpubr)
IC = merged_gene %>%
  t_test(ICAM1 ~ label, detailed = TRUE) %>%
  add_significance()

ICp = ggboxplot(
  merged_gene, x = 'label', y = 'ICAM1',
  ylab = 'expression',
  xlab = 'condition',
  title = 'ICAM1'
  ,color = 'label'
)
IC = add_xy_position(IC, x = 'label')
ICp +
  stat_pvalue_manual(IC, tip.length = 0) +
  labs(subtitle = get_test_label(IC, detailed = TRUE)) +
  theme(text = element_text(size = 16))

IL1 = merged_gene %>%
  t_test(IL1B ~ label, detailed = TRUE) %>%
  add_significance()

IL1p = ggboxplot(
  merged_gene, x = 'label', y = 'IL1B',
  ylab = 'expression',
  xlab = 'condition',
  title = 'IL1B'
  ,color = 'label'
)
IL1 = add_xy_position(IL1, x = 'label')
IL1p +
  stat_pvalue_manual(IL1, tip.length = 0) +
  labs(subtitle = get_test_label(IL1, detailed = TRUE))+
  theme(text = element_text(size = 16))

VC = merged_gene %>%
  t_test(VCAM1 ~ label, detailed = TRUE) %>%
  add_significance()

VCp = ggboxplot(
  merged_gene, x = 'label', y = 'VCAM1',
  ylab = 'expression',
  xlab = 'condition',
  title = 'VCAM1'
  ,color = 'label'
)
VC = add_xy_position(NL, x = 'label')
VCp +
  stat_pvalue_manual(NL, tip.length = 0) +
  labs(subtitle = get_test_label(NL, detailed = TRUE))+
  theme(text = element_text(size = 16))

IL6 = merged_gene %>%
  t_test(IL6 ~ label, detailed = TRUE) %>%
  add_significance()
IL6 = add_xy_position(IL6, x = 'label')

IL6p = ggboxplot(
  merged_gene, x = 'label', y = 'IL6',
  ylab = 'expression',
  xlab = 'condition',
  title = 'IL6'
  ,color = 'label'
  
)
IL6p +
  stat_pvalue_manual(IL6, tip.length = 0) +
  labs(subtitle = get_test_label(IL6, detailed = TRUE))+
  theme(text = element_text(size = 16))

FeaturePlot(lung_integrate, features = c('ICAM1', 'IL1B', 'IL6', 'VCAM1'), reduction = 'umap')
FeaturePlot(EC, features = c('ICAM1', 'IL1B', 'IL6', 'VCAM1'), reduction = 'umap.EC')