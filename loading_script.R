#script for scRNA lung autopsy, university of columbia dataset
#https://www.nature.com/articles/s41586-021-03569-1#data-availability

packages = c('dplyr', 'parallel', 'Seurat', 'SeuratObject', 'ggplot2', 'patchwork', 'ggrepel', 'SingleR', 'celldex')
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
lungtry[['RNA']] = split(lungtry[['RNA']], f = lung_try$orig.ident)
lungtry[['RNA']] = split(lung_try[['RNA']], f = lung_try$orig.ident)
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