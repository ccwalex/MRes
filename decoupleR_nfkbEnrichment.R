#supplementary data on pathway enrichment with decoupleR
library(decoupleR)
library(pheatmap)
# Extract the normalized log-transformed counts
mat <- as.matrix(EC[['RNA']]$scale.data)

#reference
net <- decoupleR::get_progeny(organism = 'human', top = 500)
# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
EC[['pathwaysmlm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames(var = 'source') %>%
  Seurat::CreateAssayObject(.)


# Change assay
Seurat::DefaultAssay(object = EC) <- "pathwaysmlm"

# Scale the data
EC <- Seurat::ScaleData(EC)
EC@assays$pathwaysmlm@data <- EC@assays$pathwaysmlm@scale.data

p1 <- Seurat::DimPlot(EC, 
                      reduction = "umap.EC", 
                      group.by = 'disease__ontology_label',
                      #label = TRUE, 
                      pt.size = 0.5) + 
  #Seurat::NoLegend() + 
  ggplot2::ggtitle('Cell types')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p2 <- Seurat::FeaturePlot(EC, features = c("NFkB"), reduction = 'umap.EC') + 
  ggplot2::scale_colour_gradient2(low = colors[1], mid = 'white', high = colors[2]) +
  ggplot2::ggtitle('NFkb activity')

p <- p1 | p2
p

p3 <- Seurat::FeaturePlot(EC, features = c("NFkB"), reduction = 'umap.EC') + 
  ggplot2::scale_colour_gradient2(low = colors[1], mid = 'white', high = colors[2]) +
  ggplot2::ggtitle('NFkb activity')
p3

EC_id = EC$disease__ontology_label
pathway_matrix = as.data.frame(acts)
pathway_matrix = pathway_matrix[pathway_matrix$source == 'NFkB',]
pathway_matrix = cbind(pathway_matrix, EC_id)
pathway_matrix = pathway_matrix[pathway_matrix$p_value <= 0.05,]


NFkb = pathway_matrix %>%
  t_test(score ~ EC_id, detailed = TRUE) %>%
  add_significance()
NFkb = add_xy_position(NFkb, x = 'condition')



NFkbp = ggboxplot(
  pathway_matrix, x = 'EC_id', y = 'score',
  ylab = 'NFkB pathway enrichment',
  xlab = 'condition',
  title = 'NFkB pathway enrichment'
  ,color = 'EC_id'
  
)

NFkbp +
  stat_pvalue_manual(NFkb, tip.length = 0) +
  labs(subtitle = get_test_label(NFkb, detailed = TRUE))+
  theme(text = element_text(size = 16))