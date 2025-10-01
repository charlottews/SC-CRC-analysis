required_packages <- c(
  "Seurat", "dplyr", "patchwork", "sctransform", "ggpubr", "celldex", "ggplot2",
  "reshape2", "clustree", "data.table", "hash", "RColorBrewer", "tidyverse",
  "DoubletFinder", "ggsci", "ggrepel", "Matrix"
)


for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

##data
setwd()
load()



colnames(immune.combinedReso1@meta.data)
pdf("TSNE.data.HumanLiverSeurat.res0.8.pdf",width=6,height=5)
DimPlot(HumanLiverSeurat, reduction = "tsne",label = "T")+scale_color_igv()
dev.off()
immune.combinedReso1 <- HumanLiverSeurat

cell_markers <- list(
  Hep_5 = c("CYP3A7", "CYP2A7", "CYP2A6"),
  Hep_3 = c("SCD", "HMGCS1", "ACSS2", "TM7SF2"),
  Hep_2 = c("SEC16B", "SLBP", "RND3"),
  Hep_1 = c("PCK1", "BCHE"),
  Hep_6 = c("G6PC", "GHR"),
  Hep_4 = c("ALDH6A1", "RPP25L", "HSD11B1"),
  Periportal_Hep = c("HAMP", "GHR"),
  Pericentral_Hep = c("HPR", "GSTA2", "AKR1C1"),
  LSEC_Central = c("MASP2", "MGP", "SPARCL1"),
  LSEC_Portal = c("TM4SF1", "CLEC14A"),
  LSEC_General = c("CCL14", "CLEC1B", "FCN2"),
  Portal_endothelial = c("S100A13", "RAMP3"),
  Cholangiocyte = c("INMT", "DNASE1L3", "LIFR", "KRT7", "KRT19", "SOX9", "EPCAM"),
  Stellate_cells = c("ACTA2", "COL1A1", "RBP1"),
  Macrophage = c("S100A8", "LYZ", "S100A9", "HLA-DPB1", "CD5L", "MARCO", "VSIG4"),
  T_cell = c("CD2", "CD3D", "TRAC", "GZMK", "GNLY", "PTGDS", "GZMB", "TRDC"),
  Cycling_cells = c("STMN1", "HMGCB", "TYMS"),
  CD3_ab_T_cells = c("CD2", "CD3D"),
  CD4_T_cells = c("CD4"),
  B_cell = c("CD19")
)
combined_markers <- unlist(cell_markers, use.names = FALSE)
DefaultAssay(immune.combinedReso1) <- "RNA"


# 动态计算尺寸参数
for(cell_type in names(cell_markers)) {
markers <- cell_markers[[cell_type]]
  num_genes <- length(markers)
  base_width_per_gene <- 4 # 每个基因分配的基础宽度（英寸）
  min_width <- 12             # 最小PDF宽度
  max_width <- 30             # 最大PDF宽度（A4纸宽度为11.7英寸）
  base_height <- 4            # 基础高度
   plot_width <- pmin(max_width,
                    pmax(min_width,
                        (num_genes * base_width_per_gene) + 0.5))
# 高度按基因数量对数缩放（保证可读性）
  plot_height <- base_height + log(num_genes + 1) * 1.2


  plot_name <- paste0("VlnPlot_", cell_type, ".pdf")
  pdf(plot_name, width = plot_width, height = plot_height)
  print(VlnPlot(immune.combinedReso1, features = cell_markers[[cell_type]],pt.size = 0,
                  group.by = "res.0.8") +
        ggtitle(paste("VlnPlotfor", cell_type)))
  dev.off()
}




#pdf("Featureplot.markerGene.final.reso0.5.pdf",width=20,height=28)
##FeaturePlot(immune.combinedReso1,features = combined_markers,reduction = "umap",cols=c("grey","yellow","red"))
#dev.off()

pdf("VlnPlot.markerGenes.final.reso_0.5.pdf",width=25,height=30)
VlnPlot(immune.combinedReso1,features = combined_markers,pt.size = 0)
dev.off()


pdf("DotPlot.markerGenes.final.reso_0.5.pdf",width=12,height=10)
DotPlot(immune.combinedReso1,features = markerGenefinal,cols=c("blue","red")) + RotatedAxis()
dev.off()




###DEG_mark
colnames(immune.combinedReso1@meta.data)
Idents(immune.combinedReso1) <- 'RNA_snn_res.0.3'


degs <- FindAllMarkers(immune.combinedReso1,logfc.threshold = log(1.25),only.pos = T,max.cells.per.ident = 100)
degs <- degs %>% filter(p_val_adj < 0.05)
write.csv(degs,file = './res.0.5.csv')
#degs <- setdiff(degs$gene,grep('^MT-',degs$gene,value = T))
#a <- grep('^MT-|^RP',degs$gene,value = T)
#degs <- degs[!degs$gene %in% a,]
top5 <- degs %>% group_by(cluster) %>% top_n(5,-p_val_adj)
top10 <- degs %>% group_by(cluster) %>% top_n(10,-p_val_adj)
write.csv(top10,file = './res.0.5.top10.csv')



immune.combinedReso1$celltype <- factor(immune.combinedReso1$celltype, levels = c("Cardiomyocytes", "Fibroblasts", "Endothelial",'Neuronal_cells','Undefined1','Undefined2'))
pdf("DEG_Feature.0.5.pdf",width=12,height=12)
DotPlot(immune.combinedReso1, features = unique(top5$gene), cols = c('grey','red'))+
  coord_flip()+theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title = element_blank())+scale_color_gradientn(colours = BlueAndRed(10)) + labs(title = 'data')
dev.off()


Idents(immune.combinedReso1) <- 'res.0.8'
levels(Idents(immune.combinedReso1))

immune.combinedReso1 <- RenameIdents(
  immune.combinedReso1,
  '1' = 'Hep 1',
  '2' = 'αβ T cells',
  '3' = 'Hep 2',
  '4' = 'Inflammatory Macs',
  '5' = 'Hep 3',
  '6' = 'Hep 4',
  '7' = 'Plasma cells',
  '8' = 'NK-like cells',
  '9' = 'γδ T cells 1',
  '10' = 'Non-inflammatory Macs',
  '11' = 'Periportal LSECs',
  '12' = 'Central venous LSECs',
  '13' = 'Portal endothelial cells',
  '14' = 'Hep 5',
  '15' = 'Hep 6',
  '16' = 'Mature B cells',
  '17' = 'Cholangiocytes',
  '18' = 'γδ T cells',
  '19' = 'Erythroid cells',
  '20' = 'Hepatic stellate cells'
)
immune.combinedReso1$celltype <- immune.combinedReso1@active.ident
saveRDS(immune.combinedReso1, file = "./analysis_celltype.V1.rds")


pdf("TSNE.celtype.pdf",width=8,height=6)
DimPlot(immune.combinedReso1, label = FALSE, repel = TRUE,group.by="celltype", reduction = "tsne")+scale_color_igv()
dev.off()





###比例
library(dplyr)
library(reshape2)
pdf("bili.celltype..pdf", width = 6, height = 6)  # 设置 PDF 输出尺寸
# 数据处理和比例计算
data_for_plot <- dplyr::select(immune.combinedReso1@meta.data, c('celltype', 'batch')) %>%
  table() %>%
  apply(., 2, function(x) x / sum(x)) %>%  # 计算每个 dataset 中各 celltype 的比例
  data.frame() %>%
  tibble::rownames_to_column('celltype') %>%
  melt(id = 'celltype') %>%
  dplyr::rename('dataset' = 'variable', 'proportion' = 'value')
# 固定 dataset 的顺序
data_for_plot$dataset <- factor(data_for_plot$dataset,
                            levels = c("control", "LIFSI_5ug", "LIFSI_50ug"))
# 绘制竖向堆叠柱状图
ggplot(data_for_plot, aes(x = dataset, y = proportion, fill = celltype)) +
  geom_bar(stat = 'identity', position = 'stack') +  # 堆叠柱状图
  theme_bw() +
  scale_fill_manual(values = cols) +  # 手动设置颜色
  labs(title = "Proportion of Cell Types in Each Dataset",
       x = "Dataset", y = "Proportion", fill = "Cell Type") +  # 修改标签
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转横坐标标签
dev.off()