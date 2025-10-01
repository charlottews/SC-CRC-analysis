# 加载必要的包
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


# 1. 数据读取函数
read_sc_data <- function() {
  file_list <- list.files(pattern = "\\.(tsv\\.gz|mtx\\.gz|tar)$")
  sample_ids <- unique(gsub("_(barcodes|features|matrix)\\.(tsv\\.gz|mtx\\.gz)$", "", file_list))
  
  seurat_list <- list()
  
  for (sample_id in sample_ids) {
    barcodes_file <- paste0(sample_id, "_barcodes.tsv.gz")
    features_file <- paste0(sample_id, "_features.tsv.gz")
    matrix_file <- paste0(sample_id, "_matrix.mtx.gz")
    
    if (all(file.exists(c(barcodes_file, features_file, matrix_file)))) {
      cat("Processing:", sample_id, "\n")
      
      counts <- ReadMtx(
        mtx = matrix_file,
        cells = barcodes_file,
        features = features_file,
        feature.column = 2
      )
      
      seurat_obj <- CreateSeuratObject(
        counts = counts,
        project = sample_id,
        min.cells = 3,
        min.features = 200
      )
      
      # 添加样本元数据
      meta_info <- strsplit(sample_id, "_")[[1]]
      seurat_obj$sample_id <- sample_id
      seurat_obj$patient <- meta_info[2]
      seurat_obj$tissue_type <- meta_info[3]
      
      seurat_list[[sample_id]] <- seurat_obj
    }
  }
  return(seurat_list)
}

# 2. 质量控制函数
perform_qc <- function(seurat_list, output_dir) {
  setwd(output_dir)
  
  qc_params <- list(
    minGene = 600,
    maxGene = 20000,
    minUMI = 6000,
    maxUMI = 50000,
    pctMT = 5
  )
  
  MT.genes <- c("ND3", "ND5", "ND2", "COX1", "ND6", "ND1", "ATP6", "COX3",
                "ND4L", "COX2", "CYTB", "ATP8", "ND4")
  
  seurat_qc_list <- list()
  
  for (sample_id in names(seurat_list)) {
    obj <- seurat_list[[sample_id]]
    cat("QC for:", sample_id, "\n")
    
    # 计算线粒体基因百分比
    mt_genes_found <- intersect(toupper(MT.genes), toupper(rownames(obj)))
    if (length(mt_genes_found) > 0) {
      obj[["percent.MT"]] <- PercentageFeatureSet(obj, features = rownames(obj)[toupper(rownames(obj)) %in% toupper(MT.genes)])
    } else {
      warning("No MT genes found in ", sample_id)
      obj[["percent.MT"]] <- 0
    }
    
    # 生成QC前图
    pdf(paste0(sample_id, ".QC.before.pdf"), width = 8, height = 6)
    print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), 
                  ncol = 3, pt.size = 0.1))
    dev.off()
    
    # 筛选数据
    obj_qc <- subset(obj, subset = 
                      nCount_RNA < qc_params$maxUMI &
                      nCount_RNA > qc_params$minUMI &
                      nFeature_RNA > qc_params$minGene &
                      nFeature_RNA < qc_params$maxGene &
                      percent.MT < qc_params$pctMT)
    
    # 生成QC后图
    pdf(paste0(sample_id, ".QC.after.pdf"), width = 8, height = 6)
    print(VlnPlot(obj_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), 
                  ncol = 3, pt.size = 0.1))
    dev.off()
    
    # 保存结果
    seurat_qc_list[[sample_id]] <- obj_qc
    saveRDS(obj_qc, file = paste0(sample_id, "_qc.rds"))
  }
  return(seurat_qc_list)
}

# 3. 细胞周期评分和双联体检测函数
process_cell_cycle_and_doublets <- function(seurat_qc_list, output_dir, cc_genes_path) {
  setwd(output_dir)
  colors_DF <- list(Doublet_hi = "red", Singlet = "blue")
  
  # 读取细胞周期基因
  cc.genes <- read.table(cc_genes_path)$V1
  s_genes_all <- cc.genes[1:43]
  g2m_genes_all <- cc.genes[44:97]
  
  seurat_singlet_list <- list()
  
  for (sample_id in names(seurat_qc_list)) {
    obj <- seurat_qc_list[[sample_id]]
    cat("\nProcessing cell cycle and doublets for:", sample_id, "\n")
    
    # 细胞周期评分
    s_genes <- intersect(s_genes_all, rownames(obj))
    g2m_genes <- intersect(g2m_genes_all, rownames(obj))
    
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst")
    obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes)
    obj <- ScaleData(obj, features = rownames(obj), vars.to.regress = c("S.Score", "G2M.Score"))
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
    
    # 可视化细胞周期
    pdf(paste0(sample_id, ".cellcycle.pdf"), width = 12, height = 6)
    p1 <- DimPlot(obj, reduction = "pca", group.by = "Phase") + ggtitle("PCA by Cell Cycle Phase")
    p2 <- DimPlot(obj, reduction = "umap", group.by = "Phase") + ggtitle("UMAP by Cell Cycle Phase")
    print(p1 + p2)
    dev.off()
    
    # 双联体检测
    sweep.res <- paramSweep_v3(obj, PCs = 1:40, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    annotations <- obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    DoubletRate <- ncol(obj) * 8 * 1e-6
    nExp_poi <- round(DoubletRate * ncol(obj))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    # 运行DoubletFinder
    obj <- doubletFinder_v3(obj, PCs = 1:40, pN = 0.25, pK = mpK,
                            nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    
    pANN_name <- grep("pANN", colnames(obj@meta.data), value = TRUE)
    obj <- doubletFinder_v3(obj, PCs = 1:40, pN = 0.25, pK = mpK,
                            nExp = nExp_poi.adj, reuse.pANN = pANN_name, sct = FALSE)
    
    # 安全获取双联体分类列名
    df_cols <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
    
    # 调试信息
    cat("找到的双联体分类列：", paste(df_cols, collapse = ", "), "\n")
    
    # 优先选择调整后的分类列
    target_pattern <- paste0("nExp", nExp_poi.adj)
    adj_col <- grep(target_pattern, df_cols, value = TRUE)
    
    if(length(adj_col) > 0) {
      df_col_use <- adj_col[1]
      cat("使用调整后的分类列：", df_col_use, "\n")
    } else if(length(df_cols) > 0) {
      df_col_use <- tail(df_cols, 1)  # 默认取最后一个
      cat("未找到调整后的分类列，使用最后一个分类列：", df_col_use, "\n")
    } else {
      warning("在样本 ", sample_id, " 中未找到任何双联体分类列，跳过双联体检测")
      next
    }
    
    # 创建新的分类列
    obj@meta.data$DF_hi.lo <- ifelse(
      obj@meta.data[[df_col_use]] == "Doublet", 
      "Doublet_hi", 
      "Singlet"
    )
    
    # 可视化双联体
    pdf(paste0(sample_id, ".doublets.pdf"), width = 8, height = 6)
    print(DimPlot(obj, reduction = "umap", group.by = "DF_hi.lo", cols = unlist(colors_DF)) +
            ggtitle("Doublet Detection"))
    dev.off()
    
    # 筛选单细胞并保存
    obj_singlet <- subset(obj, subset = DF_hi.lo == "Singlet")
    seurat_singlet_list[[sample_id]] <- obj_singlet
    saveRDS(obj_singlet, file = paste0(sample_id, "_singlet.rds"))
  }
  return(seurat_singlet_list)
}

# 主执行流程

seurat_list <- read_sc_data()
seurat_qc_list <- perform_qc(seurat_list,)
seurat_singlet_list <- process_cell_cycle_and_doublets(
  seurat_qc_list
)

# 保存最终结果
saveRDS(seurat_singlet_list, "all_samples_singlets.rds")