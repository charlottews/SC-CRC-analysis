# SC-CRC-analysis
Implementation of a single-cell RNA-seq analysis workflow using Seurat, applied to colorectal tumorigenesis dataset GSE161277. 


1.Data soure
The dataset analyzed in this project is GSE161277 (Gene Expression Omnibus, GEO), consisting of single-cell RNA-seq profiles from colorectal cancer patients.
Samples include multiple tissue conditions (normal, adenoma, carcinoma, para-cancer, and blood) from patient-matched biopsies.

2.Analysis Pipeline
The main steps are: Quality Control (QC), Cell Cycle Scoring & Doublet Detection, Clustering & Dimension Reduction, Cell Type Annotation, Differential Expression (DEG) Analysis

3.The analysis was performed using Seurat (Hao et al., Nature Biotechnology, 2023) with reference to best practices in single-cell analysis (Zheng et al., 2022).

Workflow Reference:
Zheng X, Song J, Yu C, et al. Single-cell transcriptomic profiling unravels the adenoma-initiation role of protein tyrosine kinases during colorectal tumorigenesis. Signal Transduct Target Ther. 2022;7(1):60. Published 2022 Feb 28. doi:10.1038/s41392-022-00881-8
