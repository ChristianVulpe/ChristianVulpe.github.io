``` r
library(Seurat)
library(dplyr)
library(patchwork)
library(randomcoloR)
library(ggmin)
library(pheatmap)
library(harmony)
library(tidyr)
library(scales)
library(tidyverse)
library(OmnipathR)
library(decoupleR)
library(msigdbr)
library(DESeq2)
library(biomaRt)
library(GENIE3)
library(sva)
library(gridExtra)
library(igraph)
library(tibble)
library(purrr)
library(reshape2)
library(ggraph)
library(tidygraph)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(GSVA)
library(clustifyr)
library(viridis)
library(monocle3)
library(SeuratWrappers)
library(ggpubr)
library(RColorBrewer)
library(presto)
library(scales)
```

# Files

``` r
folders <- c("Hu_NSCLC_GSE207422", "Lui_NSCLC_GSE179994", "Mathewson_GBM_GSE163108",
             "Qiu-Roth_PDAC_GSE205049", "Sade-Feldman_Melanoma_GSE120575", "Yost_BCC-SCC_GSE123813",
             "Caushi_NSCLC_GSE176021")

files <- unlist(lapply(folders, function(f) {
  list.files(f, pattern = "\\.rds$", full.names = TRUE)
}))
```

# CD8 Subsets

## According to articles annotations

``` r
# Qiu-Roth_PDAC_GSE205049

qiu_roth <- readRDS(files[5])
qiu_roth_cd8 <- subset(qiu_roth, subset = major == "CD8 T")
rm(qui_roth)

# Yost_BCC-SCC_GSE123813

yost <- readRDS(files[7])
yost_cd8 <- subset(yost, subset = grepl("^CD8", cluster))
rm(yost)

# Mathewson_GBM_GSE163108

mathewson_10X <- readRDS(files[3])
mathewson_10X_cd8 <- subset(mathewson_10X, subset = annotate_Tcelltype == "CD8+")
rm(mathewson_10X)

mathewson_smartseq2 <- readRDS(files[4])
mathewson_smartseq2_cd8 <- subset(mathewson_smartseq2, subset = annotate == "CD8+")
rm(mathewson_smartseq2)

# Lui_NSCLC_GSE179994

lui <- readRDS(files[2])
lui_cd8 <- subset(lui, subset = cluster %in% c("Non-exhausted", "Prolif.", "Tex"))
rm(lui)

# Caushi_NSCLC_GSE176021

caushi <- readRDS(files[8])
caushi_cd8 <- subset(caushi, subset = grepl("^CD8", CellType))
rm(caushi)

# Wells_healthy_donnors

health_donors <- readRDS("Wells_healthy_donnors/T_cell_only_seurat_object_HD.rds")
health_donors_cd8 <- subset(health_donors, subset = grepl("^cd8_", author_cell_type))
rm(health_donors)
```

## According to CD3 and CD8 expression

``` r
# Hu_NSCLC_GSE207422

hu <- readRDS(files[1])
hu_cd8 <- subset(hu, cells = WhichCells(hu, expression = CD3G > 0 & CD3D > 0 & CD3E > 0 & CD8A > 0 & CD8B > 0))
rm(hu)

# Sade-Feldman_Melanoma_GSE120575

sade_feldman <- readRDS(files[6])
sade_feldman_cd8 <- subset(sade_feldman, cells = WhichCells(sade_feldman, expression = CD3G > 0 & CD3D > 0 & CD3E > 0 & CD8A > 0 & CD8B > 0))
rm(sade_feldman)
```

# Application of the TRM SubPopulation Signature on the CD8 subsets

## Re-clustering with greater resolution

``` r
# Qiu-Roth_PDAC_GSE205049

qiu_roth_cd8 <- NormalizeData(qiu_roth_cd8)
qiu_roth_cd8 <- ScaleData(qiu_roth_cd8, vars.to.regress = c("nFeature_RNA", "percent.mt"))
qiu_roth_cd8 <- FindVariableFeatures(qiu_roth_cd8, nfeatures = 2000)
qiu_roth_cd8 <- RunPCA(qiu_roth_cd8, features = VariableFeatures(object = qiu_roth_cd8))
qiu_roth_cd8 <- RunHarmony(qiu_roth_cd8, c("orig.ident", "DiseaseState"))
qiu_roth_cd8 <- FindNeighbors(qiu_roth_cd8, dims = 1:30, reduction = "harmony")
qiu_roth_cd8 <- FindClusters(qiu_roth_cd8, resolution = 2)
qiu_roth_cd8 <- RunUMAP(qiu_roth_cd8, dims = 1:30, reduction = "harmony", reduction.name = "Harmony.umap")

# Yost_BCC-SCC_GSE123813

yost_cd8[["RNA"]] <- split(yost_cd8[["RNA"]], f = yost_cd8$patient)
yost_cd8 <- NormalizeData(yost_cd8)
yost_cd8 <- FindVariableFeatures(yost_cd8)
yost_cd8 <- ScaleData(yost_cd8)
yost_cd8 <- RunPCA(yost_cd8, features = VariableFeatures(object = yost_cd8))
yost_cd8 <- IntegrateLayers(
  object = yost_cd8, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = T)
yost_cd8 <- JoinLayers(yost_cd8)
yost_cd8 <- FindNeighbors(yost_cd8, dim = 1:30, reduction = "integrated.rpca")
yost_cd8 <- FindClusters(yost_cd8, resolution = 2)
yost_cd8 <- RunUMAP(yost_cd8, dim = 1:30, reduction = "integrated.rpca", reduction.name = "RPCA.umap")

# Mathewson_GBM_GSE163108

mathewson_smartseq2_cd8 <- NormalizeData(mathewson_smartseq2_cd8)
mathewson_smartseq2_cd8[["RNA"]]$data <- log1p(mathewson_smartseq2_cd8[["RNA"]]$counts) # car des TPM
mathewson_smartseq2_cd8 <- ScaleData(mathewson_smartseq2_cd8)
mathewson_smartseq2_cd8 <- FindVariableFeatures(mathewson_smartseq2_cd8, nfeatures = 2000)
mathewson_smartseq2_cd8 <- RunPCA(mathewson_smartseq2_cd8, features = VariableFeatures(object = mathewson_smartseq2_cd8))
mathewson_smartseq2_cd8 <- FindNeighbors(mathewson_smartseq2_cd8, dims = 1:30)
mathewson_smartseq2_cd8 <- FindClusters(mathewson_smartseq2_cd8, resolution = 2)
mathewson_smartseq2_cd8 <- RunUMAP(mathewson_smartseq2_cd8, dims = 1:30)

mathewson_10X_cd8[["RNA"]] <- split(mathewson_10X_cd8[["RNA"]], f = mathewson_10X_cd8$sampleid)
mathewson_10X_cd8 <- NormalizeData(mathewson_10X_cd8)
mathewson_10X_cd8 <- ScaleData(mathewson_10X_cd8)
mathewson_10X_cd8 <- FindVariableFeatures(mathewson_10X_cd8, nfeatures = 2000)
mathewson_10X_cd8 <- RunPCA(mathewson_10X_cd8, features = VariableFeatures(object = mathewson_10X_cd8))
mathewson_10X_cd8 <- IntegrateLayers(
  object = mathewson_10X_cd8, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = T)
mathewson_10X_cd8 <- JoinLayers(mathewson_10X_cd8)
mathewson_10X_cd8 <- FindNeighbors(mathewson_10X_cd8, dims = 1:30, reduction = "integrated.rpca")
mathewson_10X_cd8 <- FindClusters(mathewson_10X_cd8, resolution = 0.9)
mathewson_10X_cd8 <- RunUMAP(mathewson_10X_cd8, dims = 1:30, reduction = "integrated.rpca", reduction.name = "RPCA.umap")

# Lui_NSCLC_GSE179994

lui_cd8 <- NormalizeData(lui_cd8)
lui_cd8 <- ScaleData(lui_cd8)
lui_cd8 <- FindVariableFeatures(lui_cd8, nfeatures = 2000)
lui_cd8 <- RunPCA(lui_cd8, features = VariableFeatures(object = lui_cd8))
lui_cd8 <- RunHarmony(lui_cd8, c("patient", "sample"))
lui_cd8 <- FindNeighbors(lui_cd8, dims = 1:30, reduction = "harmony")
lui_cd8 <- FindClusters(lui_cd8, resolution = 0.8)
lui_cd8 <- RunUMAP(lui_cd8, dims = 1:30, reduction = "harmony", reduction.name = "Harmony.umap")

# Caushi_NSCLC_GSE176021

caushi_cd8 <- NormalizeData(caushi_cd8)
caushi_cd8 <- FindVariableFeatures(caushi_cd8, nfeatures = 3000)
ifn_related_genes <- read_csv("Caushi_NSCLC_GSE176021/raw_data/others/ifn.csv")
mitochondrial_related_genes <- read_csv("Caushi_NSCLC_GSE176021/raw_data/others/mito.csv")
ig_related_genes <- read_csv("Caushi_NSCLC_GSE176021/raw_data/others/Ig.csv")
blacklist <- c(ifn_related_genes$IFN_module,
               mitochondrial_related_genes$mito,
               ig_related_genes$Ig_module)
VariableFeatures(caushi_cd8) <- setdiff(VariableFeatures(caushi_cd8), blacklist)

caushi_cd8 <- ScaleData(caushi_cd8)
caushi_cd8 <- RunPCA(caushi_cd8, features = setdiff(VariableFeatures(caushi_cd8), blacklist))
caushi_cd8 <- RunBBKNN(caushi_cd8, batch_key = "sample")
caushi_cd8 <- FindClusters(caushi_cd8, resolution = 1)

min_cells <- 200

pca_data <- Embeddings(caushi_cd8, "pca")[, 1:30]
Idents(caushi_cd8) <- "RNA_snn_res.1"
clusters <- Idents(caushi_cd8)

centroids <- aggregate(pca_data, by=list(cluster=clusters), FUN=mean)
rownames(centroids) <- centroids$cluster
centroids <- centroids[, -1]  

cluster_sizes <- table(clusters)
small_clusters <- names(cluster_sizes[cluster_sizes < min_cells])
large_clusters <- names(cluster_sizes[cluster_sizes >= min_cells])

new_clusters <- as.character(clusters)  

for (cl in small_clusters) {
  small_centroid <- as.numeric(centroids[cl, ])

  large_centroids <- centroids[large_clusters, , drop=FALSE]

  dists <- apply(large_centroids, 1, function(x) sqrt(sum((x - small_centroid)^2)))

  closest_cluster <- names(which.min(dists))
  new_clusters[new_clusters == cl] <- closest_cluster
}
caushi_cd8[[]]$seurat_clusters <- ifelse(new_clusters == "45", "40", new_clusters)
caushi_cd8[[]]$seurat_clusters <- factor(new_clusters)


# Hu_NSCLC_GSE207422

hu_cd8 <- NormalizeData(hu_cd8)
hu_cd8 <- ScaleData(hu_cd8)
hu_cd8 <- FindVariableFeatures(hu_cd8, nfeatures = 2000)
hu_cd8 <- RunPCA(hu_cd8, features = VariableFeatures(object = hu_cd8))
hu_cd8 <- RunHarmony(hu_cd8, c("Sample", "Patient"))
hu_cd8 <- FindNeighbors(hu_cd8, dims = 1:30, reduction = "harmony")
hu_cd8 <- FindClusters(hu_cd8, resolution = 0.3)
hu_cd8 <- RunUMAP(hu_cd8, dims = 1:30, reduction = "harmony", reduction.name = "Harmony.umap")

# Sade-Feldman_Melanoma_GSE120575

sade_feldman_cd8 <- NormalizeData(sade_feldman_cd8) 
sade_feldman_cd8[["RNA"]]$data <- sade_feldman_cd8[["RNA"]]$counts
sade_feldman_cd8 <- ScaleData(sade_feldman_cd8)
sade_feldman_cd8 <- FindVariableFeatures(sade_feldman_cd8, nfeatures = 2000)  
sade_feldman_cd8 <- RunPCA(sade_feldman_cd8, features = VariableFeatures(object = sade_feldman_cd8))
sade_feldman_cd8 <- FindNeighbors(sade_feldman_cd8, dims = 1:30)
sade_feldman_cd8 <- FindClusters(sade_feldman_cd8, resolution = 0.6)
sade_feldman_cd8 <- RunUMAP(sade_feldman_cd8, dims = 1:30)
```

### First UMAP visualization

``` r
plotting <- function(object, red, title){
  gg <- DimPlot(object, reduction = red, group.by = "seurat_clusters", cols = distinctColorPalette(length(unique(object[[]]$seurat_clusters))))
  ggsave(title, gg, width = 10, height = 5)
}

plotting(qiu_roth_cd8, "Harmony.umap", "qiu_roth_cd8.png")
plotting(yost_cd8, "RPCA.umap", "yost_cd8.png")
plotting(mathewson_smartseq2_cd8, "umap", "mathewson_smartseq2_cd8.png")
plotting(mathewson_10X_cd8, "RPCA.umap", "mathewson_10X_cd8.png")
plotting(lui_cd8, "Harmony.umap", "lui_cd8.png")
plotting(caushi_cd8, "umap", "caushi_cd8.png")
plotting(hu_cd8, "umap", "hu_cd8.png")
plotting(sade_feldman_cd8, "umap", "sade_feldman_cd8.png")
```

## TRM and SubPop annotation

### Main function

``` r
# Loading results 

load("sig_mat.RData")

files <- list.files("CD8 subsets", pattern = "\\.rds$", full.names = TRUE)
clean_names <- sub("\\.rds$", "", basename(files))
for (i in seq_along(files)) {
  assign(clean_names[i], readRDS(files[i]), envir = .GlobalEnv)
}

# Signatures of interest

sig <- readRDS("TRM_signatures.rds")

trm_articles_signature_genes <- c(
  "CD69", "ITGAE", "ITGA1", "CD101", "CXCR6", "CXCR3", "CXCR4",
  "IL7R", "PDCD1", "ENTPD1", "CD244", "CD160", "LAG3", "HAVCR2",
  
  "ZNF683",  
  "PRDM1",   
  "RUNX3",
  "BHLHE40",
  "NR4A1", "NR4A2", "NR4A3",
  "TOX", "ID2", "BATF",
  "HOPX", "EOMES",
  
  "GZMB", "GZMA", "PRF1", "IFNG", "TNF", "CCL3", "CCL4", "CCL5",
  "XCL1", "XCL2", "CXCL13",
  "S100A4", "S100A6", "S100A10", "ANXA1", "ANXA2",
  
  "FABP5", "FABP4", "S1PR2", "CD9", "CD44",
  "LGALS1", "LGALS3", "ITGB1", "IL2RB", "RGS1",
  "KLF3", "KLF6",
  
  "MALAT1", "MT2A", "MT1E", "MT1X",
  "HSP90AB1", "HSPB1", "NFKBIA",
  
  "PLAUR", "CSF1", "CTSW", "LYAR", "IFITM1", "IFITM2", "IFITM3",
  "TNFAIP3", "TNFRSF18", "TNFRSF4", "TNFRSF9", "IL2RG", "CD2",
  "CD3D", "CD3E", "CD3G", "CD8A", "CD8B"
)

sig[["Articles_Signatures"]] <- trm_articles_signature_genes

sig_ifn = list(IFN = c("IFNAR1", "IFNAR2", "STAT1", "STAT2", "IRF9", "MX1", "MX2", 
            "OAS1", "OAS2", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IRF1", 
            "IRF7", "IRF8", "LAG3", "PDCDA", "TIGIT"))

metabolic <- readRDS("metabolism_short_version.rds")
tumor_reactivity <- readRDS("tumor_reactivity.rds")

all_sig <- c(sig, sig_ifn, metabolic, tumor_reactivity)

trm_ifn_scoring <- function(seu){
  expr <- as.matrix(GetAssayData(seu, layer = "data"))
  clusters <- seu[[]]$seurat_clusters
  scores <- gsva(ssgseaParam(expr, all_sig), verbose = T)
  
  seu[[]]$IFN <- scores[11,]
  seu[[]]$Cytotoxic <- scores[5,]
  seu[[]]$NeoTCR_spe_CD8 <- scores[9,]
  metab <- t(scores[12:20,])
  tumor_react <- t(scores[21:24,])
  
  seu[[]] <- cbind(seu[[]], metab, tumor_react)
  
  scores <- scores[-c(11, 5, 9, 4, 2, 8, 12:24),]
  
  scores_df <- as.data.frame(t(scores))
  scores_df$cluster <- clusters
  
  cluster_means <- scores_df %>%
    group_by(cluster) %>%
    summarise(across(where(is.numeric), mean, .names = "mean_{.col}"))
  
  global_medians <- scores_df %>%
    summarise(across(where(is.numeric), median))

  comparison_df <- cluster_means
  
  for (path in names(global_medians)) {
    mean_col <- paste0("mean_", path)
    if (mean_col %in% names(cluster_means)) {
      comparison_df[[path]] <- ifelse(
        cluster_means[[mean_col]] > global_medians[[path]], 
        "High", 
        "Low"
      )
    }
  }

  comparison_df <- comparison_df %>%
    select(cluster, all_of(names(global_medians)))
  
  cluster_votes <- comparison_df %>%
    mutate(TRM_status = apply(select(., -cluster), 1, function(x) {
      if (sum(x == "High") > sum(x == "Low")) "TRM_high" else "TRM_low"
    }))

  cluster_votes <- cluster_votes %>%
    select(cluster, TRM_status)
  
  low <- as.character(cluster_votes$cluster[which(cluster_votes$TRM_status == "TRM_low")])
  
  seu[[]]$TRM <- ifelse(seu[[]]$seurat_clusters %in% low, "Other", "TRM")
  
  return(seu)
}

auto_subcluster_plot <- function(seurat_obj, reduction, threshold = 1.22, final_png) {
  set.seed(123)
  
  seurat_obj <- trm_ifn_scoring(seurat_obj)
  
  final_annot <- clustify(
    seurat_obj,
    sig_mat,
    cluster_col = "seurat_clusters",
    query_genes = rownames(sig_mat),
    n_genes = 0,
    n_perm = 1000,
    threshold = threshold,
    obj_out = TRUE,
    rename_prefix = "TRM_SubPop")
  
  var_group_by <- "TRM_SubPop_type"
  annotations <- final_annot@meta.data[[var_group_by]] 
  final_annot@meta.data[[var_group_by]] <- ifelse(annotations == "unassigned-CLASH!", "unassigned",
                                                  ifelse(grepl("-CLASH!", annotations) & annotations != "unassigned-CLASH!",
                                                         "Ambiguous", annotations))
  
  final_annot@meta.data[[var_group_by]] <- ifelse(final_annot[[]]$TRM == "Other", "Other", final_annot@meta.data[[var_group_by]])
  
  colors <- distinctColorPalette(length(unique(final_annot[["TRM_SubPop_type"]][,1])))
  colors_clusters <- distinctColorPalette(length(unique(final_annot[["seurat_clusters"]][,1])))
  
  p_trmpop <- DimPlot(final_annot, reduction = reduction, group.by = "TRM", label = T, raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 20) +
    theme(legend.position = "bottom") +
    labs(title = "TRM annotation") +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 2))
  
  p_subpop <- DimPlot(final_annot, reduction = reduction, group.by = var_group_by, cols = colors,
                      label = T, repel = TRUE, raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 20) +
    theme(legend.position = "bottom") +
    labs(title = "TRM SubPopulation annotation") +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 2))
  
  p_clusters <- DimPlot(final_annot, reduction = reduction, group.by = "seurat_clusters", cols =
                          colors_clusters, label= T, repel = TRUE, raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 20) +
    theme(legend.position = "bottom") +
    labs(title = "Clusters") +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5))
  p_final <- cowplot::plot_grid(p_trmpop, p_subpop, p_clusters, ncol = 3, align = "h")
  
  ggsave(final_png, p_final, width = 15, height = 10)
  
  return(final_annot)
}

umaps <- function(final_annot, reduction, final_png){
  p_10 <- FeaturePlot(final_annot, reduction = reduction, features = "IFN", raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = "IFN") +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_11 <- FeaturePlot(final_annot, reduction = reduction, features = "Cytotoxic", raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = "Cytotoxic") +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_12 <- FeaturePlot(final_annot, reduction = reduction, features = "NeoTCR_spe_CD8", raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = "NeoTCR_spe_CD8") +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  sigs <- c("HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM",
            "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "KEGG_OXIDATIVE_PHOSPHORYLATION",
            "KEGG_FATTY_ACID_METABOLISM", "GOBP_GLUCOSE_METABOLIC_PROCESS",
            "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_FATTY_ACID_METABOLIC_PROCESS")
  
  p_1 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[1], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[1])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_2 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[2], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[2])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_3 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[3], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[3])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_4 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[4], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[4])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_5 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[5], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[5])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_6 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[6], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[6])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_7 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[7], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[7])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_8 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[8], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[8])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_9 <- FeaturePlot(final_annot, reduction = reduction, features = sigs[9], raster = F) +
    ggmin::theme_min() +
    scale_y_continuous(breaks=NULL) +
    scale_x_continuous(breaks=NULL) +
    xlab("UMAP1") + ylab("UMAP2") +
    FontSize(x.title = 20, y.title = 20, main = 16) +
    theme(legend.position = "bottom") +
    labs(title = sub("^[^_]*_", "", sigs[9])) +
    guides(colour = guide_legend(override.aes = list(size=5),
                                 title.theme = element_text(size=15, face="bold"),
                                 title.position = "top",
                                 label.theme = element_text(size=15),
                                 ncol = 5)) +
    scale_colour_viridis_c(option = "plasma")
  
  p_final <- cowplot::plot_grid(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10, p_11, p_12, ncol = 4, align = "h")
  
  ggsave(final_png, p_final, width = 15, height = 15)
}
```

### Application

``` r
qiu_roth_cd8 <- auto_subcluster_plot(qiu_roth_cd8, reduction = "Harmony.umap",
                                     final_png = "Plots/qiu_roth_cd8_CD8_final_annotation.png")
yost_cd8 <- auto_subcluster_plot(yost_cd8, reduction = "RPCA.umap",
                                     final_png = "Plots/yost_cd8_CD8_final_annotation.png")
mathewson_smartseq2_cd8 <- auto_subcluster_plot(mathewson_smartseq2_cd8, reduction = "umap",
                                     final_png = "Plots/mathewson_smartseq2_cd8_CD8_final_annotation.png")
mathewson_10X_cd8 <- auto_subcluster_plot(mathewson_10X_cd8, reduction = "RPCA.umap",
                                     final_png = "Plots/mathewson_10X_cd8_CD8_final_annotation.png")
lui_cd8 <- auto_subcluster_plot(lui_cd8, reduction = "Harmony.umap",
                                     final_png = "Plots/lui_cd8_CD8_final_annotation.png")
caushi_cd8 <- auto_subcluster_plot(caushi_cd8, reduction = "umap",
                                     final_png = "Plots/caushi_cd8_CD8_final_annotation.png")
hu_cd8 <- auto_subcluster_plot(hu_cd8, reduction = "Harmony.umap",
                                     final_png = "Plots/hu_cd8_CD8_final_annotation.png")
sade_feldman_cd8 <- auto_subcluster_plot(sade_feldman_cd8, reduction = "umap",
                                     final_png = "Plots/sade_feldman_cd8_CD8_final_annotation.png")
luiy_cd8 <- auto_subcluster_plot(luiy_cd8, reduction = "umap",
                                     final_png = "Plots/luiy_cd8_CD8_final_annotation.png")
che_cd8 <- auto_subcluster_plot(che_cd8, reduction = "umap",
                                     final_png = "Plots/che_cd8_CD8_final_annotation.png")
chen_cd8 <- auto_subcluster_plot(chen_cd8, reduction = "bbknn_umap",
                                     final_png = "Plots/chen_cd8_CD8_final_annotation.png")
in_vitro_cd8 <- FindClusters(in_vitro_cd8, resolution = 1)
in_vitro_cd8 <- auto_subcluster_plot(in_vitro_cd8, reduction = "umap",
                                final_png = "Plots/in_vitro_cd8_CD8_final_annotation.png")

qiu_roth_cd8_results <- qiu_roth_cd8[[]]
yost_cd8_results <- yost_cd8[[]]
hu_cd8_results <- hu_cd8[[]]
mathewson_smartseq2_cd8_results <- mathewson_smartseq2_cd8[[]]
sade_feldman_cd8_results <- sade_feldman_cd8[[]]
mathewson_10X_cd8_results <- mathewson_10X_cd8[[]]
lui_cd8_results <- lui_cd8[[]]
in_vitro_cd8_results <- in_vitro_cd8[[]]
caushi_cd8_results <- caushi_cd8[[]]
luiy_cd8_results <- luiy_cd8[[]]
che_cd8_results <- che_cd8[[]]
chen_cd8_results <- chen_cd8[[]]

metadata_cd8 <- list(qiu_roth_cd8 = qiu_roth_cd8_results,
                         hu_cd8 = hu_cd8_results,
                         mathewson_smartseq2_cd8 = mathewson_smartseq2_cd8_results,
                         mathewson_10X_cd8 = mathewson_10X_cd8_results,
                         yost_cd8 = yost_cd8_results,
                         sade_feldman_cd8 = sade_feldman_cd8_results,
                         lui_cd8 = lui_cd8_results,
                         in_vitro_cd8 = in_vitro_cd8_results,
                         caushi_cd8 = caushi_cd8_results,
                         luiy_cd8 = luiy_cd8_results,
                         che_cd8 = che_cd8_results,
                         chen_cd8 = chen_cd8_results)

vec <- ifelse(yost_cd8[[]]$patient %in% c("su009","su010-S","su010","su011","su012","su001","su002","su003","su004"), "Yes", "No")
metadata_cd8$yost_cd8$Response <- vec

vec <- ifelse(grepl("pre", lui_cd8[[]]$sample, ignore.case = TRUE), "Pre", 
              ifelse(grepl("post", lui_cd8[[]]$sample, ignore.case = TRUE), "Post", NA))
metadata_cd8$lui_cd8$TimePoint <- vec

vec <- ifelse(lui_cd8[[]]$sample %in% c("P1.pre", "P1.post.1", "P1.post.2", "P10.pre", "P10.post.1", "P13.pre", "P13.post.1", "P19.pre", "P19.post.1", "P29.pre", "P29.post.1", "P30.pre", "P30.post.1", "P33.pre", "P33.post.1", "P35.pre", "P35.post.1"), "Yes", "No")
metadata_cd8$lui_cd8$Response <- vec

vec <- ifelse(grepl("Pre", sade_feldman_cd8[[]]$`Patient ID`, ignore.case = TRUE), "Pre",
              ifelse(grepl("Post", sade_feldman_cd8[[]]$`Patient ID`, ignore.case = TRUE), "Post", NA))
metadata_cd8$sade_feldman_cd8$TimePoint <- vec

vec <- ifelse(mathewson_10X_cd8[[]]$sampleid %in% c("E37_GEX", "E99_GEX"), "adenovirus + IL12", "None")
metadata_cd8$mathewson_10X_cd8$ImmunoTherapy <- vec

save(metadata_cd8, file = "metadata_cd8.RData")


umaps(qiu_roth_cd8, reduction = "Harmony.umap", final_png = "Plots/qiu_roth_cd8_CD8_gene_sets_of_interest.png")
umaps(yost_cd8, reduction = "RPCA.umap", final_png = "Plots/yost_cd8_CD8_gene_sets_of_interest.png")
umaps(mathewson_smartseq2_cd8, reduction = "umap", final_png = "Plots/mathewson_smartseq2_cd8_CD8_gene_sets_of_interest.png")
umaps(mathewson_10X_cd8, reduction = "RPCA.umap", final_png = "Plots/mathewson_10X_cd8_CD8_gene_sets_of_interest.png")
umaps(lui_cd8, reduction = "Harmony.umap", final_png = "Plots/lui_cd8_CD8_gene_sets_of_interest.png")
umaps(caushi_cd8, reduction = "umap", final_png = "Plots/caushi_cd8_CD8_gene_sets_of_interest.png")
umaps(hu_cd8, reduction = "Harmony.umap", final_png = "Plots/hu_cd8_CD8_gene_sets_of_interest.png")
umaps(sade_feldman_cd8, reduction = "umap", final_png = "Plots/sade_feldman_cd8_CD8_gene_sets_of_interest.png")
umaps(luiy_cd8, reduction = "umap", final_png = "Plots/luiy_cd8_CD8_gene_sets_of_interest.png")
umaps(che_cd8, reduction = "umap", final_png = "Plots/che_cd8_CD8_gene_sets_of_interest.png")
umaps(chen_cd8, reduction = "bbknn_umap", final_png = "Plots/chen_cd8_CD8_gene_sets_of_interest.png")
umaps(in_vitro_cd8, reduction = "umap", final_png = "Plots/in_vitro_cd8_CD8_gene_sets_of_interest.png")
```

# Loading of results

``` r
files <- list.files("CD8 subsets", pattern = "\\.rds$", full.names = TRUE)
clean_names <- sub("\\.rds$", "", basename(files))
for (i in seq_along(files)) {
  assign(clean_names[i], readRDS(files[i]), envir = .GlobalEnv)
}

load("metadata_cd8.RData")

for (nm in names(metadata_cd8)) {
  if (exists(nm, envir = .GlobalEnv)) {
    seurat_obj <- get(nm, envir = .GlobalEnv)
    meta <- metadata_cd8[[nm]]

    seurat_obj@meta.data <- meta

    assign(nm, seurat_obj, envir = .GlobalEnv)

    message("Métadonnées mises à jour pour : ", nm)
  } else {
    warning("L’objet Seurat ", nm, " n’existe pas dans l’environnement global")
  }
}
```

# Independence seurat_clusters x response/non-response; pre/post treatment; or else

``` r
independece_plot <- function(prop, dataset_name = "dataset") {
  res <- suppressWarnings(chisq.test(prop))

  stdres <- res$stdres
  df <- as.data.frame(as.table(stdres))
  colnames(df) <- c("Cluster", "Condition", "StdResid")

  df$Count <- as.vector(prop)

  pvals <- 2 * (1 - pnorm(abs(df$StdResid)))
  df$FDR <- p.adjust(pvals, method = "fdr")

  p <- ggplot(df, aes(x = Condition, y = Cluster, fill = StdResid)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-max(abs(df$StdResid)), max(abs(df$StdResid))),
      name = "Résidu\nstandardisé"
    ) +
    geom_text(aes(label = Count), color = "black", size = 4) +
    geom_text(
      aes(label = ifelse(FDR < 0.05, "*", "")),
      color = "black", size = 6, hjust = -2
    ) +
    labs(
      title = paste("Test d'indépendance –", dataset_name),
      subtitle = "Couleurs = résidus standardisés\nRouge = surreprésenté | Bleu = sous-représenté\n* = FDR < 0.05",
      x = "Condition",
      y = "Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  file_name <- paste0("Plots/", gsub(" ", "_", dataset_name), "_chi2_heatmap.png")
  ggsave(file_name, plot = p, width = 8, height = 12, dpi = 300)

  return(list(result = res, data = df, plot = p))
}

# Qiu-Roth_PDAC_GSE205049

prop_qiu_roth_cd8 <- table(qiu_roth_cd8[[]]$seurat_clusters, qiu_roth_cd8[[]]$DiseaseState) 

# Yost_BCC-SCC_GSE123813

prop_yost_cd8_treatment <- table(yost_cd8[[]]$seurat_clusters, yost_cd8[[]]$treatment)
prop_yost_cd8_response <- table(yost_cd8[[]]$seurat_clusters, yost_cd8[[]]$Response) 

# Mathewson_GBM_GSE163108

prop_mathewson_10X_cd8 <- table(mathewson_10X_cd8[[]]$seurat_clusters, mathewson_10X_cd8[[]]$ImmunoTherapy)

prop_mathewson_smartseq2_cd8_therapy <- table(mathewson_smartseq2_cd8[[]]$seurat_clusters, mathewson_smartseq2_cd8[[]]$therapy)
prop_mathewson_smartseq2_cd8_cohort <- table(mathewson_smartseq2_cd8[[]]$seurat_clusters, mathewson_smartseq2_cd8[[]]$Cohort)

# Lui_NSCLC_GSE179994

prop_lui_cd8_treatment <- table(lui_cd8[[]]$seurat_clusters, lui_cd8[[]]$TimePoint) 
prop_lui_cd8_response <- table(lui_cd8[[]]$seurat_clusters, lui_cd8[[]]$Response)

# Caushi_NSCLC_GSE176021

prop_caushi_cd8_response <- table(caushi_cd8[[]]$seurat_clusters, caushi_cd8[[]]$`MPR status`) 
prop_caushi_cd8_tissue <- table(caushi_cd8[[]]$seurat_clusters, caushi_cd8[[]]$tissue) 

# Hu_NSCLC_GSE207422

prop_hu_cd8_treatment <- table(hu_cd8[[]]$seurat_clusters, hu_cd8[[]]$Resource) 
prop_hu_cd8_response <- table(hu_cd8[[]]$seurat_clusters, hu_cd8[[]]$Pathologic_Response) 

# Sade-Feldman_Melanoma_GSE120575

prop_sade_feldman_cd8_response <- table(sade_feldman_cd8[[]]$seurat_clusters, sade_feldman_cd8[[]]$Response) 
prop_sade_feldman_cd8_treatment <- table(sade_feldman_cd8[[]]$seurat_clusters, sade_feldman_cd8[[]]$TimePoint) 

# Chen

prop_chen_cd8_response <- table(chen_cd8[[]]$seurat_clusters, chen_cd8[[]]$Response) 
prop_chen_cd8_treatment <- table(chen_cd8[[]]$seurat_clusters, chen_cd8[[]]$`Treatment Stage`) 
prop_chen_cd8_tissue <- table(chen_cd8[[]]$seurat_clusters, chen_cd8[[]]$Tissue)

# Luiy

prop_luiy_cd8 <- table(luiy_cd8[[]]$seurat_clusters, luiy_cd8[[]]$tissue)

# Che 

prop_che_cd8 <- table(che_cd8[[]]$seurat_clusters, che_cd8[[]]$tissue)


independece_plot(prop_qiu_roth_cd8, "Qiu-Roth CD8")
independece_plot(prop_yost_cd8_treatment, "Yost CD8 Treatment")
independece_plot(prop_yost_cd8_response, "Yost CD8 Response")
independece_plot(prop_mathewson_smartseq2_cd8_cohort, "Mathewson CD8 Smart-Seq2 Cohort")
independece_plot(prop_mathewson_smartseq2_cd8_therapy, "Mathewson CD8 Smart-Seq2 Therapy")
independece_plot(prop_mathewson_10X_cd8, "Mathewson CD8 10X")
independece_plot(prop_lui_cd8_treatment, "Lui CD8 Treatment")
independece_plot(prop_lui_cd8_response, "Lui CD8 Response")
independece_plot(prop_caushi_cd8_tissue, "Caushi CD8 Tissue")
independece_plot(prop_caushi_cd8_response, "Caushi CD8 Response")
independece_plot(prop_hu_cd8_treatment, "Hu CD8 Treatment")
independece_plot(prop_hu_cd8_response, "Hu CD8 Response")
independece_plot(prop_sade_feldman_cd8_treatment, "Sade-Feldman CD8 Treatment")
independece_plot(prop_sade_feldman_cd8_response, "Sade-Feldman CD8 Response")
independece_plot(prop_chen_cd8_treatment, "Chen CD8 Treatment")
independece_plot(prop_chen_cd8_response, "Chen CD8 Response")
independece_plot(prop_chen_cd8_tissue, "Chen CD8 Tissue")
independece_plot(prop_luiy_cd8, "Luiy CD8 Tissue")
independece_plot(prop_che_cd8, "Che CD8 Tissue")
```

# Proportions visualizations

``` r
plot_prop <- function(tab, dataset_name = "dataset") {
  prop_df <- as.data.frame(prop.table(tab, margin = 1))
  colnames(prop_df) <- c("Cluster", "Condition", "Proportion")
  
  totals <- as.data.frame(tab) %>%
    group_by(Var2) %>%
    summarise(Total = sum(Freq)) %>%
    mutate(Label = paste0(Var2, ": ", Total)) %>%
    pull(Label)
  
  n <- ceiling(length(totals) / 2)
  col1 <- totals[1:n]
  col2 <- totals[(n+1):length(totals)]

  if(length(col2) < n) col2 <- c(col2, rep("", n - length(col2)))
  
  totals_text <- paste(paste(col1, col2, sep = "   |   "), collapse = "\n")
  
  p <- ggplot(prop_df, aes(x = Cluster, y = Proportion, fill = Condition)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste("Proportions par cluster –", dataset_name),
      subtitle = paste("Totaux par condition :\n", totals_text),
      x = "Cluster",
      y = "Proportion (%)",
      fill = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"),
      plot.margin = margin(t = 15, r = 10, b = 10, l = 10) 
    )
  
  file_name <- paste0("Plots/", gsub(" ", "_", dataset_name), "_proportions.png")
  ggsave(file_name, plot = p, width = 8, height = 5, dpi = 300)
  
  return(p)
}


plot_prop(prop_qiu_roth_cd8, "Qiu-Roth CD8")
plot_prop(prop_yost_cd8_treatment, "Yost CD8 Treatment")
plot_prop(prop_yost_cd8_response, "Yost CD8 Response")
plot_prop(prop_mathewson_smartseq2_cd8_cohort, "Mathewson CD8 Smart-Seq2 Cohort")
plot_prop(prop_mathewson_smartseq2_cd8_therapy, "Mathewson CD8 Smart-Seq2 Therapy")
plot_prop(prop_mathewson_10X_cd8, "Mathewson CD8 10X")
plot_prop(prop_lui_cd8_treatment, "Lui CD8 Treatment")
plot_prop(prop_lui_cd8_response, "Lui CD8 Reponse")
plot_prop(prop_caushi_cd8_tissue, "Caushi CD8 Tissue")
plot_prop(prop_caushi_cd8_response, "Caushi CD8 Response")
plot_prop(prop_hu_cd8_treatment, "Hu CD8 Treatment")
plot_prop(prop_hu_cd8_response, "Hu CD8 Response")
plot_prop(prop_sade_feldman_cd8_treatment, "Sade-Feldman CD8 Treatment")
plot_prop(prop_sade_feldman_cd8_response, "Sade-Feldman CD8 Response")
plot_prop(prop_chen_cd8_treatment, "Chen CD8 Treatment")
plot_prop(prop_chen_cd8_response, "Chen CD8 Response")
plot_prop(prop_chen_cd8_tissue, "Chen CD8 Tissue")
plot_prop(prop_luiy_cd8, "Luiy CD8 Tissue")
plot_prop(prop_che_cd8, "Che CD8 Tissue")
```

# Bulk data

## Loading

``` r
sample_name <- gsub("_.*", "", list.files("Counts"))
gene_id <- read.delim("Counts/A115_MKRN250026717-1A_22VTKLLT4_L6_.txt", comment = "#")$Geneid

for (i in list.files("Counts/")) {  
  print(i)
  if (i == list.files("Counts/")[1]) {
    data <- read.delim(paste0("Counts/", i), comment = "#")
    data <- data[, 7]
  } else {
    temp <- read.delim(paste0("Counts/", i), comment = "#")
    temp <- temp[, 7]
    data <- cbind(data, temp)
    rm(temp)
  }
}

colnames(data) <- sample_name
rownames(data) <- gene_id
data <- as.data.frame(data)  

colData <- read_tsv("colData.txt")
rownames(colData) <- colData$Replicate
colData$subpopulation <- factor(colData$subpopulation)
colData$sample <- factor(colData$Sample)
colData$HD <- factor(colData$HD)  

rm(list = setdiff(ls(), c("data", "colData")))  

dds_list <- list()

for (sp in unique(colData$subpopulation)) {
  colData_temp <- colData
  colData_temp$current_vs_other <- ifelse(colData_temp$subpopulation == sp, sp, "other")
  colData_temp$current_vs_other <- factor(colData_temp$current_vs_other, levels = c("other", sp))
  
  dds_temp <- DESeqDataSetFromMatrix(
    countData = data,
    colData = colData_temp,
    design = ~ HD + current_vs_other
  )
  
  dds_list[[sp]] <- dds_temp
}

dds_list <- lapply(dds_list, FUN = function(dds){
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep, ]
  dds <- collapseReplicates(dds, dds$sample, dds$Replicate)
  dds <- DESeq(dds)
  dds
})

save(dds_list, file = "dds_list.RData")

results_list <- lapply(dds_list, function(dds) {
  res <- results(dds, contrast = c("current_vs_other", setdiff(dds$current_vs_other, "other"), "other"))       
  res <- as.data.frame(res)          
  res_filtered <- res[!is.na(res$padj) & res$padj <= 0.05, ]  
  res_filtered
})

convert_gene_IDs <- function(DESeq2_results) {
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  DESeq2_results <- as.data.frame(DESeq2_results)
  
  gene_names <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = rownames(DESeq2_results),
    mart = ensembl
  )
  
  DESeq2_results$gene <- gene_names$external_gene_name[match(rownames(DESeq2_results), gene_names$ensembl_gene_id)]
  
  keep <- !is.na(DESeq2_results$gene) & DESeq2_results$gene != ""
  if(sum(!keep) > 0) {
    warning(sum(!keep), " lignes ont été supprimées car le nom de gène est NA ou vide.")
  }
  DESeq2_results <- DESeq2_results[keep, , drop = FALSE]
  
  return(DESeq2_results)
}

for (i in names(results_list)) {
  print(i)
  results_list[[i]] <- convert_gene_IDs(results_list[[i]])
}

save(results_list, file = "results_list.RData")
```

## GRN with GENIE3

``` r
load("dds_list.RData")
load("results_list.RData")

# Prepare expressions matrixes for GRN

vst_list <- lapply(names(dds_list), FUN = function(sp) {
  dds <- dds_list[[sp]]
  samples_sp <- rownames(colData(dds))
  vst <- vst(dds, blind = T)
  mat <- assay(vst)
  mat_sp <- mat[, dds$subpopulation == sp, drop = FALSE]
  return(mat_sp)
})

convert_vst_mat <- function(vst_mat) {
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  gene_names <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = rownames(vst_mat),
    mart = ensembl
  )
  
  gene_map <- setNames(gene_names$external_gene_name, gene_names$ensembl_gene_id)
  
  new_names <- gene_map[rownames(vst_mat)]

  keep <- !is.na(new_names) & new_names != ""
  if(sum(!keep) > 0) {
    warning(sum(!keep), " lignes ont été supprimées car non trouvées dans Biomart.")
  }
  
  vst_mat_filtered <- vst_mat[keep, , drop = FALSE]
  rownames(vst_mat_filtered) <- new_names[keep]
  
  return(vst_mat_filtered)
}

vst_list <- lapply(vst_list, convert_vst_mat)
names(vst_list) <- names(dds_list)

save(vst_list, file = "vst_list.RData")

load("vst_list.RData")

# GENIE3 on all genes

vst_list_unique <- lapply(vst_list, function(mat) {
  rownames(mat) <- make.unique(rownames(mat))
  mat
})

for (i in seq_along(vst_list_unique)) {
  assign(names(vst_list_unique)[i], vst_list_unique[[i]], envir = .GlobalEnv)
}

set.seed(123)
wm_neg_neg <- GENIE3(CXCR6negPD1neg, verbose = TRUE)
wm_neg_pos <- GENIE3(CXCR6negPD1pos, verbose = TRUE)     
wm_pos_neg <- GENIE3(CXCR6posPD1neg, verbose = TRUE)
wm_pos_pos <- GENIE3(CXCR6posPD1pos, verbose = TRUE)

res <- list(CXCR6negPD1neg = wm_neg_neg, CXCR6negPD1pos = wm_neg_pos,
            CXCR6posPD1neg = wm_pos_neg, CXCR6posPD1pos = wm_pos_pos)
save(res, file = "GENIE3_res.RData")

load("GENIE3_res.RData")

# GENIE3 on DEG genes

vst_list_filtered <- mapply(function(vst_mat, res_df) {
  genes_keep <- res_df$gene
  
  vst_mat_filtered <- vst_mat[rownames(vst_mat) %in% genes_keep, , drop = FALSE]
  
  vst_mat_filtered
}, vst_list, results_list, SIMPLIFY = FALSE)

vst_list_filtered_unique <- lapply(vst_list_filtered, function(mat) {
  rownames(mat) <- make.unique(rownames(mat))
  mat
})

for (i in seq_along(vst_list_filtered_unique)) {
  assign(names(vst_list_filtered_unique)[i], vst_list_filtered_unique[[i]], envir = .GlobalEnv)
}

wm_neg_neg_DEG <- GENIE3(CXCR6negPD1neg, verbose = TRUE)
wm_neg_pos_DEG <- GENIE3(CXCR6negPD1pos, verbose = TRUE)      
wm_pos_neg_DEG <- GENIE3(CXCR6posPD1neg, verbose = TRUE)
wm_pos_pos_DEG <- GENIE3(CXCR6posPD1pos, verbose = TRUE)

res_DEG <- list(CXCR6negPD1neg = wm_neg_neg_DEG, CXCR6negPD1pos = wm_neg_pos_DEG,
            CXCR6posPD1neg = wm_pos_neg_DEG, CXCR6posPD1pos = wm_pos_pos_DEG)

save(res_DEG, file = "GENIE3_res_DEG.RData")
```

## Data interpretation

### Main 2 functions

``` r
compare_GENIE3_Collectri <- function(genie3_mat, tf_pathways, top_n = 500) {
  
  # limiter les colonnes à des TFs connus
  common_tfs <- intersect(colnames(genie3_mat), names(tf_pathways))
  if (length(common_tfs) == 0) return(NULL)
  
  results <- lapply(common_tfs, function(tf) {
    # top N cibles prédictives du TF selon GENIE3
    top_predicted <- rownames(genie3_mat)[order(genie3_mat[, tf], decreasing = TRUE)][1:top_n]
    
    # cibles connues selon CoLlecTRI
    known_targets <- tf_pathways[[tf]]$target
    
    # intersection
    common_targets <- intersect(top_predicted, known_targets)
    
    # score de recouvrement (Jaccard)
    jaccard <- length(common_targets) / length(unique(c(top_predicted, known_targets)))
    
    tibble(
      TF = tf,
      n_predicted = length(top_predicted),
      n_known = length(known_targets),
      n_common = length(common_targets),
      overlap_ratio = jaccard,
      common_targets = list(common_targets)
    )
  })
  
  bind_rows(results)
}

# Fonction pour tracer un réseau des meilleurs TFs détectés
plot_topTFs_network <- function(topTF_table, genie3_mat, tf_pathways,
                                n_tfs = 10, n_predicted = 5, max_validated = 10) {
  
  # Sélection des meilleurs TFs
  top_tfs <- topTF_table %>%
    slice_max(overlap_ratio, n = n_tfs) %>%
    pull(TF)
  
  edges <- list()
  
  for (tf in top_tfs) {
    if (!tf %in% colnames(genie3_mat)) next
    
    scores <- genie3_mat[, tf]
    top_pred <- names(sort(scores, decreasing = TRUE))
    
    # Récupération des interactions connues avec signe
    known_df <- tf_pathways[[tf]]
    known <- known_df$target
    
    # Intersection prédictions / validées
    common <- intersect(top_pred, known)
    if (length(common) > max_validated) {
      common <- common[1:max_validated]
    }
    
    pred_only <- setdiff(top_pred, known)[1:n_predicted]
    
    edges_tmp <- list()
    
    # Interactions validées (avec signe)
    if (length(common) > 0) {
      edges_tmp[[length(edges_tmp) + 1]] <- data.frame(
        from = tf,
        to = common,
        type = "validated",
        regulation = sapply(common, function(g) {
          val <- known_df$mor[known_df$target == g]
          ifelse(val == 1, "Activation",
                 ifelse(val == -1, "Inhibition", "Unknown"))
        }),
        stringsAsFactors = FALSE
      )
    }
    
    # Interactions prédictives (pas de signe connu)
    if (length(pred_only) > 0) {
      edges_tmp[[length(edges_tmp) + 1]] <- data.frame(
        from = tf,
        to = pred_only,
        type = "predicted",
        regulation = "Unknown",
        stringsAsFactors = FALSE
      )
    }
    
    if (length(edges_tmp) > 0) {
      edges[[tf]] <- do.call(rbind, edges_tmp)
    }
  }
  
  edge_df <- do.call(rbind, edges)
  
  # Création du graphe
  g <- tbl_graph(
    nodes = data.frame(name = unique(c(edge_df$from, edge_df$to))),
    edges = edge_df,
    directed = TRUE
  )
  
  # --- Tracé du réseau ---
  ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(
        color = regulation,
        alpha = ifelse(type == "validated", 0.9, 0.4)
      ),
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      show.legend = TRUE
    ) +
    
    geom_node_point(aes(color = case_when(
      name %in% top_tfs ~ "TF",
      name %in% edge_df$to[edge_df$type == "validated"] ~ "Validated",
      TRUE ~ "Predicted"
    )),
    size = 3.5,
    show.legend = TRUE) +
    
    geom_node_text(
      aes(label = name, fontface = ifelse(name %in% top_tfs, "bold", "plain")),
      repel = TRUE,
      size = 2.8
    ) +
    
    # Légendes des noeuds
    scale_color_manual(
      name = "Node type",
      values = c(
        "TF" = "#C44E52",        # rouge atténué / brique rosée
        "Validated" = "#F4D35E", # doré pâle / ocre clair
        "Predicted" = "steelblue3"
      ),
      labels = c(
        "TF" = "Transcription Factor",
        "Validated" = "Validated Target (CoLlecTRI)",
        "Predicted" = "Predicted Target (GENIE3)"
      ),
      guide = guide_legend(order = 1)
    ) +
    
    # Légende des arêtes (avec lignes colorées)
    scale_edge_color_manual(
      name = "Interaction type",
      values = c(
        "Activation" = "#5CA27C", # vert doux / sauge
        "Inhibition" = "#D1615D",   # rouge atténué
        "Unknown" = "gray70"
      ),
      labels = c(
        "Activation" = "Activation",
        "Inhibition" = "Inhibition",
        "Unknown" = "Unknown"
      ),
      guide = guide_legend(
        order = 2,
        override.aes = list(linetype = 1, shape = NA, size = 1.2)
      )
    ) +
    
    guides(edge_alpha = "none") +
    
    theme_void() +
    ggtitle("Top 10 TF Regulatory Network with Activation/Inhibition")
}
```

### DEG GENIE3

``` r
# Loading

load("GENIE3_res_DEG.RData")

for (i in seq_along(res_DEG)) {
  assign(names(res_DEG)[i], res_DEG[[i]], envir = .GlobalEnv)
}

genie3_list <- list(
  CXCR6neg_PD1neg_DEG = CXCR6negPD1neg,
  CXCR6neg_PD1pos_DEG = CXCR6negPD1pos,
  CXCR6pos_PD1neg_DEG = CXCR6posPD1neg,
  CXCR6pos_PD1pos_DEG = CXCR6posPD1pos
)

net <- get_collectri(organism = "human", split_complexes = FALSE)

tf_pathways <- split(net[, c("target", "mor")], net$source)
tf_list <- names(tf_pathways)

cat("Nombre de TFs connus :", length(tf_list), "\n")

# GENIE3 x CollecTRI

tf_overlap_list <- map(genie3_list, ~compare_GENIE3_Collectri(.x, tf_pathways, top_n = 500))

for (sp in names(tf_overlap_list)) {
  if (!is.null(tf_overlap_list[[sp]])) {
    tf_overlap_list[[sp]] <- tf_overlap_list[[sp]] %>%
      mutate(subpopulation = sp)
  }
}

tf_overlap_all <- bind_rows(tf_overlap_list)

top_TFs <- tf_overlap_all %>%
  group_by(subpopulation) %>%
  arrange(desc(overlap_ratio)) %>%
  slice_head(n = 10)

top_TFs

# Visualization

for (subpop in names(genie3_list)) {
  tf_data <- top_TFs %>% filter(subpopulation == subpop)
  if (nrow(tf_data) == 0) next
  file_name <- paste0("Plots/", subpop, "_network.png")
  
  png(file_name, width = 2500, height = 2000, res = 300)
  print(plot_topTFs_network(
    tf_data,
    genie3_list[[subpop]],
    tf_pathways,
    n_tfs = 10,
    n_predicted = 5
  ))
  dev.off()
}
```

### All genes GENIE3

``` r
# Loading

load("GENIE3_res.RData")

for (i in seq_along(res)) {
  assign(names(res)[i], res[[i]], envir = .GlobalEnv)
}

genie3_list <- list(
  CXCR6neg_PD1neg = CXCR6negPD1neg,
  CXCR6neg_PD1pos = CXCR6negPD1pos,
  CXCR6pos_PD1neg = CXCR6posPD1neg,
  CXCR6pos_PD1pos = CXCR6posPD1pos
)

net <- get_collectri(organism = "human", split_complexes = FALSE)

tf_pathways <- split(net[, c("target", "mor")], net$source)
tf_list <- names(tf_pathways)

cat("Nombre de TFs connus :", length(tf_list), "\n")

# CENIE3 x CollecTRI

tf_overlap_list <- map(genie3_list, ~compare_GENIE3_Collectri(.x, tf_pathways, top_n = 500))

for (sp in names(tf_overlap_list)) {
  if (!is.null(tf_overlap_list[[sp]])) {
    tf_overlap_list[[sp]] <- tf_overlap_list[[sp]] %>%
      mutate(subpopulation = sp)
  }
}

tf_overlap_all <- bind_rows(tf_overlap_list)

top_TFs <- tf_overlap_all %>%
  group_by(subpopulation) %>%
  arrange(desc(overlap_ratio)) %>%
  slice_head(n = 10)

top_TFs

# Visualization

for (subpop in names(genie3_list)) {
  tf_data <- top_TFs %>% filter(subpopulation == subpop)
  if (nrow(tf_data) == 0) next
  file_name <- paste0("Plots/", subpop, "_network.png")
  
  png(file_name, width = 2500, height = 2000, res = 300)
  print(plot_topTFs_network(
    tf_data,
    genie3_list[[subpop]],
    tf_pathways,
    n_tfs = 10,
    n_predicted = 5
  ))
  dev.off()
}
```

### TF modules creation

``` r
load("GENIE3_res.RData")
load("GENIE3_res_DEG.RData")

genie3_comparison <- list()

for (subpop in names(res)) {
  genie3_comparison[[paste0(subpop, "_all")]] <- res[[subpop]]
  genie3_comparison[[paste0(subpop, "_DEG")]] <- res_DEG[[subpop]]
}

robust_TFs_modules <- list()

for (subpop in c("CXCR6negPD1neg", "CXCR6negPD1pos", "CXCR6posPD1neg", "CXCR6posPD1pos")) {
  
  # TF top-N sur chaque approche
  tf_all <- compare_GENIE3_Collectri(genie3_comparison[[paste0(subpop, "_all")]], tf_pathways) %>%
    slice_max(overlap_ratio, n = 500) %>% pull(TF)
  
  tf_deg <- compare_GENIE3_Collectri(genie3_comparison[[paste0(subpop, "_DEG")]], tf_pathways) %>%
    slice_max(overlap_ratio, n = 500) %>% pull(TF)
  
  # TF robustes = présents dans les deux approches
  tf_robust <- intersect(tf_all, tf_deg)
  
  modules <- lapply(tf_robust, function(tf) {
    
    # cibles GENIE3
    top_all <- names(sort(genie3_comparison[[paste0(subpop, "_all")]][, tf], decreasing = TRUE))
    top_deg <- names(sort(genie3_comparison[[paste0(subpop, "_DEG")]][, tf], decreasing = TRUE))
    
    # cibles validées CoLlecTRI
    validated <- tf_pathways[[tf]]
    
    # cibles robustes ET validées
    common_targets <- intersect(top_all, top_deg) %>% intersect(validated$target)
    
    # cibles seulement prédictives
    predicted_only <- setdiff(intersect(top_all, top_deg), validated$target)
    
    # cibles validées non retrouvées par GENIE3
    validated_only <- setdiff(validated$target, intersect(top_all, top_deg))
    
    list(
      TF = tf,
      common_targets = common_targets,        # robustes et validées
      predicted_only = predicted_only,        # robustes mais non validées
      validated_only = validated_only,        # validées mais non robustes
      all_targets = top_all,
      DEG_targets = top_deg,
      validated_targets = validated
    )
  })
  
  names(modules) <- tf_robust
  robust_TFs_modules[[subpop]] <- modules
}

save(robust_TFs_modules, file = "robust_TFs_modules.RData")

# --- Construire un data.frame à partir des modules robustes ---
flatten_results <- function(data_list) {
  bind_rows(lapply(names(data_list), function(subpop) {
    modules <- data_list[[subpop]]
    
    bind_rows(lapply(names(modules), function(tf) {
      item <- modules[[tf]]
      
      # ajouter la régulation si disponible
      regulation_info <- tf_pathways[[tf]]
      mor_map <- setNames(ifelse(regulation_info$mor == 1, "Activation", "Inhibition"), regulation_info$target)
      
      # annotation des gènes cibles avec signe
      annotated_targets <- sapply(item$common_targets, function(g) {
        if (g %in% names(mor_map)) {
          paste0(g, " (", mor_map[g], ")")
        } else {
          g
        }
      })
      
      data.frame(
        Sous_population = subpop,
        Facteur_Transcription = tf,
        Nb_genes_cibles = length(item$common_targets),
        Genes_cibles = paste(annotated_targets, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
  }))
}


robust_TFs_df <- flatten_results(robust_TFs_modules)

# ---- Export Excel ----
xlsx_file <- "Robust_TF_Modules.xlsx"
wb <- createWorkbook()

for (sp in unique(robust_TFs_df$Sous_population)) {
  
  df_sub <- robust_TFs_df %>%
    filter(Sous_population == sp) %>%
    arrange(desc(Nb_genes_cibles))  
  
  addWorksheet(wb, sp)
  writeData(wb, sp, df_sub, headerStyle = createStyle(textDecoration = "bold"))
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)
```

# Single sample Enrichment analysis using ssGSEA

## MSigDB Pathways: Hallmarks + C2 + C7

``` r
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
c2 <- msigdbr(species = "Homo sapiens", category = "C2")
c7 <- msigdbr(species = "Homo sapiens", category = "C7")

hallmark_list <- hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
c2_list <- c2 %>% split(x = .$gene_symbol, f = .$gs_name)
c7_list <- c7 %>% split(x = .$gene_symbol, f = .$gs_name)

all_pathways <- c(hallmark_list, c2_list, c7_list)

all_pathways <- lapply(all_pathways, unique)
all_pathways <- all_pathways[lengths(all_pathways) >= 5]
```

## TF patwhays

``` r
load("robust_TFs_modules.RData")

flatten_results <- function(data_list) {
  bind_rows(lapply(names(data_list), function(subpop) {
    modules <- data_list[[subpop]]
    
    bind_rows(lapply(names(modules), function(tf) {
      item <- modules[[tf]]
      data.frame(
        Sous_population = subpop,
        Facteur_Transcription = tf,
        Nb_genes_cibles = length(item$common_targets),
        Genes_cibles = paste(item$common_targets, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

robust_TFs_df <- flatten_results(robust_TFs_modules)

filtered_df <- robust_TFs_df %>%
  filter(Nb_genes_cibles >= 5)

modules_gene_sets_by_subpop <- lapply(
  split(filtered_df, filtered_df$Sous_population),
  function(df_sub) {
    lapply(
      split(df_sub, df_sub$Facteur_Transcription),
      function(x) unique(unlist(strsplit(x$Genes_cibles, ";\\s*")))
    )
  }
)

tf_pathways <- unlist(
  lapply(names(modules_gene_sets_by_subpop), function(sp) {
    tf_list <- modules_gene_sets_by_subpop[[sp]]

    setNames(tf_list, paste0(sp, "_", names(tf_list)))
  }),
  recursive = FALSE
)
```

## Main function

``` r
ssgsea_seurat_pipeline <- function(seurat_obj,
                                   all_pathways = all_pathways,
                                   tf_pathways = tf_pathways,
                                   output_dir = "results_ssgsea",
                                   ssgseaParam_value = 0.25,
                                   prefix = "ssgsea") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Utilisation log-normalisée
  assay_name <- "RNA"
  expr_mat <- as.matrix(GetAssayData(seurat_obj, assay = assay_name, layer = "data"))

  # Clusters
  Idents(seurat_obj) <- "seurat_clusters"
  clusters <- as.factor(Idents(seurat_obj))

  run_ssgsea_and_plot <- function(gene_sets, tag) {
    
    show <- ifelse(tag == "ALL", F, T)

    message("Running ssGSEA for ", tag, " pathways...")
    ssgsea_scores <- gsva(ssgseaParam(expr_mat, gene_sets), verbose = TRUE)

    # Palette
    vals <- as.vector(ssgsea_scores)
    col_fun <- colorRamp2(
      quantile(vals, c(0.01, 0.50, 0.99), na.rm=TRUE),
      c("#313695", "#ffffbf", "#a50026")
    )

    png_file_cells <- file.path(output_dir, paste0(prefix, "_", tag, "_per_cell.png"))
    png(png_file_cells, width=3300, height=2000, res=300)
    col_order <- order(clusters)
    ht1 <- Heatmap(
      ssgsea_scores[, col_order, drop=FALSE],
      name="ssGSEA score",
      col=col_fun,
      show_column_names=FALSE,
      show_row_names = show,
      row_names_gp=gpar(fontsize=9),
      column_split=clusters[col_order],
      cluster_rows=TRUE,
      cluster_columns=TRUE,
      column_title=paste0(tag, " pathways ssGSEA per cell"),
      column_title_gp=gpar(fontsize=14, fontface="bold"),
      heatmap_legend_param=list(title="score")
    )
    draw(ht1, heatmap_legend_side="right")
    dev.off()
    message("Saved: ", png_file_cells)

    # Moyenne par cluster
    unique_clusters <- levels(clusters)
    ssgsea_by_cluster <- sapply(unique_clusters, function(cl) {
      cols <- which(clusters == cl)
      if (length(cols) == 1) ssgsea_scores[, cols]
      else rowMeans(ssgsea_scores[, cols, drop=FALSE], na.rm=TRUE)
    })

    vals2 <- as.vector(ssgsea_by_cluster)
    col_fun2 <- colorRamp2(
      quantile(vals2, c(0.01, 0.50, 0.99), na.rm=TRUE),
      c("#313695", "#ffffbf", "#a50026")
    )

    png_file_cluster <- file.path(output_dir, paste0(prefix, "_", tag, "_per_cluster.png"))
    png(png_file_cluster, width=3300, height=2000, res=300)
    ht2 <- Heatmap(
      ssgsea_by_cluster,
      name="ssGSEA score",
      col=col_fun2,
      show_column_names=TRUE,
      show_row_names = show,
      column_names_gp=gpar(fontsize=12, fontface="bold", rot=45),
      row_names_gp=gpar(fontsize=9),
      cluster_rows=TRUE,
      cluster_columns=TRUE,
      column_title=paste0(tag, " pathways ssGSEA per cluster"),
      column_title_gp=gpar(fontsize=14, fontface="bold")
    )
    draw(ht2, heatmap_legend_side="right")
    dev.off()
    message("Saved: ", png_file_cluster)
  }

  run_ssgsea_and_plot(all_pathways, "ALL")
  run_ssgsea_and_plot(tf_pathways, "TF")

  message("Analyse ssGSEA terminée avec succès.")
}

object_list <- list(qiu_roth_cd8, yost_cd8, mathewson_smartseq2_cd8, mathewson_10X_cd8, lui_cd8, caushi_cd8, hu_cd8,
                    sade_feldman_cd8)
names(object_list) <- c("qiu_roth_cd8", "yost_cd8", "mathewson_smartseq2_cd8", "mathewson_10X_cd8", "lui_cd8",
                        "caushi_cd8", "hu_cd8", "sade_feldman_cd8")

lapply(object_list, ssgsea_seurat_pipeline)
```

# Different plots

``` r
# DotPlot of genes of interest

genes_of_interest <-  list(CXCR6posPD1neg = c("RUNX2", "RUNX3", "TOX", "TIGIT", "LAG3", "TOX2", "CXCR6", "PDCD1", "ENTPD1",  "CD101",  "HAVCR2",  "KLRB1", "CD74", "CD226",  "ID2",  "HIC1", "CD244", "SLAMF7",  "CREM"),
CXCR6negPD1pos = c("PRDM1", "IL7R", "ZNF683", "TCF7", "SELL","BACH2", "CCR7", "TNFRSF9", "CD28", "EOMES", "ICOS", "ID3",  "CD27", "CD28", "NR4A3", "NR4A1", "NR4A2", "ITGA1")
)

files <- list.files("CD8 subsets", pattern = "\\.rds$", full.names = TRUE)
clean_names <- sub("\\.rds$", "", basename(files))

for(name in clean_names){
  seu <- get(name)
  gg <- DotPlot(seu, features = lapply(genes_of_interest, unique), group.by = "seurat_clusters", cluster.idents = T) +
    scale_colour_viridis_c(option = "plasma") + 
    ggmin::theme_min() + 
    theme(
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(margin = margin(t = 15), size = 16),
        axis.title.y = element_text(margin = margin(r = 15), size = 16),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16, face = "bold")
    )
  ggsave(paste0("Plots/",name, "_DotPlot_gene_of_interest.png"), gg, width = 15, height = 10)
}

# CXCR6 and PD1 densities

datasets <- list(
  GBM_Smart_Seq2 = mathewson_smartseq2_cd8,
  NSCLC_1 = hu_cd8,
  PDAC = qiu_roth_cd8,
  GBM_10X = mathewson_10X_cd8,
  NSCLC_2 = lui_cd8,
  Melanoma = sade_feldman_cd8,
  SCC_BCC = yost_cd8,
  CRC_1 = che_cd8,
  CRC_2 = chen_cd8,
  CRC_3 = luiy_cd8,
  NSCLC_3 = caushi_cd8
)

extract_expr <- function(seu_obj, dataset_name, genes = c("CXCR6", "PDCD1")) {
  df <- FetchData(seu_obj, vars = genes, layer = "data")
  df$dataset <- dataset_name
  return(df)
}

expr_list <- purrr::map2(datasets, names(datasets), extract_expr)
expr_all <- bind_rows(expr_list)

expr_long_density <- expr_all %>%
  pivot_longer(cols = c("CXCR6", "PDCD1"),
               names_to = "Gene",
               values_to = "expr") %>%
  filter(expr > 0)

expr_long_prop <- expr_all %>%
  pivot_longer(cols = c("CXCR6", "PDCD1"),
               names_to = "Gene",
               values_to = "expr")

prop_df <- expr_long_prop %>%
  group_by(dataset, Gene) %>%
  summarise(
    prop_expressing = mean(expr > 0) * 100,
    n_cells = n(),
    .groups = "drop"
  )

p_density <- ggplot(expr_long_density, aes(x = expr, fill = Gene)) +
  geom_density(alpha = 0.5, adjust = 1.2) +
  facet_wrap(~dataset, scales = "free", ncol = 4) +
  scale_fill_viridis_d(option = "E") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Density of expressions of CXCR6 and PD1 by cancer type",
    x = "Expression",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

p_prop <- ggplot(prop_df, aes(x = dataset, y = prop_expressing, fill = Gene)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = paste0(round(prop_expressing, 1), "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 4.5) +
  scale_fill_viridis_d(option = "E") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Proportion of cells expressing CXCR6 or PD1",
    x = "Cancer type",
    y = "Proportion"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

final_plot <- p_density | p_prop

ggsave("Plots/CXCR6 and PD1 distributions.png", final_plot, width = 25, height = 15, dpi = 300)
```

# DE genes per cluster + Excel table

``` r
markers <- function(object, file_path){
  Idents(object) <- "seurat_clusters"
  object.markers <- FindAllMarkers(object, only.pos = TRUE)    
  object.markers <- object.markers %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.05)
  
  top_genes_per_cluster_object <- object.markers %>%    
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 200) %>%
    ungroup()
  
  wb <- createWorkbook()

  clusters <- unique(top_genes_per_cluster_object$cluster)
  for (clust in clusters) {
    cluster_data <- top_genes_per_cluster_object %>%
      filter(cluster == clust)
  
    addWorksheet(wb, paste0("Cluster_", clust))
  
    writeData(wb, sheet = paste0("Cluster_", clust), 
              x = data.frame(
                Gene = cluster_data$gene,
                avg_log2FC = cluster_data$avg_log2FC,
                p_val_adj = cluster_data$p_val_adj
              ))
  }
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

markers(qiu_roth_cd8, "200_DEG_Qiu_Roth_PDAC.xlsx")
markers(yost_cd8, "200_DEG_Yost_SCC_BCC.xlsx")
markers(caushi_cd8, "200_DEG_Caushi_NSCLC.xlsx")
markers(che_cd8, "200_DEG_Che_CRC.xlsx")
markers(chen_cd8, "200_DEG_Chen_CRC.xlsx")
markers(lui_cd8, "200_DEG_Lui_NSCLC.xlsx")
markers(hu_cd8, "200_DEG_Hu_NSCLC.xlsx")
markers(sade_feldman_cd8, "200_DEG_Sade_Feldman_Melanoma.xlsx")
markers(mathewson_smartseq2_cd8, "200_DEG_Mathewson_smartseq2_GBM.xlsx")
markers(mathewson_10X_cd8, "200_DEG_Mathewson_10X_GBM.xlsx")
markers(in_vitro_cd8, "200_DEG_In_Vitro.xlsx")
markers(luiy_cd8, "200_DEG_Luiy_CRC.xlsx")
```

# Pseudotime: DRAFT!!!

``` r
hu_cd8 <- readRDS("CD8 subsets/hu_cd8.rds")
gene_annotation <- as.data.frame(rownames(hu_cd8@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(hu_cd8@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

New_matrix <- as.matrix(GetAssayData(hu_cd8, layer = "data"))
New_matrix <- New_matrix[rownames(hu_cd8@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = hu_cd8[[]],
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


Idents(hu_cd8) <- "seurat_clusters"
list_cluster <- hu_cd8@active.ident
names(list_cluster) <- colnames(New_matrix)

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-hu_cd8@reductions[["Harmony.umap"]]@cell.embeddings


cds_from_seurat@principal_graph_aux$gene_loadings <- hu_cd8@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

plot_cells(cds_from_seurat, 
                   color_cells_by = 'cluster',
                   label_groups_by_cluster=TRUE,
                   label_leaves=FALSE,
                   label_branch_points=TRUE,
           graph_label_size=4)

cds_from_seurat <- order_cells(cds_from_seurat)

plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```

# Other requests

## BoxPlot

``` r
plot_seurat_box <- function(obj, patterns, char1, char2, output_path, char3) {
  obj$test <- paste(obj$seurat_clusters, obj[[char1]][, 1], sep = "_")

  tab <- table(obj[[]]$test, obj[[]][[char2]])
  tab <- as.data.frame.matrix(tab)

  tab <- sweep(tab, 2, colSums(tab), FUN = "/")  
  
  pattern_exact <- paste0("^(", paste(patterns, collapse = "|"), ")_")
  keep <- grepl(pattern_exact, rownames(tab))
  tab_filtered <- tab[keep, , drop = FALSE]

  df_long <- reshape2::melt(tab_filtered, id.vars = NULL, variable.name = "Condition", value.name = "Count")
  df_long$Condition <- factor(df_long$Condition, levels = colnames(tab_filtered))
  df_long$Count <- asin(sqrt(df_long$Count))

  stat_results <- compare_means(Count ~ Condition, data = df_long, method = "wilcox.test")

  stat_results$p.adj <- p.adjust(stat_results$p, method = "BH")

  stat_results$p.signif <- symnum(
    stat_results$p.adj,
    corr = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", "ns")
  )
  
  y_max <- max(df_long$Count, na.rm = TRUE)
  step <- (y_max * 0.1)
  stat_results <- stat_results %>%
    mutate(y.position = y_max + step * (1:n()))

  n_cols <- length(unique(df_long$Condition))
  if (n_cols == 2) {
    colors <- c("#FFB3BA", "#BAE1FF") 
  } else {
    colors <- scales::hue_pal(l = 80, c = 30)(n_cols)
  }

  p <- ggplot(df_long, aes(x = Condition, y = Count, fill = Condition)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "black", outlier.size = 2, alpha = 0.8) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste0(char2 , " in ", char3),
      x = char2,
      y = "Arcsin(sqrt(Proportion))"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    stat_pvalue_manual(
      stat_results,
      label = "p.signif",
      tip.length = 0.01,
      step.increase = 0.08,
      size = 5
    )

  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300)
}

plot_seurat_box(hu_cd8, c("2", "5"), "Patient", "Pathologic_Response", "Plots/box_plot_hu_cd8_Response_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(hu_cd8, c("2", "5"), "Patient", "Resource", "Plots/box_plot_hu_cd8_Treatment_Pos_Neg.png", "CXCR6posPD1neg")
# plot_seurat_box(yost_cd8, c("2", "5"), "Patient", "Pathologic_Response", "test.png")  pas de clusters d'interêt
plot_seurat_box(sade_feldman_cd8, c("3", "5"), "Patient ID", "Response", "Plots/box_plot_sade_feldman_cd8_Response_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(sade_feldman_cd8, c("3", "5"), "Patient ID", "TimePoint", "Plots/box_plot_sade_feldman_cd8_TimePoint_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(lui_cd8, c("9"), "patient", "TimePoint", "Plots/box_plot_lui_cd8_TimePoint_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(lui_cd8, c("10", "1", "13"), "patient", "TimePoint", "Plots/box_plot_lui_cd8_TimePoint_Neg_Pos.png", "CXCR6negPD1pos")
plot_seurat_box(lui_cd8, c("9"), "patient", "Response", "Plots/box_plot_lui_cd8_Response_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(lui_cd8, c("10", "1", "13"), "patient", "Response", "Plots/box_plot_lui_cd8_Response_Neg_Pos.png", "CXCR6negPD1pos")
plot_seurat_box(luiy_cd8, c("21"), "patient", "tissue", "Plots/box_plot_luiy_cd8_Tissue_Neg_Pos.png", "CXCR6negPD1pos")
plot_seurat_box(che_cd8, c("17"), "patient", "tissue", "Plots/box_plot_che_cd8_Tissue_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(che_cd8, c("19"), "patient", "tissue", "Plots/box_plot_che_cd8_Tissue_Neg_Pos.png", "CXCR6negPD1pos")
# plot_seurat_box(chen_cd8, c("2", "5"), "Patient", "Pathologic_Response", "test.png") pas de clusters d'interêt
plot_seurat_box(mathewson_10X_cd8, c("3"), "sampleid", "ImmunoTherapy", "Plots/box_plot_mathewson_10X_cd8_Therapy_Pos_Neg.png", "CXCR6posPD1neg")
# plot_seurat_box(mathewson_smartseq2_cd8, c("2", "5"), "Patient", "Pathologic_Response", "test.png") pas de clusters d'interêt
# plot_seurat_box(in_vitro_cd8, c("3"), "Patient", "Pathologic_Response", "in_vitro_cd8_Therapy_Pos_Neg.png") pas de var
plot_seurat_box(qiu_roth_cd8, c("6"), "orig.ident", "DiseaseState", "Plots/box_plot_qiu_roth_cd8_DiseaseState_Pos_Neg.png", "CXCR6posPD1neg")
plot_seurat_box(qiu_roth_cd8, c("9"), "orig.ident", "DiseaseState", "Plots/box_plot_qiu_roth_cd8_DiseaseState_Neg_Pos.png", "CXCR6negPD1pos")
plot_seurat_box(caushi_cd8, c("26"), "patient", "tissue", "Plots/box_plot_caushi_cd8_Tissue_Neg_Pos.png", "CXCR6negPD1pos")
plot_seurat_box(caushi_cd8, c("26"), "patient", "MPR status", "Plots/box_plot_caushi_cd8_Response_Neg_Pos.png", "CXCR6negPD1pos")
```

## BarPlot

``` r
plot_group_bar <- function(obj, clusters_vec, char1, output_path, char2) {
  meta <- obj[[]]
  meta$seurat_clusters <- as.character(meta$seurat_clusters)

  # Filtrer uniquement les clusters sélectionnés
  sub_meta <- meta %>% filter(seurat_clusters %in% clusters_vec)

  # Comptage par catégorie dans le sous-groupe
  df_plot <- sub_meta %>%
    group_by(!!sym(char1)) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      pct = n / nrow(meta) * 100  # proportion par rapport à toutes les cellules
    )

  # Palette de couleurs pastel
  pastel_palette <- c(
    "#FFB3BA", "#BAE1FF", "#BFFCC6", "#FFFFBA", "#FFDFBA", "#D5AAFF", "#C7CEEA"
  )

  # Générer le barplot
  p <- ggplot(df_plot, aes(x = !!sym(char1), y = pct, fill = !!sym(char1))) +
    geom_bar(stat = "identity", width = 0.7, color = "white") +
    scale_fill_manual(values = pastel_palette[seq_len(nrow(df_plot))]) +
    labs(
      title = paste0(char1 , " in ", char2),
      x = NULL,         
      y = "Percentage of total cells"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  # Sauvegarder le plot
  ggsave(output_path, plot = p, width = 7, height = 5, dpi = 300)
}


plot_group_bar(hu_cd8, c("2", "5"), "Pathologic_Response", "Plots/bar_plot_hu_cd8_Response_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(hu_cd8, c("2", "5"), "Resource", "Plots/bar_plot_hu_cd8_Treatment_Pos_Neg.png", "CXCR6posPD1neg")
# plot_group_bar(yost_cd8, c("2", "5"), "Pathologic_Response", "test.png")  pas de clusters d'interêt
plot_group_bar(sade_feldman_cd8, c("3", "5"), "Response", "Plots/bar_plot_sade_feldman_cd8_Response_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(sade_feldman_cd8, c("3", "5"), "TimePoint", "Plots/bar_plot_sade_feldman_cd8_TimePoint_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(lui_cd8, c("9"), "TimePoint", "Plots/bar_plot_lui_cd8_TimePoint_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(lui_cd8, c("10", "1", "13"), "TimePoint", "Plots/bar_plot_lui_cd8_TimePoint_Neg_Pos.png", "CXCR6negPD1pos")
plot_group_bar(lui_cd8, c("9"), "Response", "Plots/bar_plot_lui_cd8_Response_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(lui_cd8, c("10", "1", "13"), "Response", "Plots/bar_plot_lui_cd8_Response_Neg_Pos.png", "CXCR6negPD1pos")
plot_group_bar(luiy_cd8, c("21"), "tissue", "Plots/bar_plot_luiy_cd8_Tissue_Neg_Pos.png", "CXCR6negPD1pos")
plot_group_bar(che_cd8, c("17"), "tissue", "Plots/bar_plot_che_cd8_Tissue_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(che_cd8, c("19"), "tissue", "Plots/bar_plot_che_cd8_Tissue_Neg_Pos.png", "CXCR6negPD1pos")
# plot_group_bar(chen_cd8, c("2", "5"), "Pathologic_Response", "test.png") pas de clusters d'interêt
plot_group_bar(mathewson_10X_cd8, c("3"), "ImmunoTherapy", "Plots/bar_plot_mathewson_10X_cd8_Therapy_Pos_Neg.png", "CXCR6posPD1neg")
# plot_group_bar(mathewson_smartseq2_cd8, c("2", "5"), "Pathologic_Response", "test.png") pas de clusters d'interêt
# plot_group_bar(in_vitro_cd8, c("3"), "Pathologic_Response", "in_vitro_cd8_Therapy_Pos_Neg.png") pas de var
plot_group_bar(qiu_roth_cd8, c("6"), "DiseaseState", "Plots/bar_plot_qiu_roth_cd8_DiseaseState_Pos_Neg.png", "CXCR6posPD1neg")
plot_group_bar(qiu_roth_cd8, c("9"), "DiseaseState", "Plots/bar_plot_qiu_roth_cd8_DiseaseState_Neg_Pos.png", "CXCR6negPD1pos")
plot_group_bar(caushi_cd8, c("26"), "tissue", "Plots/bar_plot_caushi_cd8_Tissue_Neg_Pos.png", "CXCR6negPD1pos")
plot_group_bar(caushi_cd8, c("26"), "MPR status", "Plots/bar_plot_caushi_cd8_Response_Neg_Pos.png", "CXCR6negPD1pos")
```

## ssGSEA Comparison

``` r
calculate_ssgsea_comparison <- function(obj, clusters_vec1, clusters_vec2, group_name1, group_name2, output_path) {
  vars <- c("Caushi, Smith, Nature 2021",
            "Oliveira, Wu, Nature 2021",
            "Lowery, F.J,  Science 2022",
            "Hanada, K.i, Cancer Cell 2022")

  ssgsea_scores <- t(obj[[vars]])

  meta <- obj[[]]
  meta$seurat_clusters <- as.character(meta$seurat_clusters)

  meta$TRM_flag <- meta$TRM == "TRM"
  
  if (!is.null(clusters_vec1) && !is.null(clusters_vec2)) {
    meta$group <- ifelse(meta$seurat_clusters %in% clusters_vec1, group_name1,
                         ifelse(meta$seurat_clusters %in% clusters_vec2, group_name2, NA))
    
  } else if (!is.null(clusters_vec1) && is.null(clusters_vec2)) {
    meta$group <- ifelse(meta$seurat_clusters %in% clusters_vec1, group_name1,
                         ifelse(meta$TRM_flag & !(meta$seurat_clusters %in% clusters_vec1), group_name2, NA))
    
  } else if (is.null(clusters_vec1) && !is.null(clusters_vec2)) {
    meta$group <- ifelse(meta$seurat_clusters %in% clusters_vec2, group_name2,
                         ifelse(meta$TRM_flag & !(meta$seurat_clusters %in% clusters_vec2), group_name1, NA))
    
  }
  
  keep_cells <- !is.na(meta$group)
  meta <- meta[keep_cells, , drop = FALSE]
  ssgsea_scores <- ssgsea_scores[, keep_cells, drop = FALSE]

  df_long <- as.data.frame(t(ssgsea_scores))
  df_long$Group <- meta$group
  df_long <- reshape2::melt(df_long, id.vars = "Group", variable.name = "Pathway", value.name = "Score")

  stat_results <- df_long %>%
    group_by(Pathway) %>%
    summarise(
      p = tryCatch(
        wilcox.test(Score[Group == group_name1],
                    Score[Group == group_name2])$p.value,
        error = function(e) NA
      ),
      med1 = median(Score[Group == group_name1], na.rm = TRUE),
      med2 = median(Score[Group == group_name2], na.rm = TRUE)
    ) %>%
    mutate(
      medianDiff = med1 - med2,
      p.adj = p.adjust(p, method = "BH"),
      p.signif = symnum(
        p.adj,
        corr = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("***", "**", "*", ".", "ns")
      ),
      group1 = group_name1,
      group2 = group_name2
    )

  colors <- c("#FFB3BA", "#BAE1FF")
  names(colors) <- c(group_name1, group_name2)

  y_positions <- df_long %>%
    group_by(Pathway) %>%
    summarise(y.position = max(Score, na.rm = TRUE) * 1.1)
  
  stat_results_filtered <- stat_results %>%
    left_join(y_positions, by = "Pathway")
  
  p <- ggplot(df_long, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "black", outlier.size = 1.8,
                 alpha = 0.85, width = 0.6) +
    scale_fill_manual(values = colors) +
    facet_wrap(~ Pathway, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = NULL, x = NULL, y = "ssGSEA score") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey80"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11)
    ) +
    ggpubr::stat_pvalue_manual(
      stat_results_filtered,
      label = "p.signif",
      tip.length = 0.005,
      size = 4.5
    )

  ggsave(output_path, plot = p, width = 12, height = 6, dpi = 300)
  
  return(stat_results)
}



res_hu_cd8 <- calculate_ssgsea_comparison(hu_cd8, c("2", "5"), NULL, "CXCR6posPD1neg", "Other TRM", "Plots/ssgsea_hu_cd8_Pos_Neg_vs_other.png")
# res_ <- calculate_ssgsea_comparison(yost_cd8, c("2", "5"), "Pathologic_Response", "test.png")  pas de clusters d'interêt
res_sade_feldman_cd8 <- calculate_ssgsea_comparison(sade_feldman_cd8, c("3", "5"), NULL, "CXCR6posPD1neg", "Other TRM", "Plots/ssgsea_sade_feldman_cd8_Pos_Neg_vs_other.png")
res_lui_cd8 <- calculate_ssgsea_comparison(lui_cd8, c("10", "1", "13"), c("9"), "CXCR6negPD1pos", "CXCR6posPD1neg", "Plots/ssgsea_lui_cd8_Neg_Pos_vs_Pos_Neg.png")
res_luiy_cd8 <- calculate_ssgsea_comparison(luiy_cd8, c("21"), NULL, "CXCR6negPD1pos", "Other TRM", "Plots/ssgsea_luiy_cd8_Neg_Pos_vs_other.png")
res_che_cd8 <- calculate_ssgsea_comparison(che_cd8, c("19"), c("17"), "CXCR6negPD1pos", "CXCR6posPD1neg", "Plots/ssgsea_che_cd8_Neg_Pos_vs_Pos_Neg.png")
# res_ <- calculate_ssgsea_comparison(chen_cd8, c("2", "5"), "Pathologic_Response", "test.png") pas de clusters d'interêt
res_mathewson_10X_cd8 <- calculate_ssgsea_comparison(mathewson_10X_cd8, c("3"), NULL, "CXCR6posPD1neg", "Other TRM", "Plots/ssgsea_mathewson_10X_cd8_Pos_Neg_vs_other.png")
# res_ <- calculate_ssgsea_comparison(mathewson_smartseq2_cd8, c("2", "5"), "Pathologic_Response", "test.png") pas de clusters d'interêt
# res_ <- calculate_ssgsea_comparison(in_vitro_cd8, c("3"), "Pathologic_Response", "in_vitro_cd8_Therapy_Pos_Neg.png") pas de var
res_qiu_roth_cd8 <- calculate_ssgsea_comparison(qiu_roth_cd8, c("9"),  c("6"), "CXCR6negPD1pos", "CXCR6posPD1neg", "Plots/ssgsea_qiu_roth_cd8_Neg_Pos_vs_Pos_Neg.png")
res_caushi_cd8 <- calculate_ssgsea_comparison(caushi_cd8, c("26"), NULL, "CXCR6negPD1pos", "Other TRM", "Plots/ssgsea_caushi_cd8_Neg_Pos_vs_other.png")
```

## Split.by UMAPS

``` r
plot_umap_grid_facets <- function(obj, char1 = NULL, char2 = NULL, red, output_path) {
  meta <- obj[[]]

  if (!(red %in% names(obj@reductions))) stop(paste0(red, " not found"))
  if (is.null(char1) & is.null(char2)) stop("At least one of char1 or char2 must be provided")
  if (!is.null(char1) && !(char1 %in% colnames(meta))) stop(paste0(char1, " not found"))
  if (!is.null(char2) && !(char2 %in% colnames(meta))) stop(paste0(char2, " not found"))
  
  umap_coords <- Embeddings(obj, red)
  xlim <- range(umap_coords[,1])
  ylim <- range(umap_coords[,2])
  
  if (!is.null(char1) && !is.null(char2)) {
    levels_char1 <- unique(meta[[char1]])
    levels_char2 <- unique(meta[[char2]])
  } 
  else if (!is.null(char1)) {
    levels_char1 <- unique(meta[[char1]])
    levels_char2 <- "All cells"
    meta$`All cells` <- "All cells"
    obj$`All cells` <- "All cells"
    char2 <- "All cells"
  } else if (!is.null(char2)) {
    levels_char1 <- "All cells"
    levels_char2 <- unique(meta[[char2]])
    meta$`All cells` <- "All cells"
    obj$`All cells` <- "All cells"
    char1 <- "All cells"
  }
  
  n_colors <- length(unique(meta[[char1]]))
  pastel_palette <- colorRampPalette(brewer.pal(min(8, n_colors), "Set3"))(n_colors)
  names(pastel_palette) <- unique(meta[[char1]])
  
  plots_list <- list()
  
  for (i in seq_along(levels_char1)) {      
    lvl1 <- levels_char1[i]
    for (j in seq_along(levels_char2)) {   
      lvl2 <- levels_char2[j]
      cells <- rownames(meta)[meta[[char1]] == lvl1 & meta[[char2]] == lvl2]
      
      if (length(cells) > 0) {
        sub_obj <- subset(obj, cells = cells)
        p <- DimPlot(sub_obj, reduction = red, group.by = char1, raster = FALSE) +
          scale_color_manual(values = pastel_palette) +
          labs(title = NULL) +
          coord_cartesian(xlim = xlim, ylim = ylim)
      } else {
        df_empty <- data.frame(UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2])[0, ]
        p <- ggplot(df_empty, aes(x = UMAP_1, y = UMAP_2)) +
          geom_blank() +
          coord_cartesian(xlim = xlim, ylim = ylim)
      }
      
      p <- p +
        ggmin::theme_min() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm")
        )
      
      if (i == 1) {
        p <- p + ggtitle(lvl2) + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      }
      
      if (j == 1) {
        p <- p + ylab(lvl1) + theme(axis.title.y = element_text(size = 14, face = "bold"))
      }
      
      plots_list[[paste0("line", i, "_col", j)]] <- p
    }
  }
  
  final_plot <- wrap_plots(plots_list, nrow = length(levels_char1), ncol = length(levels_char2), byrow = TRUE) +
    plot_annotation(title = paste("UMAPs splited by", char1, "and", char2, sep = " "),
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  ggsave(output_path, plot = final_plot,
         width = 4*length(levels_char2), height = 4*length(levels_char1), dpi = 300)
}

plot_umap_grid_facets(hu_cd8, "Resource", "Pathologic_Response", "Harmony.umap", "Plots/hu_cd8_UMAP_split.png")
plot_umap_grid_facets(sade_feldman_cd8, "Response", "TimePoint", "umap", "Plots/sade_feldman_cd8_UMAP_split.png")
plot_umap_grid_facets(lui_cd8, "Response", "TimePoint", "Harmony.umap", "Plots/lui_cd8_UMAP_split.png")
plot_umap_grid_facets(luiy_cd8, NULL, "tissue", "umap", "Plots/luiy_cd8_UMAP_split.png")
plot_umap_grid_facets(mathewson_10X_cd8, NULL, "ImmunoTherapy", "RPCA.umap", "Plots/mathewson_10X_cd8_UMAP_split.png")
plot_umap_grid_facets(mathewson_smartseq2_cd8, "Cohort", "therapy", "umap", "Plots/mathewson_smartseq2_cd8_UMAP_split.png")
plot_umap_grid_facets(che_cd8, NULL, "tissue", "umap", "Plots/che_cd8_UMAP_split.png")
plot_umap_grid_facets(chen_cd8, "Response", "Treatment Stage", "bbknn_umap", "Plots/chen_cd8_UMAP_split.png")
plot_umap_grid_facets(caushi_cd8, "MPR status", "tissue", "umap", "Plots/caushi_cd8_cd8_UMAP_split.png")
plot_umap_grid_facets(yost_cd8, "Response", "treatment", "RPCA.umap", "Plots/yost_cd8_UMAP_split.png")
plot_umap_grid_facets(qiu_roth_cd8, NULL, "DiseaseState", "Harmony.umap", "Plots/qiu_roth_cd8_UMAP_split.png")
```

# Test : BPCells

``` r
library(Azimuth)  
library(BPCells)

file.path <- "Wells_healthy_donnors/7c98dadf-b249-4873-83be-51776674e9e1.h5ad"

# Ouvrir la matrice depuis le h5ad
data <- open_matrix_anndata_hdf5(file.path)

# Écrire les matrices BP sur disque
output.dir <- paste0(gsub(".h5ad", "", file.path), "_BP")
write_matrix_dir(
  mat = data,
  dir = output.dir
)

# Charger la matrice BP depuis le disque
mat <- open_matrix_dir(dir = output.dir)

# Convertir des IDs de gènes Ensembl en symboles (humain)
mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")

# Charger les métadonnées du fichier h5ad
metadata <- LoadH5ADobs(path = file.path)
```

# Test : scPred

``` r
reference <- scPred::pbmc_1
query <- scPred::pbmc_2

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

load("dds_list.RData")

vst <- vst(dds_list[[1]], blind = T)

coldata <- dds_list[[1]]$subpopulation

sce_ref <- SingleCellExperiment(
  assays = list(
    expression = assay(vst)   
  ),
  colData = DataFrame(
    cell_type = coldata    
  )
)


res <- CHETAHclassifier(sce, ref_cells = sce_ref, ref_ct = "cell_type", ref_c = "expression", input_c = "logcounts")

res_seu <- as.Seurat(res)

head(res_seu[[]])

load("metadata_cd8.RData")

res_seu[[]]$clustify <- metadata_cd8$hu_cd8$TRM_SubPop_type

DimPlot(res_seu, reduction = "HARMONY.UMAP", group.by = c("clustify", "celltype_CHETAH"))
```
