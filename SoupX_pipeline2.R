sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /home/relgamal/.conda/envs/reticulate/lib/libmkl_rt.so.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] SoupX_1.6.1               knitr_1.40               
 [3] forcats_0.5.1             stringr_1.4.1            
 [5] purrr_0.3.5               readr_2.1.2              
 [7] tidyr_1.2.1               tibble_3.1.8             
 [9] tidyverse_1.3.2           ggpubr_0.4.0             
[11] data.table_1.14.4         harmony_0.1.0            
[13] Rcpp_1.0.9                Matrix_1.5-1             
[15] ggplot2_3.4.0             dplyr_1.0.10             
[17] EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.18.3         
[19] AnnotationFilter_1.18.0   GenomicFeatures_1.46.4   
[21] AnnotationDbi_1.56.2      Biobase_2.54.0           
[23] Signac_1.7.0              sp_1.5-0                 
[25] SeuratObject_4.1.2        Seurat_4.2.0             
[27] hdf5r_1.3.5               GenomicRanges_1.46.1     
[29] GenomeInfoDb_1.30.1       IRanges_2.28.0           
[31] S4Vectors_0.32.4          BiocGenerics_0.40.0      

#######################################################################################################################################################

install.packages("SoupX")
library(SoupX)
library(Seurat)
library(dplyr)
library(tidyverse)
library(data.table)

samples <- c("HPAP-019","HPAP-020","HPAP-021","HPAP-022","HPAP-023","HPAP-024","HPAP-026","HPAP-027","HPAP-028","HPAP-029","HPAP-032","HPAP-034","HPAP-035","HPAP-036","HPAP-037","HPAP-038", "HPAP-039", "HPAP-040", "HPAP-042","HPAP-043","HPAP-044", "HPAP-045", "HPAP-047","HPAP-049", "HPAP-050","HPAP-052", "HPAP-053", "HPAP-054", "HPAP-055", "HPAP-056","HPAP-059","HPAP-063","HPAP-064","HPAP-071","HPAP-072","HPAP-074","HPAP-075","HPAP-077","HPAP-080","HPAP-082","HPAP-084","HPAP-087","HPAP-092","HPAP-093","HPAP-099","HPAP-101","HPAP-103","HPAP-104","HPAP-105", "HPAP-107")

for (sample in samples){
    wd <- sprintf('/path/to/data/here/Upenn_scRNAseq/cellranger_RME/%s/outs', samples)
    }

############## Automated Version of SoupX ##############

for (x in wd){
    sample_name <- str_split_fixed(x, "/", n=12)[9] #Adjust this to output your sample name
    rna_counts <- Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
    data <- CreateSeuratObject(counts=rna_counts)
    data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
    data <- subset(x = data, subset = nFeature_RNA > 500)
    data <- subset(x = data, subset = nFeature_RNA < 4000)
    
    #Running sctransform takes into account sequencing depth at each cell
    #data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes = FALSE)
    #data <- SCTransform(data, verbose = FALSE)
    
    #Log normalization alternative to sctransform
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    
    data <- RunPCA(data, verbose = FALSE)
    data <- RunUMAP(data, dims = 1:30, verbose = FALSE)
    data <- FindNeighbors(data, dims = 1:30, verbose = FALSE)
    data <- FindClusters(data, algorithm=4, resolution = 1, verbose=FALSE)
    
    #Read in RNA assay counts from our filtered seurat object
    DefaultAssay(data) <- 'RNA'
    toc <- GetAssayData(object = data, slot = "counts") #with nFeature >500 filter
    tod <- Seurat::Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5')) #raw count matrix
    
    #Pull out the required metadata from the clustered filtered adata object
    #We need the UMAP coordinates (RD1 and RD2) and the cluster assignments at minimum
    metadata <- (cbind(as.data.frame(data[["umap"]]@cell.embeddings),
                   as.data.frame(Idents(data)),
                   as.data.frame(Idents(data))))
    colnames(metadata) <- c("RD1","RD2","Cluster","Annotation")
    
    #Create the SoupChannel Object
    sc <- SoupChannel(tod,toc)
    
    #Add in the metadata (dimensionality reduction UMAP coords and cluster assignments)
    sc <- setDR(sc,metadata[colnames(sc$toc),c("RD1","RD2")])
    sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc)
    
    #Create post-SoupX Seurat Object
    data2 <- CreateSeuratObject(out)
    data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
    data2 <- NormalizeData(data2, normalization.method = "LogNormalize", scale.factor = 10000)  #Can be changed to sctransform
    data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data2)
    data2 <- ScaleData(data2, features = all.genes)
    data2 <- RunPCA(data2, verbose = FALSE)
    data2 <- RunUMAP(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindNeighbors(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindClusters(data2, algorithm=4, resolution = 1, verbose=FALSE)
    saveRDS(data2, file = sprintf("~/soupX_%s_automated.rds",sample_name))
    }

############## Gene Selection Version of SoupX ##############
#This method requires you to input marker genes that should be cluster/cell-type specific ie. INS in beta cells
 
contam_frac_results0 <- NULL
contam_frac_results <- data.frame() #contam_frac_results will be a dataframe with all samples and their contamination fraction estimates

for (x in wd){
    sample_name <- str_split_fixed(x, "/", n=12)[9] #Adjust this to output your sample name
    inputdata.10x <- Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
    rna_counts <- inputdata.10x
    data <- CreateSeuratObject(counts=rna_counts)
    data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
    data <- subset(x = data, subset = nFeature_RNA > 500)
    data <- subset(x = data, subset = nFeature_RNA < 4000)
    
    #Running sctransform takes into account sequencing depth at each cell
    #data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes = FALSE)
    #data <- SCTransform(data, verbose = FALSE)
    
    #Log normalization alternative to sctransform
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    
    data <- RunPCA(data, verbose = FALSE)
    data <- RunUMAP(data, dims = 1:30, verbose = FALSE)
    data <- FindNeighbors(data, dims = 1:30, verbose = FALSE)
    data <- FindClusters(data, algorithm=4, resolution = 1, verbose=FALSE)
    
    DefaultAssay(data) <- 'RNA'
    toc <- GetAssayData(object = data, slot = "counts") #with nFeature filter
    tod <- Seurat::Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5')) #raw count matrix
    
    metadata <- (cbind(as.data.frame(data[["umap"]]@cell.embeddings),
                   as.data.frame(Idents(data)),
                   as.data.frame(Idents(data))))
    colnames(metadata) <- c("RD1","RD2","Cluster","Annotation")
    
    sc <- SoupChannel(tod,toc)
    sc <- setDR(sc,metadata[colnames(sc$toc),c("RD1","RD2")])
    sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
    
    alphaGenes = c("GCG")
    betaGenes = c("INS")
    acinarGenes = c("REG1A")
    deltaGenes = c("SST")
    
    useToEst <- estimateNonExpressingCells(sc, nonExpressedGeneList = list(B=betaGenes, C=acinarGenes, A = alphaGenes, D= deltaGenes),)
    sc <- calculateContaminationFraction(sc,list(B=betaGenes, C=acinarGenes, A = alphaGenes, D = deltaGenes),useToEst=useToEst)
    contamination_fraction <- mean(sc$metaData$rho*100)
    
    contam_frac_results0$Sample <- sample_name
    contam_frac_results0$ContaminationFraction <- contamination_fraction
    contam_frac_results <- rbind(contam_frac_results, contam_frac_results0)
    
    out <- adjustCounts(sc)
    data2 <- CreateSeuratObject(out)
    data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
    data2 <- NormalizeData(data2, normalization.method = "LogNormalize", scale.factor = 10000) #Can be changed to sctransform
    data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data2)
    data2 <- ScaleData(data2, features = all.genes)
    data2 <- RunPCA(data2, verbose = FALSE)
    data2 <- RunUMAP(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindNeighbors(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindClusters(data2, algorithm=4, resolution = 1, verbose=FALSE)
    saveRDS(data2, file = sprintf("~/soupX_%s_gene_selection.rds",sample_name))
    }
    
############## Merging Samples to Single RDS Post-SoupX ##############
    
setwd("~/SoupX/")
soupx_files <- list.files("~/SoupX/", pattern="automated.rds")
soupx_data <- list()

for (x in soupx_files){
    sample_name <- str_split_fixed(x, "_", n=4)[2]  #adjust to output your sample name
    tmp <- readRDS(x)
    soupx_data[[sample_name]] <- tmp
}
   
soupx_merged_data <- merge(soupx_data[[samples[[1]]]], y=soupx_data[samples[2:length(samples)]], add.cell.ids=samples, project='HPAP')
soupx_merged_data$library <- substr(rownames(soupx_merged_data@meta.data),1,8)
soupx_merged_data
