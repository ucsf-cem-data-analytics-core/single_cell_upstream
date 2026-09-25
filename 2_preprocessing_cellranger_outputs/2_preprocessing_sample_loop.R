#!/usr/bin/env Rscript

# START ----
# Workflow for creating analysis objects

library(Seurat)
library(rjson)
library(tidyverse)
library(stringr)
library(knitr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(kableExtra)
library(devtools)
library(parallel)
library(data.table)
library(harmony)
library(DoubletFinder)
library(impactSingleCellToolkit)

library(SoupX)
library(dsb)
library(DropletUtils)

# INPUT and OUTPUT Directories
#  - refer to template param file 
# 2_preprocessing_sample_loop.R  [Param File]

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("ERROR: This code needs 1 arguements: param file with [PATH copied cellranger counts directory] [hd5 files PATH] [PATH output directory]", call.=FALSE)
} 
param_file_fh = args[1]
##param_file_fh = "../input/input_full_cohort_analysis.json"
params        = fromJSON(file = param_file_fh)


# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
plot_dir      = params$plot_directory
THREADS       = params$threads

# Functions ----
# These functions could be incorporated into the R package

### runAmbientRNAfilterAndOuputFiles
runAmbientRNAfilterAndOuputFiles <- function(raw_hd5_fh,filt_hd5_fh,output_dir) {
  # Load and only focus on Gene Expression
  raw_mtx  <- Read10X_h5(raw_hd5_fh,use.names = T)
  raw_mtx  <- raw_mtx$`Gene Expression`
  filt_mtx <- Read10X_h5(filt_hd5_fh,use.names = T)
  filt_mtx <- filt_mtx$`Gene Expression`
  
  # make object
  seurat_obj   <- CreateSeuratObject(counts = filt_mtx)
  soup_channel <- SoupChannel(raw_mtx, filt_mtx)
  
  # process object
  seurat_obj    <- SCTransform(seurat_obj, verbose = F)
  seurat_obj    <- RunPCA(seurat_obj, verbose = F)
  seurat_obj    <- RunUMAP(seurat_obj, dims = 1:30, verbose = F)
  seurat_obj    <- FindNeighbors(seurat_obj, dims = 1:30, verbose = F)
  seurat_obj    <- FindClusters(seurat_obj, verbose = T)
  
  # Add clusters to channel
  meta_obj    <- seurat_obj@meta.data
  umap_obj    <- seurat_obj@reductions$umap@cell.embeddings
  soup_channel  <- setClusters(soup_channel, setNames(meta_obj$seurat_clusters, rownames(meta_obj)))
  soup_channel  <- setDR(soup_channel, umap_obj)
  
  # Calculate ambient RNA
  soup_channel  <- autoEstCont(soup_channel)
  # Convert to integers
  adj_matrix  <- adjustCounts(soup_channel, roundToInt = T)
  
  # Output
  DropletUtils:::write10xCounts(output_dir, adj_matrix, version = "3",overwrite = TRUE)
}


### preprocessCountsUsingMetadata ----
preprocessCountsUsingMetadata <- function(sample,metadata_counts,threads = 1,
                                          ambient.rna.filter = TRUE,output.type = "gex",
                                          min.genes.gex=400,multi.status=TRUE) {
  metadata_sample <- metadata_counts[metadata_counts$sample %in% sample,]
  counts_dir      <- metadata_sample$results_directory_path
  hd5_dir         <- metadata_sample$path_hd5
  counts_filt_dir <- metadata_sample$path_counts_filt
  
  # Completion status
  entry <- data.frame(matrix(ncol = 2, nrow = 1)) %>% setNames(c("sample","passed_preprocessing"))
  entry$sample               <- sample
  entry$passed_preprocessing <- "no"
  
  # Change multi status based on count type
  count_type      <- metadata_sample$type
  if(output.type == "gex"){
    if(count_type == "counts_gex"){multi.status <- FALSE}
  }
  
  # Directory to read in and apply filters and doublet removal
  input_dir <- counts_dir
  if(ambient.rna.filter){input_dir    <- counts_filt_dir}
  if(ambient.rna.filter){multi.status <- FALSE} # Post ambient RNA filter, the output is in a non-multi format
  
  # If performing ambient RNA filter output results to filtered counts folder
  gex.matrix <- NULL
  tryCatch({
    if(ambient.rna.filter){
      raw_mtx_fh  <- file.path(hd5_dir,"raw_feature_bc_matrix.h5")
      filt_mtx_fh <- file.path(hd5_dir,"sample_filtered_feature_bc_matrix.h5")
      
      runAmbientRNAfilterAndOuputFiles(raw_hd5_fh = raw_mtx_fh,
                                       filt_hd5_fh = filt_mtx_fh,
                                       output_dir = counts_filt_dir)
    }
    
    # Gene count filter on input sample matrix
    dataset_loc <- tstrsplit(input_dir,sample)[[1]]
    gex.matrix <- generateCombinedMatrix(dataset_loc, samples.vec = sample,THREADS = 1,
                                         multi.results = multi.status,assay = output.type,
                                         min.genes.per.cell = min.genes.gex,max.genes.per.cell = NULL)
  }, error = function(e) {
    gex.matrix <- NULL
  })
  
  # Only continue pre-processing if there are at least 100 cells
  if(length(colnames(gex.matrix)) > 100){
    # Create Seurat Object
    seurat.obj <- createSeuratObjectFromMatrix(
      sc.data      = gex.matrix,
      project.name = "GEX_one_sample",
      npca         = 20, min.genes = min.genes.gex,
      normalize = F,dim.reduction = F
    )
    
    # Add samples column
    cells      <- row.names(seurat.obj@meta.data)
    sample_col <- cells
    sample_col <- tstrsplit(sample_col,"-")[[2]]
    seurat.obj@meta.data$sample <- sample_col
    
    # Doublet removal
    seurat.obj <- findDoublets(seurat.obj,sample.col = "sample",threads = threads)
    
    # Omit Doublets 
    cells_to_keep  <- names(seurat.obj$doublet_finder[seurat.obj$doublet_finder %in% "Singlet"])
    seurat.obj     <- seurat.obj[,colnames(seurat.obj) %in% cells_to_keep]
    
    # Write output files
    counts_matrix_processed <- gex.matrix[,colnames(gex.matrix) %in% colnames(seurat.obj)]
    DropletUtils:::write10xCounts(input_dir, counts_matrix_processed, version = "3",overwrite = TRUE)
    
    # If it runs through all of this then it passed pre-processing
    entry$passed_preprocessing <- "yes"
  }
  
  
  return(entry)
}


# 1. list samples ----
metadata_counts <- metadata[metadata$type %in% c("counts","counts_gex"),]
metadata_counts$path_hd5 <- str_replace_all(metadata_counts$results_directory_path,"\\/multi_counts\\/","/multi_counts_hd5/")
metadata_counts$path_counts_filt <- str_replace_all(metadata_counts$results_directory_path,"\\/multi_counts\\/","/multi_counts_filt/")

samples <- metadata_counts$sample

# 2. Initialize variables used in nested functions ----
dataset_loc        <- ""
samples.vec        <- c()
multi.results      <- NULL
assay              <- NULL
min.genes.per.cell <- NULL
max.genes.per.cell <- NULL

# 2. loop through samples and pre-process ----
df_list <- lapply(as.list(samples), preprocessCountsUsingMetadata,
                  metadata_counts=metadata_counts,threads = THREADS,ambient.rna.filter=TRUE,
                  output.type="gex",min.genes.gex=400,multi.status=TRUE)

preprocessing_status_df <- do.call("rbind",df_list)

# Write a CSV of which samples passed pre-processing
write.csv(preprocessing_status_df,file = file.path(plot_dir,"results_preprocessing_status.csv"),row.names = F,quote = F)

