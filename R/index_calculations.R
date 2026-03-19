


#' General function for index calculation based on a known dataset (i.e. microglia development index)
#'
#' @param mat A gene expression matrix
#' @param index_df Dataframe with columns including gene names, valence, and directionality
#' @param gene_column Name of the column in the dataframe that contains the gene names for indexing
#' @param direction_column Dataframe with columns including gene names, valence, and directionality

#' @return Dataframe with index score for each sample (row from initial matrix)
#' @export
#'

calculate_index <- function(mat, index_df, gene_column = 'gene', direction_column = 'direction') {

  up_genes   <- stringr::str_to_title(index_df[[gene_column]][index_df[[direction_column]] == "UP"])
  down_genes <- stringr::str_to_title(index_df[[gene_column]][index_df[[direction_column]] == "DOWN"])

  # intersecting the two datasets to make sure we only have the genes that are present in the receiving dataset
  up_genes_present   <- intersect(up_genes, rownames(mat))
  down_genes_present <- intersect(down_genes, rownames(mat))

  up_mat   <- mat[up_genes_present, , drop = FALSE]
  down_mat <- mat[down_genes_present, , drop = FALSE]

  up_mean   <- colMeans(up_mat)
  down_mean <- colMeans(down_mat)

  idx_score <- up_mean / down_mean

  return(data.frame(index = idx_score))
}



#' Function wrapper for calculating index using general calculate_index function on seurat object
#'
#' @param seurat_obj A seurat object already created
#' @param index_df Dataframe with columns including gene names, valence, and directionality
#' @param gene_column Name of the column in the dataframe that contains the gene names for indexing
#' @param direction_column Dataframe with columns including gene names, valence, and directionality
#' @param column_name What to name the output metadata column that is contained within the seurat object
#'
#' @return A seurat obj that was the input object + a metadata column containing index values accessible by "column_name" argument in metadata
#' @export
#'

#
calculate_index_seurat <- function(seurat_obj, index_df, gene_column = 'gene', direction_column = 'direction', column_name = "developmental_index") {
  mat <- Seurat::GetAssayData(seurat_obj, layer = "data")
  seurat_obj[[column_name]] <- calculate_index(mat, index_df)
  return(seurat_obj)
}

#' Function for calculating gene variance and adding it as a metadata column to the seurat object
#'
#' @param s_obj A seurat object already created
#' @param assay What to name the output metadata column that is contained within the seurat object
#' @param slot to name the output metadata column that is contained within the seurat object
#' @param column_name What to name the output metadata column that is contained within the seurat object
#'
#' @return A seurat obj that was the input object + a metadata column containing index values accessible by "column_name" argument in metadata
#' @export
#'

#
seurat_gene_var <- function(s_obj, assay = "RNA", slot = "data", column_name = "gene_variance") {
  # extract expression matrix
  mat <- Seurat::GetAssayData(s_obj, layer = "data")

  # calculate per gene variance
  cell_var <- sparseMatrixStats::colVars(mat)

  # add to metadata
  s_obj[[column_name]] <- cell_var

  return(s_obj)
}



#' Function for calculating rolling variance averaged across a specified gene window size
#'
#' @param s_obj A seurat object already created
#' @param index_col Name of the metadata column in the seurat object that contains index calculation
#' @param variance_col Name of the metadata column in the seurat object that contains gene variance calculation
#' @param window_size Window size (number of genes) to bin
#'
#' @return A dataframe containing the rolling index and rolling variance for comparison
#' @export
#'

rolling_variance_index <- function(s_obj, index_col = "developmental_index", variance_col = "gene_variance", window_size = 200) {

  meta <- s_obj@meta.data

  if (!(index_col %in% colnames(meta))) {
    stop("Index column not found in seurat metadata. Please calculate it with calculate_index_seurat(sobject, index_df)")
  }

  if (!(variance_col %in% colnames(meta))) {
    stop("Variance column not found in seurat metadata. Please calculate it with the function total_gene_var(sobject)")
  }

  df <- data.frame(index = meta[[index_col]],variance = meta[[variance_col]])
  df <- df[order(df$index), ]
  n <- nrow(df)

  rolling_index <- numeric(n - window_size + 1)
  rolling_var <- numeric(n - window_size + 1)

  for (i in seq_len(n - window_size + 1)) {

    window <- df[i:(i + window_size - 1), ]

    rolling_index[i] <- mean(window$index)
    rolling_var[i] <- mean(window$variance)

  }

  data.frame(index = rolling_index, variance = rolling_var)
}


#' Function to find genes that are involved in transitions between cells across index values
#'
#' @param s_obj A seurat object already created
#' @param index_name Which index to use for the calculations, contained in seurat object
#' @param chunk_size For chunking by x number of genes to do calculation without overloading memory, default of 500 genes works just fine
#' @param min_pct Minimum percent of cells that express a given gene. Default is 5%, if a gene is expressed in less that 5% of cells, it is dropped
#'
#' @return A dataframe containing the correlation between each gene and overall developmental index
#' @export
#'

find_transition_genes <- function(s_obj, index_name = "developmental_index", chunk_size = 500, min_pct = 0.05) {
  expr <- Seurat::GetAssayData(s_obj, layer = "data")
  pct_expr <- Matrix::rowMeans(expr > 0)
  expr <- expr[pct_expr >= min_pct, ]
  meta <- s_obj@meta.data

  dev_index <- meta[[index_name]]

  ri <- rank(dev_index)
  ri_c <- ri - mean(ri)
  denom_i <- sqrt(sum(ri_c^2))

  n_genes <- nrow(expr)
  cors <- numeric(n_genes)
  n_chunks <- ceiling(n_genes / chunk_size)

  for (chunk in seq_len(n_chunks)) {
    i <- (chunk - 1) * chunk_size + 1
    idx <- i:min(i + chunk_size - 1, n_genes)

    message(sprintf("Processing chunk %d / %d (genes %d-%d)",
                    chunk, n_chunks, min(idx), max(idx)))

    re <- matrixStats::rowRanks(as.matrix(expr[idx, ]), ties.method = "average")
    re_c <- re - rowMeans(re)
    cors[idx] <- (re_c %*% ri_c) / (sqrt(rowSums(re_c^2)) * denom_i)
  }

  names(cors) <- rownames(expr)
  data.frame(cors)
}


#' Function wrapper for calculating index using general calculate_index function on deseq2 object
#'
#' @param dds A deseq2 dataset object already created
#' @param index_df Dataframe with columns including gene names, valence, and directionality
#' @param gene_column Name of the column in the dataframe that contains the gene names for indexing
#' @param direction_column Dataframe with columns including gene names, valence, and directionality
#' @param use_vst Conditional to use variance stabilizing transformation data or not. If not, it will default to calculating index on the raw gene counts
#' @param column_name What to name the output metadata column that is contained within the seurat object
#'
#' @return A deseq2 dataset obj that was the input object + a metadata column containing index values accessible by "column_name" argument in metadata
#' @export



calculate_index_deseq <- function(dds, index_df, gene_column = 'gene', direction_column = 'direction', use_vst = TRUE, column_name = "developmental_index") {
  if (use_vst) {
    mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  } else {
    mat <- DESeq2::counts(dds, normalized = TRUE)
  }

  dds[[column_name]] <- calculate_index(mat, index_df)
  return(dds)
}


#' Function for calculating TPM on any gene matrix
#'
#' @param counts_mat Matrix of gene counts
#' @param length_vector Labelled character vector containing gene names and their lengths
#' @param log Whether to log transform the output or not (I usually don't)
#'
#' @return A matrix of TPM values
#' @export
#'

calculate_tpm.matrix <- function(counts_mat, length_vector, log = FALSE) {

  if (is.null(names(length_vector))) {
    stop("length_vector must be named by gene IDs.")
  }

  gene_lengths <- length_vector[rownames(counts_mat)]

  if (any(is.na(gene_lengths))) {
    warning(paste("Removing", sum(is.na(gene_lengths)), "genes out of" , length(gene_lengths)," with missing lengths."))
    keep <- !is.na(gene_lengths)
    counts_mat <- counts_mat[keep, , drop = FALSE]
    gene_lengths <- gene_lengths[keep]
  }

  if (any(gene_lengths <= 0)) {
    stop("Gene lengths must be positive.")
  }

  gene_lengths_kb <- gene_lengths / 1000

  # RPK
  rpk <- counts_mat / gene_lengths_kb
  # Per-sample scaling
  scaling_factors <- colSums(rpk)
  tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6

  if (log == TRUE) {
    tpm <- log2(tpm + 1)
  }

  return(tpm)
}

#' Function for calculating TPM on any gene matrix
#'
#' @param object DESeq2 object
#' @param length_vector Labelled character vector containing gene names and their lengths
#' @param assay_name What to name the resulting assay that is saved in the deseq2 object
#' @param log Whether to log transform the output or not (I usually don't)
#' @param add_to_object Conditional for adding to object. You should almost always do it
#'
#' @return A matrix of TPM values
#' @export
#'

calculate_tpm.DESeqDataSet <- function(object, length_vector, assay_name = "TPM", log = FALSE, add_to_object = TRUE) {

  counts_mat <- DESeq2::counts(object, normalized = FALSE)

  tpm <- calculate_tpm.matrix(counts_mat, length_vector = length_vector, log = log)

  if (add_to_object) {
    SummarizedExperiment::assay(object, assay_name) <- tpm
    return(object)
  } else {
    return(tpm)
  }
}
