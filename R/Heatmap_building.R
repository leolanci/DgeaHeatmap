#' Builds a matrix from an input csv file, with one column as rownames
#'
#' @param counts_data A data frame with floats or integers that
#' has strings in one column.
#' @param x Number of column which contains the rownames description.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
build_matrix <- function(counts_data, x) {
  list_rownames <- as.list(counts_data[, x])
  counts_data <- as.matrix(counts_data[, -x])
  rownames(counts_data) <- list_rownames
  return(counts_data)
}

#' Creates a matrix only containing chosen columns of an
#' original matrix with more data.
#'
#' @param factors_for_matrix_devision A list containing variables used to
#' extract specific columns from the original matrix.
#' @param mmatrix A matrix used as input, from which specifically chosen columns
#' are extracted.
#'
#' @return A new matrix containing fewer columns than before,
#' that have previously been specified.
#' @export
#'
#' @examples
#' factors_for_individual_matrix <- list("DKD", "glomerulus")
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' indiMatrix <- individual_matrix(factors_for_individual_matrix, matrixCounts)
#'
individual_matrix <- function(factors_for_matrix_devision, mmatrix) {
  for (i in factors_for_matrix_devision) {
    chosen_columns <- BiocGenerics::grep(i, BiocGenerics::colnames(mmatrix))
    rmatrix <- BiocGenerics::cbind(mmatrix[, chosen_columns])
    mmatrix <- rmatrix
    rmatrix <- matrix()
  }
  return(mmatrix)
}

#' Function to filter a matrix to extract a chosen number of most variable rows
#' through calculation of the variance.
#'
#' @param counts_data An input matrix from which the most variable rows are to
#' be extracted.
#' @param top_number_of_genes An integer to set the number of how many rows
#' should be extracted from the original matrix by variance.
#'
#' @return A new matrix containing only the top_number_of_genes count of rows
#' with the highest variance from input matrix.
#' @export
#'
#' @examples
#' x <- 1
#' top_number_of_genes <- 20
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' counts_data <- matrixCounts
#' varGenesMatrix <- filtering_for_top_exprGenes(counts_data,
#' top_number_of_genes)
filtering_for_top_exprGenes <- function(counts_data, top_number_of_genes) {
  # estimating the variance of each gene
  var_genes <- apply(counts_data, 1, stats::var)
  # Sorting the genes by their variance and creating a new object
  #with chosen number of most variable genes
  select_var <- names(sort(var_genes,
                           decreasing = TRUE))[seq_len(top_number_of_genes)]
  highly_variable_genes <- counts_data[select_var, ]
  dim(highly_variable_genes)

  return(highly_variable_genes)
}


#' Function to Z-count scale the values of a matrix.
#'
#' @param countsmatrix An input matrix whose values are scaled by Z-count.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A new matrix with Z-count scaled values.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' countsmatrix <- build_matrix(input_data, x)
#' scaled_counts <- scale_counts(countsmatrix)
scale_counts <- function(countsmatrix) {
  scaled_counts <-
    countsmatrix %>%
    t(.) %>%
    scale() %>%
    t(.)
  return(scaled_counts)
}

#' Function to visualize the data distribution within the data of a matrix.
#'
#' @importFrom rlang .data
#'
#' @param scaled_counts An input matrix with Z-count scaled values.
#'
#' @return A plot visualizing the data distribution of Z-count scaled data.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' scaled_counts <- scale_counts(matrixCounts)
#' show_data_distribution(scaled_counts)
show_data_distribution <- function(scaled_counts) {
  apply(scaled_counts, MARGIN = 1, mean) %>% # calculate the mean per row
    graphics::hist(., main = "", xlab = "Z-score values", col = "dodgerblue2")
  # build histogram to see data distribution
}



#' Function to create an elbow plot to choose k for clustering by k-Means.
#'
#' @param top_genes_matrix A matrix for which the best number of clusters in
#' k-Means Clustering is supposed to be calculated.
#' @param maxK An integer determining the maximum number of possible clusters.
#'
#' @return An elbow plot to choose k.
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' scaled_counts <- scale_counts(matrixCounts)
#' top_genes_matrix <- scaled_counts
#' set.seed(seed)
#' elbow_plot(top_genes_matrix)
elbow_plot <- function(top_genes_matrix, # input matrix.
                       maxK = 15 # maximum number of clusters
) {
  maxK <- maxK
  wss <- function(k) {
    stats::kmeans(top_genes_matrix,
                  algorithm = "Lloyd",
                  k, nstart = 10,
                  iter.max = 50)$tot.withinss
  }
  k.values <- seq_len(maxK) # compute and plot wws for k <- 1 to k <- 15
  wss_values <- purrr::map_dbl(k.values, wss)
  plot(k.values, wss_values, type = "b", pch = 19,
       frame = FALSE, xlab = "Number of clusters K",
       ylab = "Total within-clusters sum of squares")
}


#' Function to set color scheme for a heatmap.
#'
#' @param colorPalette Chosen color palette.
#'
#' @return No return.
#' @export
#'
#' @examples
#' colorPalette <- "RdBu"
#' color_setting(colorPalette)
color_setting <- function(colorPalette) {
  my_colors <- RColorBrewer::brewer.pal(n = 11, name = colorPalette)
  my_colors <- grDevices::colorRampPalette(my_colors)(50)
  my_colors <- rev(my_colors)
  my_colors
}

#' Summarizes columns biological replicates of a matrix into one.
#'
#' @param top_genes_matrix An input matrix with columns representing biological
#' replicates.
#' @param probes List of strings by which columns of biological replicates can
#' be identified.
#'
#' @return A new matrix where biological replicated are summarized
#' as their mean.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' top_genes_matrix <- build_matrix(input_data, x)
#' probes <- list("disease3", "disease4", "disease1B")
#' newMatrix <- summarise_bio_replicates(top_genes_matrix, probes)
summarise_bio_replicates <- function(top_genes_matrix, probes) {
  m_top_genes_matrix <- top_genes_matrix
  number_columns <- ncol(m_top_genes_matrix)
  # loop to grep columns for each mean group

  for (p in probes) {
    m_group <- BiocGenerics::grep(p, BiocGenerics::colnames(m_top_genes_matrix))
    # find columns containing "string"
    group_mean <- rowMeans(subset(m_top_genes_matrix, select = m_group),
                           na.rm = TRUE)
    # get means of each sample and safe as a numeric values
    m_top_genes_matrix <- BiocGenerics::cbind(m_top_genes_matrix, group_mean)
  }
  m_top_genes_matrix <- m_top_genes_matrix[, -seq_len(number_columns)]
  if (ncol(m_top_genes_matrix) == length(probes)) {
    BiocGenerics::colnames(m_top_genes_matrix) <- probes
    return(m_top_genes_matrix)
  } else {
    stop("Number of probes does not match number of columns.")
  }
}

#' Generates K-means for the columns and rows of a matrix.
#'
#' @param m_top_genes_matrix An input matrix for which they K-Means
#' are calculated next.
#' @param k Number k of how many clusters are created.
#'
#' @return All calculated K-Means of the matrix.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' m_top_genes_matrix <- scale_counts(matrixCounts)
#' seed <- 1
#' k <- 1
#' set.seed(seed)
#' k_means <- Kmean_generation(m_top_genes_matrix, k)
Kmean_generation <- function(m_top_genes_matrix, k) {
  # Creating a copy of the matrix with the summarized biological replicates
  y <- m_top_genes_matrix
  # Generating the k-means
  km <- stats::kmeans(y, k)
  m_kmeans <- BiocGenerics::cbind(y, km$cluster)
  last_column <- ncol(m_kmeans)
  # Ordering the genes by their k-means
  o <- BiocGenerics::order(m_kmeans[, last_column])
  m_kmeans <- m_kmeans[o, ]
  return(m_kmeans)
}

#' Function to determine the most variable genes of each cluster to
#' enable annotation..
#'
#' @param m_kmeans Matrix of the k-Means for each row in a gene set.
#' @param number_of_annotations_per_cluster Number of wanted annotations
#' per cluster.
#' @param k Number of clusters.
#'
#' @return A list of the most variable rows (genes) of each cluster.
#' @export
#'
#' @examples
#' x <- 1
#' number_of_annotations_per_cluster <- 5
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' m_top_genes_matrix <- scale_counts(matrixCounts)
#' seed <- 1
#' k <- 1
#' set.seed(seed)
#' m_kmeans <- Kmean_generation(m_top_genes_matrix, k)
#' top_x_variable_genes <- most_variable_genes(m_kmeans,
#' number_of_annotations_per_cluster, k)
most_variable_genes <- function(m_kmeans,
                                number_of_annotations_per_cluster, k) {
  last_column <- ncol(m_kmeans)
  top_x_variable_genes <- list()
  for (i in seq_len(k)) {
    cluster_matrix <- m_kmeans[(m_kmeans[, last_column]) == i, ]
    cluster_matrix <- cluster_matrix[, -last_column]
    # estimates the variance for each row
    variance_row <- apply(cluster_matrix, 1, stats::var)

    # binds a column with the estimated variance of each gene to clustermatrix
    cluster_matrix <- BiocGenerics::cbind(cluster_matrix, variance_row)


    # get number of last column
    number_of_last_column <- ncol(cluster_matrix)

    # orders the variances from highest to lowest
    o <- BiocGenerics::order(cluster_matrix[, number_of_last_column],
                             decreasing = TRUE)
    # orders matrix according to order of variances ( highest to lowest)
    cluster_matrix <- cluster_matrix[o, ]
    # makes a list of the gene names with the highest variance
    cluster_genes_highest_var <- as.list(BiocGenerics::rownames(cluster_matrix)
                                         [seq_len(
                                           number_of_annotations_per_cluster)])

    # creates list with most_variable_genes of all cluster
    top_x_variable_genes <- append(top_x_variable_genes,
                                   cluster_genes_highest_var)
  }
  return(top_x_variable_genes)
}

#' Function to set row annotation for a heatmap.
#'
#' @param m_top_genes_matrix An input matrix to used for heatmap generation and
#' choosing annotation for a heatmap.
#' @param top_x_genes_cluster A list of the most variable x genes of each
#' cluster in a heatmap.
#' @param fontsize_rowAnnotation An integer defining the font size of the
#' row annotation
#'
#' @return A numeric index from the orginal matrix.
#' @export
#'
#' @examples
#' x <- 1
#' number_of_annotations_per_cluster <- 5
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' m_top_genes_matrix <- scale_counts(matrixCounts)
#' seed <- 1
#' k <- 1
#' fontsize_rowAnnotation <- 8
#' set.seed(seed)
#' m_kmeans <- Kmean_generation(m_top_genes_matrix, k)
#' top_x_genes_cluster <- most_variable_genes(m_kmeans,
#' number_of_annotations_per_cluster, k)
#' anno <- set_annotation(m_top_genes_matrix, top_x_genes_cluster,
#' fontsize_rowAnnotation)
#'
set_annotation <- function(m_top_genes_matrix, top_x_genes_cluster,
                           fontsize_rowAnnotation) {
  # get numeric indices of top_x_genes_clusters
  top_x_genes_clusters <- which(BiocGenerics::rownames(m_top_genes_matrix)
                                %in% top_x_genes_cluster)
  labeling <- BiocGenerics::rownames(m_top_genes_matrix)[top_x_genes_clusters]
  anno <- ComplexHeatmap::anno_mark(at = top_x_genes_clusters,
                                    labels = labeling,
                                    labels_gp = grid::gpar(
                                      fontsize = fontsize_rowAnnotation),
                                    which = "row")

  return(anno)
}

#' Function to perform k-Means clustering for a matrix and setting split to
#' split a heatmap into clusters.
#'
#' @param m_top_genes_matrix A matrix for which K-Means are calculated.
#' @param k Number of clusters for k-Means clustering (chosen by elbow-plot).
#'
#' @return A new dataframe containing information assigning each rows of the
#' input matrix to a cluster.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' m_top_genes_matrix <- scale_counts(matrixCounts)
#' k <- 1
#' split <- performing_kMeans(m_top_genes_matrix, k)
performing_kMeans <- function(m_top_genes_matrix, k) {
  # performance of kmean clustering
  kclus <- stats::kmeans(m_top_genes_matrix, k)
  kclus$cluster # check on kmean clusters
  # set split as kmean clusters for heatmap
  split <- paste0("Cluster\n", kclus$cluster)
  return(split)
}

#' Function to build a heatmap using other functions.
#'
#' @param m_top_genes_matrix An input matrix for which the heatmap is created.
#' @param title A string used to set the title of the heatmap.
#' @param split A dataframe containing information about the clustering of
#' the rows of the input matrix.
#' @param anno A numeric index containing row informations for annotation.
#' @param fontsize_columnNames An integer used to set the font size of the
#' columns in the heatmap.
#' @param fontsize_rowNames An integer used to set the font size of the
#' rownames in the heatmap.
#' @param title_heatmapLegend A string setting the title of the
#' legend of the heatmap.
#' @param WidthNum A float setting the width of the heatmap.
#' @param HeightNum A float setting the height of the heatmap.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of
#' HeightNum and WidthNum.
#' @param color_Palette Name of the colorPalette used for the heatmap.
#'
#' @return A plotted heatmap.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' m_top_genes_matrix <- scale_counts(matrixCounts)
#' k <- 1
#' split <- performing_kMeans(m_top_genes_matrix, k)
#' seed <- 1
#' fontsize_rowAnnotation <- 8
#' set.seed(seed)
#' m_kmeans <- Kmean_generation(m_top_genes_matrix, k)
#' number_of_annotations_per_cluster <- 5
#' top_x_genes_cluster <- most_variable_genes(m_kmeans,
#' number_of_annotations_per_cluster, k)
#' anno <- set_annotation(m_top_genes_matrix, top_x_genes_cluster,
#' fontsize_rowAnnotation)
#' title <- "Heatmap of Data"
#' fontsize_columnNames <- 6
#' fontsize_rowNames <- 4
#' title_heatmapLegend <- "Expression"
#' WidthNum <- 4.5
#' HeightNum <- 3
#' UnitSize <- "cm"
#' color_Palette <- "RdBu"
#' set.seed(seed)
#' print_heatmap(m_top_genes_matrix, title, split, anno, fontsize_columnNames,
#' fontsize_rowNames, title_heatmapLegend, WidthNum, HeightNum,
#' UnitSize, color_Palette)
print_heatmap <- function(m_top_genes_matrix, title, split, anno,
                          fontsize_columnNames, fontsize_rowNames,
                          title_heatmapLegend, WidthNum, HeightNum,
                          UnitSize, color_Palette) {
  color_setting(color_Palette)
  ht <- ComplexHeatmap::Heatmap(m_top_genes_matrix,
    name = "mat", split = split,
    column_title = title,
    use_raster = FALSE, cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = fontsize_columnNames),
    row_names_gp = grid::gpar(fontsize = fontsize_rowNames),
    heatmap_legend_param = list(title = title_heatmapLegend),
    width = grid::unit(WidthNum, UnitSize),
    height = grid::unit(HeightNum, UnitSize)
  ) + ComplexHeatmap::rowAnnotation(mark = anno)
  hm <- ComplexHeatmap::draw(ht)
  return(hm)
}

#' Creating a heatmap with annotation of x most variable rows(genes).
#'
#' @param topGenes_matrix An input matrix to create the heatmap.
#' @param probes A list of strings to summarize biological
#' replicated in the columns.
#' @param number_of_annotations_per_cluster An integer used to set the number
#' of annotations per cluster in the heatmap.
#' @param k An integer used to set number of cluster in the heatmap.
#' @param Title A string to set the title of the heatmap.
#' @param fontsize_rowAnnotation An integer used to set the font size of the
#' row annotation in the heatmap
#' @param fontsize_columnNames An integer used to set the font size of the
#' columns in the heatmap.
#' @param fontsize_rowNames An integer used to set the font size of the
#' rownames in the heatmap.
#' @param title_heatmapLegend A string setting the title of the legend
#' of the heatmap.
#' @param WidthNum A float setting the width of the heatmap.
#' @param HeightNum A float setting the height of the heatmap.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of
#' HeightNum and WidthNum.
#' @param color_Palette Name of the colorPalette used for the heatmap.
#'
#' @return A plotted heatmap of the input data.
#' @export
#'
#' @examples
#' x <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' topGenes_matrix <- scale_counts(matrixCounts)
#' probes <- list("disease3", "disease4", "disease1B")
#' k <- 1
#' seed <- 1
#' fontsize_rowAnnotation <- 8
#' number_of_annotations_per_cluster <- 5
#' Title <- "Heatmap of Data"
#' fontsize_columnNames <- 6
#' fontsize_rowNames <- 4
#' title_heatmapLegend <- "Expression"
#' WidthNum <- 4.5
#' HeightNum <- 3
#' UnitSize <- "cm"
#' color_Palette <- "RdBu"
#' set.seed(seed)
#' function_complexHeatmap_var(topGenes_matrix, probes,
#' number_of_annotations_per_cluster, k, Title, fontsize_rowAnnotation,
#' fontsize_columnNames, fontsize_rowNames, title_heatmapLegend, WidthNum,
#' HeightNum, UnitSize, color_Palette)
function_complexHeatmap_var <- function(topGenes_matrix, probes,
                                        number_of_annotations_per_cluster, k,
                                        Title, fontsize_rowAnnotation,
                                        fontsize_columnNames, fontsize_rowNames,
                                        title_heatmapLegend, WidthNum,
                                        HeightNum, UnitSize, color_Palette) {
  m_top_genes_matrix <- summarise_bio_replicates(topGenes_matrix, probes)
  m_kmeans <- Kmean_generation(m_top_genes_matrix, k)
  top_x_genes_cluster <- most_variable_genes(m_kmeans,
                                             number_of_annotations_per_cluster,
                                             k)
  anno <- DgeaHeatmap::set_annotation(m_top_genes_matrix, top_x_genes_cluster,
                                      fontsize_rowAnnotation)
  split <- performing_kMeans(m_top_genes_matrix, k)
  hm <- print_heatmap(m_top_genes_matrix, Title, split, anno,
                      fontsize_columnNames, fontsize_rowNames,
                      title_heatmapLegend, WidthNum, HeightNum,
                      UnitSize, color_Palette)
  return(hm)
}

#' Creating a color scheme based on the available color palettes of
#' RColorBrewer for the heatmap.
#'
#' @param colorPalette An input string defining the color palette
#' (available from RColorBrewer).
#'
#' @return The colors used in heatmap based on either ComplexHeatmap
#' default or RColorBrewer.
#' @export
#'
#' @examples
#' colorPalette <- "RdBu"
#' get_heatmap_colors(colorPalette)
get_heatmap_colors <- function(colorPalette) {
  if (is.null(colorPalette)) {
    # Use ComplexHeatmap default color function
    return(NULL)
  } else {
    # Generate color vector from RColorBrewer palette
    heatmap_color_scheme <- color_setting(colorPalette)
    return(heatmap_color_scheme)
  }
}

#' Creating a color scheme based on the available color palettes of
#' RColorBrewer for the heatmap.
#'
#' @param ncounts_matrix An input matrix to create the heatmap.
#' @param column_name A string to set the title of the heatmap,
#' default <- "Heatmap".
#' @param colorPalette Name of the colorPalette used for the heatmap,
#' default <- NULL.
#' @param cluster_method A string setting the cluster method for the heatmap,
#' default <- "hierarchical".
#' @param distance_method A string setting the distance method for clustering,
#' default <- "euclidean".
#' @param cluster_rows A Boolean switching the optional clustering of the rows
#' on and off, default <- TRUE.
#' @param cluster_columns A Boolean switching the optional clustering of the
#' columns on and off, default <- FALSE.
#' @param k_row An integer used to set number of clusters for row clustering
#' in the heatmap, default <- NULL.
#' @param k_col An integer used to set number of clusters for column clustering
#' in the heatmap, default <- NULL.
#' @param sample_metadata A dataframe containing the metadata information for
#' the grouping of the columns, default <- NULL.
#' @param annotation_colors A list assigning choosen colors to the corresponding
#' groups, default <- NULL.
#' @param annotation_name_side The side of the column annotation description,
#' default <- "right".
#' @param show_row_names A Boolean switching rownames on the heatmap on and off,
#' default <- FALSE.
#' @param show_column_names A Boolean switching colum names on the heatmap on
#' and off, default <- TRUE.
#' @param row_annotation A Boolean switching row annotation on the heatmap on
#' and off, default <- FALSE.
#' @param row_annotation_method A string setting the annotation method of the
#' heatmap, default <- "auto".
#' @param row_anno_names A list containing choosen rownames to use for the row
#' annotation, default <- NULL.
#' @param row_anno_number An integer setting the number of automatic annotations
#' assigned per cluster, default <- 5.
#' @param fontsize_title An integer setting the font size of the heatmap title.
#' @param fontsize_rowAnnotation An integer setting the font size of the
#' optional row annotation, default <- 10.
#' @param fontsize_columnNames An integer setting the font size of the
#' column names, default <- 6.
#' @param fontsize_rowNames An integer setting the font size of the
#' row names, default <- 4.
#' @param fontsize_cluster_labels An integer setting the font size of the
#' cluster labels, default <- 8.
#' @param fontsize_group_annotation An integer setting the font size of the
#' group annotation title, default <- 8.
#' @param fontsize_group_annotation_legend An integer setting the font size of
#' the group annotation legend name, default <- 10.
#' @param fontsize_group_annotation_labels An integer setting the font size of
#' the group annotation labels in the legend, default <- 8.
#' @param fontsize_heatmap_legend An integer setting the font size of the
#' heatmap legend title, default <- 10.
#' @param fontsize_heatmap_legend_labels An integer setting the font size of
#' the heatmap legend labels, default <- 8.
#' @param title_heatmapLegend A string setting the changeable title of the
#' legend, default "Expression".
#' @param WidthNum A float setting the width of the heatmap, default <- 4.5.
#' @param HeightNum A float setting the height of the heatmap, default <- 3.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of HeightNum
#' and WidthNum, default <- "cm".
#'
#' @return An advanced and customizable heatmap.
#' @export
#'
#' @examples
#' colorPalette <- "RdBu"
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed
#' adv_Heatmap(ncounts_matrix)
adv_Heatmap <- function(ncounts_matrix, column_name = "Heatmap",
                        colorPalette = NULL, cluster_method = "hierarchical",
                        distance_method = "euclidean", cluster_rows = TRUE,
                        cluster_columns = FALSE, k_row = NULL,
                        k_col = NULL, sample_metadata = NULL,
                        annotation_colors = NULL,
                        annotation_name_side = "right", show_row_names = FALSE,
                        show_column_names = TRUE, row_annotation = FALSE,
                        row_annotation_method = "auto", row_anno_names = NULL,
                        row_anno_number = 5, fontsize_title = 15,
                        fontsize_rowAnnotation = 10, fontsize_columnNames = 6,
                        fontsize_rowNames = 4, fontsize_cluster_labels = 8,
                        fontsize_group_annotation = 8,
                        fontsize_group_annotation_legend = 10,
                        fontsize_group_annotation_labels = 8,
                        fontsize_heatmap_legend = 10,
                        fontsize_heatmap_legend_labels = 8,
                        title_heatmapLegend = "Expression", WidthNum = 4.5,
                        HeightNum = 3, UnitSize = "cm") {
  # Clustering logic
  row_split <- col_split <- NULL
  row_dend <- col_dend <- TRUE # default TRUE if unspecified

  # --- ROW CLUSTERING ---
  rows_clustered <- row_clustering(ncounts_matrix,
                                   cluster_rows = cluster_rows,
                                   cluster_method = cluster_method,
                                   distance_method = distance_method,
                                   k_row = k_row, row_split = row_split,
                                   row_dend = row_dend)
  row_split <- rows_clustered$row_split
  row_dend <- rows_clustered$dend
  # --- COLUMN CLUSTERING ---
  columns_clustered <- column_clustering(ncounts_matrix,
                                         cluster_columns = cluster_columns,
                                         cluster_method = cluster_method,
                                         distance_method = distance_method,
                                         k_col = k_col, col_split = col_split,
                                         col_dend = col_dend)
  col_split <- columns_clustered$col_split
  col_dend <- columns_clustered$dend
  # --- Sample Annotation ---
  col_ha <- NULL
  col_ha <- set_sample_annotation(
    sample_metadata = sample_metadata,
    annotation_colors = group_colors,
    annotation_name_side = annotation_name_side,
    fontsize_group_annotation = fontsize_group_annotation,
    fontsize_group_annotation_legend = fontsize_group_annotation_legend,
    fontsize_group_annotation_labels = fontsize_group_annotation_labels)

  # --- Row Annotation ---
  annotation_for_rows <- NULL
  annotation_for_rows <- set_row_annotation(
    ncounts_matrix = ncounts_matrix,
    k_row = k_row, row_annotation = row_annotation,
    row_annotation_method = row_annotation_method,
    row_anno_names = row_anno_names,
    row_anno_number = row_anno_number,
    fontsize_rowAnnotation = fontsize_rowAnnotation)

  # --- Final Ordering Check --- #
  # Reset column split if clustering is not active
  if (!cluster_columns) {
    col_split <- NULL
    col_dend <- FALSE
  }
  if (!cluster_rows) {
    row_split <- NULL
    row_dend <- FALSE
  }
  # --- Draw Heatmap ---
  ht <- draw_adv_heatmap(
    ncounts_matrix, column_name = column_name,
    colorPalette = colorPalette, show_row_names = show_row_names,
    show_column_names = show_column_names, fontsize_title = fontsize_title,
    fontsize_columnNames = fontsize_columnNames,
    fontsize_rowNames = fontsize_rowNames,
    fontsize_cluster_labels = fontsize_cluster_labels,
    fontsize_heatmap_legend = fontsize_heatmap_legend,
    fontsize_heatmap_legend_labels = fontsize_heatmap_legend_labels,
    title_heatmapLegend = title_heatmapLegend, WidthNum = WidthNum,
    HeightNum = HeightNum, UnitSize = UnitSize, row_annotation = row_annotation,
    annotation_for_rows = annotation_for_rows, row_split = row_split,
    col_split = col_split, row_dend = row_dend, col_dend = col_dend,
    col_ha = col_ha)
  hm <- ComplexHeatmap::draw(ht)
  return(hm)
}


#' Setting the row clustering
#'
#' @param ncounts_matrix An input matrix to create the clustering.
#' @param cluster_rows A Boolean switching the optional clustering of the rows
#' on and off, default <- TRUE.
#' @param cluster_method A string setting the cluster method for the heatmap.
#' @param distance_method A string setting the distance method for clustering.
#' @param k_row An integer used to set number of clusters for row clustering.
#' @param row_split An integer seeting the column clustering by kmeans.
#' @param row_dend A Boolean switching column clustering by hierarchical
#' clustering on (default).
#'
#' @return row_split where the rows are split into clusters
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed
#' row_split <- NULL
#' row_dend <- TRUE
#' rows_clustered <- row_clustering(ncounts_matrix, cluster_rows = TRUE,
#' cluster_method = "hierarchical", distance_method = "euclidean", k_row = NULL,
#' row_split = row_split, row_dend = row_dend)
row_clustering <- function(ncounts_matrix, cluster_rows = TRUE,
                           cluster_method = "hierarchical",
                           distance_method = "euclidean",
                           k_row = NULL, row_split = NULL, row_dend = TRUE) {
  row_split <- row_split
  row_dend <- row_dend # default TRUE if unspecified
  if (cluster_rows) {
    row_data <- ncounts_matrix

    if (cluster_method == "hierarchical") {
      row_dist <- get_dist(row_data, distance_method)
      row_clust <- stats::hclust(row_dist, method = "complete")
      row_dend <- stats::as.dendrogram(row_clust)
      if (!is.null(k_row)) {
        row_split <- stats::cutree(row_clust, k = k_row)
        row_split <- as.factor(row_split)
        row_dend <- TRUE
      }
    } else if (cluster_method == "kmeans") {
      row_split <- performing_kMeans(row_data, k_row)
    }
  }
  return(list(row_split = row_split, dend = row_dend))
}

#' Calculates the distance matrix
#'
#' @param x An input matrix to create the distance matrix.
#' @param method Astring setting the distance method for clustering.

#'
#' @return A distance matrix.
#' @export
#'
#' @examples
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' x <- 1
#' matrixCounts <- build_matrix(input_data, x)
#' matrix <- scale_counts(matrixCounts)
#' method <- "euclidean"
#' distance_matrix <- get_dist(matrix, method)
get_dist <- function(x, method) {
  if (method == "correlation") {
    stats::as.dist(1 - stats::cor(t(x)))
  } else {
    stats::dist(x, method = method)
  }
}

#' Setting the column clustering
#'
#' @param ncounts_matrix An input matrix to create the clustering.
#' @param cluster_columns A Boolean switching the optional clustering of the
#' rows on and off, default <- TRUE.
#' @param cluster_method A string setting the cluster method for the heatmap.
#' @param distance_method A string setting the distance method for clustering.
#' @param k_col An integer used to set number of clusters for row clustering.
#' @param col_split An integer seeting the column clustering by kmeans.
#' @param col_dend A Boolean switching column clustering by hierarchical
#' clustering on (default).
#'
#' @return col_split where the rows are split into clusters
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed
#' col_split <- NULL
#' col_dend <- TRUE
#' columns_clustered <- column_clustering(ncounts_matrix,
#' cluster_columns = TRUE, cluster_method = "hierarchical",
#' distance_method = "euclidean", k_col = NULL, col_split = col_split,
#' col_dend = col_dend)
column_clustering <- function(ncounts_matrix, cluster_columns = TRUE,
                              cluster_method = "hierarchical",
                              distance_method = "euclidean", k_col = NULL,
                              col_split = NULL, col_dend = TRUE) {
  col_split <- col_split
  col_dend <- col_dend
  if (cluster_columns) {
    col_data <- t(ncounts_matrix)

    if (cluster_method == "hierarchical") {
      col_dist <- get_dist(col_data, distance_method)
      col_clust <- stats::hclust(col_dist, method = "average")
      col_dend <- stats::as.dendrogram(col_clust)
      if (!is.null(k_col)) {
        col_split <- stats::cutree(col_clust, k <- k_col)
        col_split <- as.factor(col_split)
        col_dend <- TRUE
      }
    } else if (cluster_method == "kmeans") {
      col_split <- performing_kMeans(col_data, k_col)
    }
  }
  return(list(col_split = col_split, dend = col_dend))
}

#' Setting the sample annotation
#'
#' @param sample_metadata A dataframe containing the metadata information for
#' the grouping of the columns, default <- NULL.
#' @param annotation_colors A list assigning choosen colors to the corresponding
#' groups, default <- NULL.
#' @param annotation_name_side The side of the column annotation description,
#' default <- "right".
#' @param fontsize_group_annotation An integer setting the font size of the
#' group annotation title, default <- 8.
#' @param fontsize_group_annotation_legend An integer setting the font size of
#' the group annotation legend name, default <- 10.
#' @param fontsize_group_annotation_labels An integer setting the font size of
#' the group annotation labels in the legend, default <- 8.
#'
#' @return A "HeamapAnnotation" object if "sample_metadata" is provided,
#' otherwise NULL.
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed
#' groups <- c("3_DKD_glomerulus_Geometric_S", "1B_DKD_glomerulus_Geometric_S",
#' "2B_DKD_glomerulus_WT")
#' sample_names <- c(colnames(ncounts_matrix))
#' group_assignment <- sapply(sample_names, function(sample) {
#'   matched <- groups[sapply(groups, function(g) grepl(g, sample))]
#'   if (length(matched) > 0) matched[1] else NA
#' })
#' stopifnot(length(sample_names) == length(group_assignment))
#' sample_metadata <- data.frame(Group = group_assignment,
#' row.names = sample_names)
#' all(colnames(ncounts_matrix) == rownames(sample_metadata))
#' group_colors <- list(Group = c("3_DKD_glomerulus_Geometric_S" = "#1b9e77",
#' "1B_DKD_glomerulus_Geometric_S" = "#7570b3",
#' "2B_DKD_glomerulus_WT" = "#e7298a"))
#' col_ha <- set_sample_annotation(sample_metadata = sample_metadata,
#' annotation_colors = group_colors)
set_sample_annotation <- function(sample_metadata = NULL,
                                  annotation_colors = NULL,
                                  annotation_name_side = "right",
                                  fontsize_group_annotation = 8,
                                  fontsize_group_annotation_legend = 10,
                                  fontsize_group_annotation_labels = 8
) {
  col_ha <- NULL
  if (!is.null(sample_metadata)) {
    col_ha <- ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(sample_metadata),
      col = annotation_colors,
      annotation_name_side = annotation_name_side,
      annotation_name_gp = grid::gpar(fontsize = fontsize_group_annotation),
      annotation_legend_param = list(
        title_gp = grid::gpar(fontsize = fontsize_group_annotation_legend),
        labels_gp = grid::gpar(fontsize = fontsize_group_annotation_labels)
      )
    )
  }
  return(col_ha)
}

#' Setting the row annotation
#'
#' @param ncounts_matrix An input matrix to create the heatmap.
#' @param row_annotation A Boolean switching row annotation on the heatmap on
#' and off, default <- FALSE.
#' @param row_annotation_method A string setting the annotation method of the
#' heatmap, default <- "auto".
#' @param row_anno_names A list containing choosen rownames to use for the row
#' annotation, default <- NULL.
#' @param row_anno_number An integer setting the number of automatic annotations
#' assigned per cluster, default <- 5.
#' @param k_row An integer used to set number of clusters for row clustering in
#' the heatmap, default <- NULL.
#' @param fontsize_rowAnnotation An integer setting the font size of the group
#' annotation labels in the legend, default <- 8.
#'
#' @return Numeric index from matrix for the row annotation
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv",
#' package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed(seed)
#' annotation_for_rows <- set_row_annotation(ncounts_matrix)
set_row_annotation <- function(ncounts_matrix,
                               k_row = NULL,
                               row_annotation = FALSE, # FALSE or TRUE
                               row_annotation_method = "auto",
                               row_anno_names = NULL,
                               row_anno_number = 5,
                               fontsize_rowAnnotation = 10) {
  annotation_for_rows <- NULL
  if (row_annotation) {
    if (row_annotation_method == "auto") {
      if (!is.null(k_row)) {
        m_kmeans <- Kmean_generation(ncounts_matrix, k_row)
        top_x_genes_cluster <- most_variable_genes(m_kmeans,
                                                   row_anno_number, k_row)
        annotation_for_rows <- set_annotation(ncounts_matrix,
                                              top_x_genes_cluster,
                                              fontsize_rowAnnotation)
      } else {
        y <- ncounts_matrix
        top_x_genes_cluster <- list()
        # estimates the variance for each row
        variance_row <- apply(y, 1, stats::var)
        # binds column with the estimated variance of each gene to clustermatrix
        cluster_matrix <- BiocGenerics::cbind(y, variance_row)
        # get number of last column
        number_of_last_column <- ncol(y)
        # orders the variances from highest to lowest
        o <- BiocGenerics::order(y[, number_of_last_column], decreasing = TRUE)
        # orders matrix according to order of variances ( highest to lowest)
        cluster_matrix <- y[o, ]
        # makes a list of the gene names with the highest variance
        cluster_genes_highest_var <- as.list(
          BiocGenerics::rownames(y)[seq_len(row_anno_number)])
        # creates list with most_variable_genes of all cluster
        top_x_genes_cluster <- append(top_x_genes_cluster,
                                      cluster_genes_highest_var)
        # set annotation for rows
        annotation_for_rows <- set_annotation(ncounts_matrix,
                                              top_x_genes_cluster,
                                              fontsize_rowAnnotation)
      }
    } else if (row_annotation_method == "specific") {
      annotation_for_rows <- set_annotation(ncounts_matrix, row_anno_names,
                                            fontsize_rowAnnotation)
    }
  }
  return(annotation_for_rows)
}

#' Draw the advanced heatmap
#'
#' @param ncounts_matrix An input matrix to create the heatmap.
#' @param column_name A string to set the title of the heatmap.
#' @param colorPalette Name of the colorPalette used for the heatmap, default <- NULL.
#' @param show_row_names A Boolean switching rownames on the heatmap on and off, default <- FALSE.
#' @param show_column_names A Boolean switching colum names on the heatmap on and off, default <- TRUE.
#' @param fontsize_title An integer setting the font size of the heatmap title, default <- 15.
#' @param fontsize_columnNames An integer setting the font size of the column names, default <- 6.
#' @param fontsize_rowNames An integer setting the font size of the row names, default <- 4.
#' @param fontsize_cluster_labels An integer setting the font size of the cluster labels, default <- 8.
#' @param fontsize_heatmap_legend An integer setting the font size of the heatmap legend title, default <- 10.
#' @param fontsize_heatmap_legend_labels An integer setting the font size of the heatmap legend labels, default <- 8.
#' @param title_heatmapLegend A string setting the changeable title of the legend, default "Expression".
#' @param WidthNum A float setting the width of the heatmap, default <- 4.5.
#' @param HeightNum A float setting the height of the heatmap, default <- 3.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of HeightNum and WidthNum, default <- "cm".
#' @param annotation_for_rows Numeric index from matrix for the row annotation
#' @param row_annotation A Boolean switching row annotation on the heatmap on and off, default <- FALSE.
#' @param row_split An integer seeting the row clustering by kmeans.
#' @param col_split An integer seeting the column clustering by kmeans.
#' @param col_ha A "HeamapAnnotation" object if "sample_metadata" is provided, otherwise NULL.
#' @param row_dend A Boolean switching row clustering by hierarchical clustering on (default).
#' @param col_dend A Boolean switching column clustering by hierarchical clustering on (default).
#'
#'
#' @return Row annotation
#' @export
#'
#' @examples
#' x <- 1
#' seed <- 1
#' input_data <- read.csv(system.file("extdata/testfile_counts.csv", package = "DgeaHeatmap"))
#' matrixCounts <- build_matrix(input_data, x)
#' ncounts_matrix <- scale_counts(matrixCounts)
#' set.seed(seed)
#' row_split <- NULL
#' row_dend <- TRUE
#' col_split <- NULL
#' col_dend <- TRUE
#' groups <- c("3_DKD_glomerulus_Geometric_S", "1B_DKD_glomerulus_Geometric_S", "2B_DKD_glomerulus_WT")
#' sample_names <- c(colnames(ncounts_matrix))
#' group_assignment <- sapply(sample_names, function(sample) {
#'   matched <- groups[sapply(groups, function(g) grepl(g, sample))]
#'   if (length(matched) > 0) matched[1] else NA
#' })
#' stopifnot(length(sample_names) == length(group_assignment))
#' sample_metadata <- data.frame(Group = group_assignment, row.names = sample_names)
#' all(colnames(ncounts_matrix) == rownames(sample_metadata))
#' group_colors <- list(Group = c("3_DKD_glomerulus_Geometric_S" = "#1b9e77", "1B_DKD_glomerulus_Geometric_S" = "#7570b3", "2B_DKD_glomerulus_WT" = "#e7298a"))
#' col_ha <- set_sample_annotation(sample_metadata = sample_metadata, annotation_colors = group_colors)
#' ht <- draw_adv_heatmap(ncounts_matrix)
draw_adv_heatmap <- function(ncounts_matrix, column_name = "Heatmap", colorPalette = NULL, col_ha = NULL, show_row_names = FALSE, show_column_names = TRUE, fontsize_title = 15, fontsize_columnNames = 6, fontsize_rowNames = 4, fontsize_cluster_labels = 8, fontsize_heatmap_legend = 10, fontsize_heatmap_legend_labels = 8, title_heatmapLegend = "Expression", WidthNum = 4.5, HeightNum = 3, UnitSize = "cm", row_annotation = FALSE, annotation_for_rows = NULL, row_split = NULL, col_split = NULL, row_dend = TRUE, col_dend = TRUE) {
  heatmap_color_scheme <- get_heatmap_colors(colorPalette) # sets the color Palette for the heatmap
  ht <- ComplexHeatmap::Heatmap(ncounts_matrix,
    cluster_rows = if (is.logical(row_dend)) row_dend else row_dend,
    cluster_columns = if (is.logical(col_dend)) col_dend else col_dend,
    col = heatmap_color_scheme,
    row_split = row_split,
    column_split = col_split,
    top_annotation = col_ha,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    column_title = column_name,
    column_title_gp = grid::gpar(fontsize = fontsize_title),
    use_raster = FALSE,
    column_names_gp = grid::gpar(fontsize = fontsize_columnNames),
    row_names_gp = grid::gpar(fontsize = fontsize_rowNames),
    row_title_gp = grid::gpar(fontsize = fontsize_cluster_labels),
    heatmap_legend_param = list(title = title_heatmapLegend, title_gp = grid::gpar(fontsize = fontsize_heatmap_legend), labels_gp = grid::gpar(fontsize = fontsize_heatmap_legend_labels)),
    width = grid::unit(WidthNum, UnitSize),
    height = grid::unit(HeightNum, UnitSize)
  )

    if (isTRUE(row_annotation)) {
    ht <- ht + ComplexHeatmap::rowAnnotation(mark = annotation_for_rows)
  }
  hm <- ComplexHeatmap::draw(ht)
  return(hm)
}
