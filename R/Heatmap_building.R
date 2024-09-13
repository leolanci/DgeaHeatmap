#' Build a matrix from an input .csv file, with one column as rownames
#'
#' @param counts_data A data frame with floats or integers that has strings in one column.
#' @param x Number of column which contains the rownames description.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' counts_data <- read.csv("Counts_data.csv")
#' x <- 1

build_matrix <- function(counts_data, x){
  BiocGenerics::colnames(counts_data)[x] <- "genes"
  list_rownames <- as.list(counts_data[,x])
  print(list_rownames)
  counts_data <- as.matrix(counts_data[,-x])
  rownames(counts_data) <-list_rownames$genes
  return(counts_data)
}

#' Creates a matrix only containing chosen columns of an original matrix with more data.
#'
#' @param factors_for_matrix_devision A list containing variables used to extract specific columns from the original matrix.
#' @param mmatrix A matrix used as input, from which specifically chosen columns are extracted.
#'
#' @return A new matrix containing fewer columns than before, that have previously been specified.
#' @export
#'
#' @examples
#' factors_for_matrix <- list("Female", "6_month")
#' mmatrix <- matrix
individual_matrix <- function(factors_for_matrix_devision, mmatrix) {
  for (i in factors_for_matrix_devision) {
    chosen_columns <- BiocGenerics::grep(i, BiocGenerics::colnames(mmatrix))
    rmatrix <- BiocGenerics::cbind(mmatrix[, chosen_columns])
    mmatrix <- rmatrix
    rmatrix <- matrix()
  }
  return(mmatrix)
}

#' Function to filter a matrix to extract a chosen number of most variable rows through calculation of the variance.
#'
#' @param counts_data An input matrix from which the most variable rows are to be extracted.
#' @param top_number_of_genes An integer to set the number of how many rows should be extracted from the original matrix by variance.
#'
#' @return A new matrix containing only the top_number_of_genes count of rows with the highest variance from input matrix.
#' @export
#'
#' @examples
#' top_number_of_genes <- 500
#' counts_data <- matrix
filtering_for_top_exprGenes <- function(counts_data, top_number_of_genes){
  var_genes <- apply(counts_data, 1, var) #estimating the variance of each gene
  # Sorting the genes by their variance and creating a new object with chosen number of most variable genes
  select_var <- names(sort(var_genes, decreasing = TRUE))[1:top_number_of_genes]
  anaMatrix <- counts_data
  highly_variable_genes <- anaMatrix[select_var,]
  dim(highly_variable_genes)

  return(highly_variable_genes)
}

#' Function to Z-count scale the values of a matrix.
#'
#' @param countsmatrix An input matrix whose values are scaled by Z-count.
#'
#' @return A new matrix with Z-count scaled values.
#' @export
#'
#' @examples
#' countsmatrix <- matrix

scale_counts <- function(countsmatrix){
  scaled_counts =
    countsmatrix%>%
    t(.) %>% #transpose to have genes in columns
    scale() %>%  #scale (x, center = TRUE, scale = TRUE)
    t(.)  #'transpose back in original shape
  return(scaled_counts)
}

#' Function to visualize the data distribution within the data of a matrix.
#'
#' @param scaled_counts An input matrix with Z-count scaled values.
#'
#' @return A plot visualizing the data distribution of Z-count scaled data.
#' @export
#'
#' @examples
#' scaled_counts <- scale_counts(matrix)
show_data_distribution <- function(scaled_counts){
  apply(scaled_counts, MARGIN = 1, mean) %>%  #calculate the mean per row
    graphics::hist(., main = "", xlab = "Z-score values", col = "dodgerblue2") #build histogram to see data distribution
}

#' Function to create an elbow plot to choose k for clustering by k-Means.
#'
#' @param seed An integer to set seed for the k-Mean generation, allows reproducibility.
#' @param top_genes_matrix A matrix for which the best number of clusters in k-Means Clustering is supposed to be calculated.
#'
#' @return An elbow plot to choose k.
#' @export
#'
#' @examples
#' seed <- 1
#' top_genes_matrix <- matrix
elbow_plot <- function(seed, top_genes_matrix){
  set.seed(seed) # reproducible outcome
  wss <- function(k) {
    stats::kmeans(top_genes_matrix, k, nstart = 10)$tot.withinss
  }
  k.values <- 1:15 # compute and plot wws for k = 1 to k = 15
  wss_values <- purrr::map_dbl(k.values, wss)
  plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")

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
color_setting <- function(colorPalette){
  my_colors = RColorBrewer::brewer.pal(n=11, name = colorPalette)
  my_colors = grDevices::colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  my_colors
}

#' Summarizes columns biological replicates of a matrix into one.
#'
#' @param top_genes_matrix An input matrix with columns representing biological replicates.
#' @param probes List of strings by which columns of biological replicates can be identified.
#'
#' @return A new matrix where biological replicated are summarized as their mean.
#' @export
#'
#' @examples
#' top_genes_matrix <- matrix
#' probes <- list("Region1", "Region2", "Region3")


summarise_bio_replicates <- function(top_genes_matrix, probes){
  m_top_genes_matrix <- top_genes_matrix
  number_columns <- ncol(m_top_genes_matrix)
  # loop to grep columns for each mean group

  for (p in probes){
    m_group <- BiocGenerics::grep(p, BiocGenerics::colnames(m_top_genes_matrix))
    # find columns containing "string"
    group_mean <- rowMeans(subset(m_top_genes_matrix, select = m_group), na.rm = TRUE)
    # get means of each sample and safe as a numeric values
    m_top_genes_matrix <- BiocGenerics::cbind(m_top_genes_matrix, group_mean)
  }
  m_top_genes_matrix <- m_top_genes_matrix[,-1:-number_columns]
  if (ncol(m_top_genes_matrix) == length(probes)) {
    BiocGenerics::colnames(m_top_genes_matrix) <- probes
    return(m_top_genes_matrix)
  } else{
    print("Error number of probes does not match number of columns!")
  }
}

#' Generates K-means for the columns and rows of a matrix.
#'
#' @param m_top_genes_matrix An input matrix for which they K-Means are calculated next.
#' @param seed An integer to set the seed for the K-Mean calculation, allows reproducibility.
#' @param k Number k of how many clusters are created.
#'
#' @return All calculated K-Means of the matrix.
#' @export
#'
#' @examples
#' m_top_genes_matrix <- matrix
#' seed <- 1
#' k <- 1

Kmean_generation <- function(m_top_genes_matrix, seed,k){
  # Creating a copy of the matrix with the summarized biological replicates
  y <- m_top_genes_matrix
  # Setting a random factor for reproducibility
  set.seed(seed)
  # Generating the k-means
  km <- stats::kmeans(y, k)
  m_kmeans <- BiocGenerics::cbind(y, km$cluster)
  last_column <- ncol(m_kmeans)
  # Ordering the genes by their k-means
  o <- BiocGenerics::order(m_kmeans[, last_column])
  m_kmeans <- m_kmeans[o,]
  return(m_kmeans)
}

#' Function to determine the most variable genes of each cluster to enable annotation..
#'
#' @param m_kmeans Matrix of the k-Means for each row in a gene set.
#' @param number_of_annotations_per_cluster Number of wanted annotations per cluster.
#' @param k Number of clusters.
#'
#' @return A list of the most variable rows (genes) of each cluster.
#' @export
#'
#' @examples
#' m_kmeans <- matrix_with_kmeans
#' number_of_annotations_per_cluster <- 10
#' k <- 1
most_variable_genes <- function(m_kmeans,number_of_annotations_per_cluster, k){
  last_column <- ncol(m_kmeans)
  top_x_variable_genes <- list()
  for (i in 1:k){
    cluster_matrix <- m_kmeans[(m_kmeans[,last_column]) == i,]
    cluster_matrix <- cluster_matrix[,-last_column]
    #estimates the variance for each row
    variance_row <- apply(cluster_matrix, 1, var)

    # binds a column with the estimated variance of each gene to clustermatrix
    cluster_matrix <- BiocGenerics::cbind(cluster_matrix, variance_row)


    #get number of last column
    number_of_last_column <- ncol(cluster_matrix)

    # orders the variances from highest to lowest
    o <- BiocGenerics::order(cluster_matrix[, number_of_last_column], decreasing = TRUE)
    # orders matrix according to order of variances ( highest to lowest)
    cluster_matrix <- cluster_matrix[o,]
    # makes a list of the gene names with the highest variance
    cluster_genes_highest_var <- as.list(BiocGenerics::rownames(cluster_matrix)[1:number_of_annotations_per_cluster])

    # creates list with most_variable_genes of all cluster
    top_x_variable_genes <- append (top_x_variable_genes, cluster_genes_highest_var)

  }
  return(top_x_variable_genes)
}

#' Function to set annotation for a heatmap.
#'
#' @param m_top_genes_matrix An input matrix to used for heatmap generation and choosing annotation for a heatmap.
#' @param top_x_genes_cluster A list of the most variable x genes of each cluster in a heatmap.
#'
#' @return A numeric index from the orginal matrix.
#' @export
#'
#' @examples
#' m_top_genes_matrix <- matrix
#' top_x_genes_cluster <- list(genes)
set_annotation <- function(m_top_genes_matrix, top_x_genes_cluster){
  #get numeric indices of top_x_genes_clusters
  top_x_genes_clusters <- which(BiocGenerics::rownames(m_top_genes_matrix) %in% top_x_genes_cluster)
  #set row annotation: at = numeric indices of wanted labels, labels = names of wanted labels, which = rows or columns
  labeling = BiocGenerics::rownames(m_top_genes_matrix)[top_x_genes_clusters]
  anno = ComplexHeatmap::anno_mark(at = top_x_genes_clusters, labels = labeling, which = "row")

  return(anno)
}

#' Function to perform k-Means clustering for a matrix and setting split to split a heatmap into clusters.
#'
#' @param m_top_genes_matrix A matrix for which K-Means are calculated.
#' @param k Number of clusters for k-Means clustering (chosen by elbow-plot).
#'
#' @return A new dataframe containing information assigning each rows of the input matrix to a cluster.
#' @export
#'
#' @examples
#' m_top_genes_matrix <- matrix
#' k <-- 3
performing_kMeans <- function(m_top_genes_matrix, k){
  kclus <- stats::kmeans(m_top_genes_matrix, k)  #performance of kmean clustering
  kclus$cluster #check on kmean clusters
  split <- paste0("Cluster\n", kclus$cluster) #set split as kmean clusters for heatmap
  return(split)
}

#' Function to build a heatmap using other functions.
#'
#' @param seed An integer used to set seed for k-Means clustering.
#' @param m_top_genes_matrix An input matrix for which the heatmap is created.
#' @param title A string used to set the title of the heatmap.
#' @param split A dataframe containing information about the clustering of the rows of the input matrix.
#' @param fontsize_columnNames An integer used to set the font size of the columns in the heatmap.
#' @param fontsize_rowNames An integer used to set the font size of the rownames in the heatmap.
#' @param title_heatmapLegend A string setting the title of the legend of the heatmap.
#' @param WidthNum A float setting the width of the heatmap.
#' @param HeightNum A float setting the height of the heatmap.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of HeightNum and WidthNum.
#'
#' @return A plotted heatmap.
#' @export
#'
#' @examples
#' seed <- 1
#' m_top_genes_matrix <- matrix
#' title <- "Heatmap of Data"
#' split <- performing_kMeans(m_top_genes_matrix)
#' fontsize_columnNames <- 6
#' fontsize_rowNames <- 4
#' title_heatmaoLegend <- "Expression"
#' WidthNum <- 4.5
#' HeightNum <- 3
#' UnitSize <- "cm"
print_heatmap <- function(seed,m_top_genes_matrix, title, split, anno,fontsize_columnNames,fontsize_rowNames, title_heatmapLegend,WidthNum, HeightNum, UnitSize){
  color_setting()
  set.seed(seed)
  ht = ComplexHeatmap::Heatmap(m_top_genes_matrix, name = "mat", split = split,
               column_title = title ,
               use_raster = FALSE, cluster_columns = FALSE,
               column_names_gp = grid::gpar(fontsize = fontsize_columnNames),
               row_names_gp = grid::gpar(fontsize = fontsize_rowNames), heatmap_legend_param = list(title = title_heatmapLegend),
               width = grid::unit(WidthNum, UnitSize),
               height = grid::unit(HeightNum, UnitSize)) + ComplexHeatmap::rowAnnotation(mark = anno)
  hm = ComplexHeatmap::draw(ht)
  return(hm)
}


function_complexHeatmap_var <- function(topGenes_matrix, probes, number_of_annotations_per_cluster, k, seed, Title){
  m_top_genes_matrix <- summarise_bio_replicates(topGenes_matrix, probes)

  m_kmeans <- Kmean_generation(m_top_genes_matrix, seed, k)
  top_x_genes_cluster <- most_variable_genes(m_kmeans, number_of_annotations_per_cluster, k)
  anno <- set_annotation(m_top_genes_matrix, top_x_genes_cluster)
  split <- performing_kMeans(m_top_genes_matrix, k)
  hm <- print_heatmap(seed, m_top_genes_matrix, Title, split, anno)
  return(hm)
}

function_complexHeatmap_woSum <- function(topGenes_matrix, probes, number_of_annotations_per_cluster, k, seed, Title){

  m_kmeans <- Kmean_generation(topGenes_matrix, seed, k)
  top_x_genes_cluster <- most_variable_genes(m_kmeans, number_of_annotations_per_cluster, k)
  anno <- set_annotation(topGenes_matrix, top_x_genes_cluster)
  split <- performing_kMeans(topGenes_matrix, k)
  hm <- print_heatmap(seed, topGenes_matrix, Title, split, anno)
  return(hm)
}

