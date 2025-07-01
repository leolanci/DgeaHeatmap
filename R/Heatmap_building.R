#' Builds a matrix from an input csv file, with one column as rownames
#'
#' @param counts_data A data frame with floats or integers that has strings in one column.
#' @param x Number of column which contains the rownames description.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' counts_data <- file
#' x <- 1

build_matrix <- function(counts_data, x){
  #BiocGenerics::colnames(counts_data)[x] <- "genes"
  list_rownames <- as.list(counts_data[,x])
  print(list_rownames)
  counts_data <- as.matrix(counts_data[,-x])
  rownames(counts_data) <-list_rownames
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
  var_genes <- apply(counts_data, 1, stats::var) #estimating the variance of each gene
  # Sorting the genes by their variance and creating a new object with chosen number of most variable genes
  select_var <- names(sort(var_genes, decreasing = TRUE))[1:top_number_of_genes]
  highly_variable_genes <- counts_data[select_var,]
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
#' @importFrom rlang .data
#'
#' @param scaled_counts An input matrix with Z-count scaled values.
#'
#' @return A plot visualizing the data distribution of Z-count scaled data.
#' @export
show_data_distribution <- function(scaled_counts){
  apply(scaled_counts, MARGIN = 1, mean) %>%  #calculate the mean per row
    graphics::hist(., main = "", xlab = "Z-score values", col = "dodgerblue2") #build histogram to see data distribution
}



#' Function to create an elbow plot to choose k for clustering by k-Means.
#'
#' @param seed An integer to set seed for the k-Mean generation, allows reproducibility.
#' @param top_genes_matrix A matrix for which the best number of clusters in k-Means Clustering is supposed to be calculated.
#' @param maxK An integer determining the maximum number of possible clusters.
#'
#' @return An elbow plot to choose k.
#' @export
#'
#' @examples
#' seed <- 1
#' top_genes_matrix <- matrix
elbow_plot <- function(
    seed,                             # seed for start of clustering.
    top_genes_matrix,                 # input matrix for which the number of clusters is calculated.
    maxK = 15                         # Integer determining the maximum number of possible clusters, default = 15.
    ){
  maxK <- maxK
  set.seed(seed) # reproducible outcome
  wss <- function(k) {
    stats::kmeans(top_genes_matrix, algorithm = "Lloyd", k, nstart = 10, iter.max = 50)$tot.withinss
  }
  k.values <- 1:maxK # compute and plot wws for k = 1 to k = 15
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
#' m_kmeans <- matrix
#' number_of_annotations_per_cluster <- 10
#' k <- 1
most_variable_genes <- function(m_kmeans,number_of_annotations_per_cluster, k){
  last_column <- ncol(m_kmeans)
  top_x_variable_genes <- list()
  for (i in 1:k){
    cluster_matrix <- m_kmeans[(m_kmeans[,last_column]) == i,]
    cluster_matrix <- cluster_matrix[,-last_column]
    #estimates the variance for each row
    variance_row <- apply(cluster_matrix, 1, stats::var)

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

#' Function to set row annotation for a heatmap.
#'
#' @param m_top_genes_matrix An input matrix to used for heatmap generation and choosing annotation for a heatmap.
#' @param top_x_genes_cluster A list of the most variable x genes of each cluster in a heatmap.
#' @param fontsize_rowAnnotation An integer defining the font size of the row annotation
#'
#' @return A numeric index from the orginal matrix.
#' @export

set_annotation <- function(m_top_genes_matrix, top_x_genes_cluster, fontsize_rowAnnotation){
  #get numeric indices of top_x_genes_clusters
  top_x_genes_clusters <- which(BiocGenerics::rownames(m_top_genes_matrix) %in% top_x_genes_cluster)
  #set row annotation: at = numeric indices of wanted labels, labels = names of wanted labels, which = rows or columns
  labeling = BiocGenerics::rownames(m_top_genes_matrix)[top_x_genes_clusters]
  anno = ComplexHeatmap::anno_mark(at = top_x_genes_clusters, labels = labeling, labels_gp = grid::gpar(fontsize = fontsize_rowAnnotation), which = "row")

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
#' @param anno A numeric index containing row informations for annotation.
#' @param fontsize_columnNames An integer used to set the font size of the columns in the heatmap.
#' @param fontsize_rowNames An integer used to set the font size of the rownames in the heatmap.
#' @param title_heatmapLegend A string setting the title of the legend of the heatmap.
#' @param WidthNum A float setting the width of the heatmap.
#' @param HeightNum A float setting the height of the heatmap.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of HeightNum and WidthNum.
#' @param color_Palette Name of the colorPalette used for the heatmap.
#'
#' @return A plotted heatmap.
#' @export
#'
#' @examples
#' seed <- 1
#' m_top_genes_matrix <- matrix
#' title <- "Heatmap of Data"
#' split <- #split
#' fontsize_columnNames <- 6
#' fontsize_rowNames <- 4
#' title_heatmaoLegend <- "Expression"
#' WidthNum <- 4.5
#' HeightNum <- 3
#' UnitSize <- "cm"
#' color_Palette <- "RdBu"
print_heatmap <- function(seed,m_top_genes_matrix, title, split, anno, fontsize_columnNames, fontsize_rowNames, title_heatmapLegend,WidthNum, HeightNum, UnitSize, color_Palette){
  color_setting(color_Palette)
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

#' Creating a heatmap with annotation of x most variable rows(genes).
#'
#' @param topGenes_matrix An input matrix to create the heatmap.
#' @param probes A list of strings to summarize biological replicated in the columns.
#' @param number_of_annotations_per_cluster An integer used to set the number of annotations per cluster in the heatmap.
#' @param k An integer used to set number of cluster in the heatmap.
#' @param seed An integer used to set seed for reproducibility and k-Mean clustering.
#' @param Title A string to set the title of the heatmap.
#' @param fontsize_rowAnnotation An integer used to set the font size of the row annotation in the heatmap
#' @param fontsize_columnNames An integer used to set the font size of the columns in the heatmap.
#' @param fontsize_rowNames An integer used to set the font size of the rownames in the heatmap.
#' @param title_heatmapLegend A string setting the title of the legend of the heatmap.
#' @param WidthNum A float setting the width of the heatmap.
#' @param HeightNum A float setting the height of the heatmap.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of HeightNum and WidthNum.
#' @param color_Palette Name of the colorPalette used for the heatmap.
#'
#' @return A plotted heatmap of the input data.
#' @export
#'
#' @examples
#' topGenes_matrix <- matrix
#' probes <- list("Female", "Male")
#' number_of_annotations_per_cluster <- 10
#' k <- 3
#' seed <- 1
#' Title <- "Heatmap with annotation of most variable genes"
#' fontsize_rowAnnotation <- 8
#' color_Palette <- "RdBu"
function_complexHeatmap_var <- function(topGenes_matrix, probes, number_of_annotations_per_cluster, k, seed, Title, fontsize_rowAnnotation, fontsize_columnNames,fontsize_rowNames, title_heatmapLegend,WidthNum, HeightNum, UnitSize, color_Palette){
  m_top_genes_matrix <- summarise_bio_replicates(topGenes_matrix, probes)

  m_kmeans <- Kmean_generation(m_top_genes_matrix, seed, k)
  top_x_genes_cluster <- most_variable_genes(m_kmeans, number_of_annotations_per_cluster, k)
  anno <- DgeaHeatmap::set_annotation(m_top_genes_matrix, top_x_genes_cluster, fontsize_rowAnnotation)
  split <- performing_kMeans(m_top_genes_matrix, k)
  hm <- print_heatmap(seed,m_top_genes_matrix, Title, split, anno, fontsize_columnNames, fontsize_rowNames, title_heatmapLegend,WidthNum, HeightNum, UnitSize, color_Palette)
  return(hm)
}

#' Creating a color scheme based on the available color palettes of RColorBrewer for the heatmap.
#'
#' @param colorPalette An input string defining the color palette (available from RColorBrewer).
#'
#' @return The colors used in heatmap based on either ComplexHeatmap default or RColorBrewer.
#' @export
#'
#' @examples
#' colorPalette <- "RdBu"
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

#' Creating a color scheme based on the available color palettes of RColorBrewer for the heatmap.
#'
#' @param ncounts_matrix An input matrix to create the heatmap.
#' @param seed An integer used to set seed for reproducibility and k-Mean clustering, default = 1.
#' @param column_name A string to set the title of the heatmap, default = "Heatmap".
#' @param colorPalette Name of the colorPalette used for the heatmap, default = NULL.
#' @param cluster_method A string setting the cluster method for the heatmap, default = "hierarchical".
#' @param distance_method A string setting the distance method for clustering, default = "euclidean".
#' @param cluster_rows A Boolean switching the optional clustering of the rows on and off, default = TRUE.
#' @param cluster_columns A Boolean switching the optional clustering of the columns on and off, default = FALSE.
#' @param k_row An integer used to set number of clusters for row clustering in the heatmap, default = NULL.
#' @param k_col An integer used to set number of clusters for column clustering in the heatmap, default = NULL.
#' @param sample_metadata A dataframe containing the metadata information for the grouping of the columns, default = NULL.
#' @param annotation_colors A list assigning choosen colors to the corresponding groups, default = NULL.
#' @param annotation_name_side The side of the column annotation description, default = "right".
#' @param show_row_names A Boolean switching rownames on the heatmap on and off, default = FALSE.
#' @param show_column_names A Boolean switching colum names on the heatmap on and off, default = TRUE.
#' @param row_annotation A Boolean switching row annotation on the heatmap on and off, default = FALSE.
#' @param row_annotation_method A string setting the annotation method of the heatmap, default = "auto".
#' @param row_anno_names A list containing choosen rownames to use for the row annotation, default = NULL.
#' @param row_anno_number An integer setting the number of automatic annotations assigned per cluster, default = 5.
#' @param fontsize_title An integer setting the font size of the heatmap title, default = 15.
#' @param fontsize_rowAnnotation An integer setting the font size of the optional row annotation, default = 10.
#' @param fontsize_columnNames An integer setting the font size of the column names, default = 6.
#' @param fontsize_rowNames An integer setting the font size of the row names, default = 4.
#' @param fontsize_cluster_labels An integer setting the font size of the cluster labels, default = 8.
#' @param fontsize_group_annotation An integer setting the font size of the group annotation title, default = 8.
#' @param fontsize_group_annotation_legend An integer setting the font size of the group annotation legend name, default = 10.
#' @param fontsize_group_annotation_labels An integer setting the font size of the group annotation labels in the legend, default = 8.
#' @param fontsize_heatmap_legend An integer setting the font size of the heatmap legend title, default = 10.
#' @param fontsize_heatmap_legend_labels An integer setting the font size of the heatmap legend labels, default = 8.
#' @param title_heatmapLegend A string setting the changeable title of the legend, default "Expression".
#' @param WidthNum A float setting the width of the heatmap, default = 4.5.
#' @param HeightNum A float setting the height of the heatmap, default = 3.
#' @param UnitSize A string such as "cm" or "inch" to set the unit of HeightNum and WidthNum, default = "cm".
#'
#' @return An advanced and customizable heatmap.
#' @export
#'
#' @examples
#' colorPalette <- "RdBu"
adv_Heatmap <- function(
    ncounts_matrix,
    seed = 1,                                  # default seed = 1, changeable and suggested to be changed
    column_name = "Heatmap",                   # changeable title of the heatmap, default "Heatmap
    colorPalette = NULL,                       # "RdBu" ir any other color Palette available through RColorBrewer, default is blue to red over white from ComplexHeatmaps
    cluster_method = "hierarchical",           # "hierarchical", "kmeans"
    distance_method = "euclidean",             # or "correlation"
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    k_row = NULL,
    k_col = NULL,
    sample_metadata = NULL,
    annotation_colors = NULL,
    annotation_name_side = "right",            # side of annotation name, default = "right"
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_annotation = FALSE,                    # FALSE or TRUE
    row_annotation_method = "auto",            # "auto", "specific", "none"
    row_anno_names = NULL,                     # option to set list of specific genes as row annotation
    row_anno_number = 5,                       # default number of annotations per cluster
    fontsize_title = 15,                       # fontsize of the heatmap title
    fontsize_rowAnnotation = 10,               # fontsize of the optional row annotation
    fontsize_columnNames = 6,                  # fontsize of column names
    fontsize_rowNames = 4,                     # fontsize of row names
    fontsize_cluster_labels = 8,               # fontsize of cluster labels
    fontsize_group_annotation = 8,             # fontsize of optional group annotation
    fontsize_group_annotation_legend = 10,     # fontsize of optional group annotation legend title
    fontsize_group_annotation_labels = 8,      # fontsize of optional group annotation legend labels
    fontsize_heatmap_legend = 10,              # fontsize of heatmap legend
    fontsize_heatmap_legend_labels = 8,        # fontsize of heatmap legend labels
    title_heatmapLegend = "Expression",        # changeable title of the legend, default "Expression"
    WidthNum = 4.5,
    HeightNum = 3,
    UnitSize = "cm"
)
{
  # Distance matrix helper
  get_dist <- function(x, method) {
    if (method == "correlation") stats::as.dist(1 - stats::cor(t(x)))
    else stats::dist(x, method = method)
  }

  # Clustering logic
  row_split <- col_split <- NULL
  row_dend <- col_dend <- TRUE  # default TRUE if unspecified

  # --- ROW CLUSTERING ---
  if (cluster_rows) {
    row_data <- ncounts_matrix

    if (cluster_method == "hierarchical") {
      row_dist <- get_dist(row_data, distance_method)
      row_clust <- stats::hclust(row_dist, method = "complete")
      row_dend <- stats::as.dendrogram(row_clust)
      if (!is.null(k_row)) {
        row_split <- stats::cutree(row_clust, k = k_row)
        row_split <- as.factor(row_split)
        row_dend = TRUE
      }

    } else if (cluster_method == "kmeans") {
      row_split <- performing_kMeans(row_data, k_row)

    }
  }

  # --- COLUMN CLUSTERING ---
  if (cluster_columns) {
    col_data <- t(ncounts_matrix)

    if (cluster_method == "hierarchical") {
      col_dist <- get_dist(col_data, distance_method)
      col_clust <- stats::hclust(col_dist, method = "average")
      col_dend <- stats::as.dendrogram(col_clust)
      if (!is.null(k_col)) {
        col_split <- stats::cutree(col_clust, k = k_col)
        col_split <- as.factor(col_split)
        col_dend = TRUE
      }

    } else if (cluster_method == "kmeans") {
      col_split <- performing_kMeans(col_data, k_col)

    }
  }
  # --- Sample Annotation ---
  col_ha <- NULL
  if (!is.null(sample_metadata)) {
    col_ha <- ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(sample_metadata),
      col = as.list(annotation_colors),
      annotation_name_side = annotation_name_side,
      annotation_name_gp = grid::gpar(fontsize = fontsize_group_annotation),
      annotation_legend_param = list(
                                     title_gp = grid::gpar(fontsize = fontsize_group_annotation_legend),   # font size for legend title
                                     labels_gp = grid::gpar(fontsize = fontsize_group_annotation_labels)   # font size for legend labels
      )
    )
  }

  # --- Row Annotation ---
  annotation_for_rows <- NULL
  if(row_annotation) {
    if (row_annotation_method == "auto") {

      if (!is.null(k_row)) {
        m_kmeans <- Kmean_generation(ncounts_matrix, seed, k_row)
        top_x_genes_cluster <- most_variable_genes(m_kmeans, row_anno_number, k_row) #get most variable expressed genes of each cluster
        annotation_for_rows <- set_annotation(ncounts_matrix, top_x_genes_cluster, fontsize_rowAnnotation)

      } else {
        y <- ncounts_matrix
        top_x_genes_cluster <- list()
        #estimates the variance for each row
        variance_row <- apply(y, 1, stats::var)
        # binds a column with the estimated variance of each gene to clustermatrix
        cluster_matrix <- BiocGenerics::cbind(y, variance_row)
        #get number of last column
        number_of_last_column <- ncol(y)
        # orders the variances from highest to lowest
        o <- BiocGenerics::order(y[, number_of_last_column], decreasing = TRUE)
        # orders matrix according to order of variances ( highest to lowest)
        cluster_matrix <- y[o,]
        # makes a list of the gene names with the highest variance
        cluster_genes_highest_var <- as.list(BiocGenerics::rownames(y)[1:row_anno_number])
        # creates list with most_variable_genes of all cluster
        top_x_genes_cluster <- append (top_x_genes_cluster, cluster_genes_highest_var)
        # set annotation for rows
        annotation_for_rows <- set_annotation(ncounts_matrix, top_x_genes_cluster, fontsize_rowAnnotation)
      }

    } else if (row_annotation_method == "specific"){
      annotation_for_rows <- set_annotation(ncounts_matrix, row_anno_names, fontsize_rowAnnotation)
    }
  }



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
  heatmap_color_scheme <- get_heatmap_colors(colorPalette)            # sets the color Palette for the heatmap
  set.seed(seed)                          # sets a seed for randomizer to generate reproducable heatmaps


  ht = ComplexHeatmap::Heatmap(
    ncounts_matrix,
    cluster_rows = if (is.logical(row_dend)) row_dend else row_dend,
    cluster_columns = if (is.logical(col_dend)) col_dend else col_dend,
    col = heatmap_color_scheme,
    row_split = row_split,
    column_split = col_split,
    top_annotation = col_ha,
    #right_annotation = right_anno,
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
    height = grid::unit(HeightNum, UnitSize))

  # Wrap addition of row annotation
  #right_anno <- NULL
  if (isTRUE(row_annotation)) {
    ht <- ht + ComplexHeatmap::rowAnnotation(mark = annotation_for_rows)
    #right_anno <- ComplexHeatmap::rowAnnotation(
    #mark = annotation_for_rows,
    #annotation_name_gp = grid::gpar(fontsize = fontsize_rowAnnotation)

  }

  hm = ComplexHeatmap::draw(ht)
  return(hm)


}
