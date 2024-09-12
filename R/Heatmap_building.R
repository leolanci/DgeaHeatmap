build_matrix <- function(counts_data){
  BiocGenerics::colnames(counts_data)[1] <- "genes"
  list_rownames <- as.list(counts_data[,1])
  print(list_rownames)
  counts_data <- as.matrix(counts_data[,-1])
  rownames(counts_data) <-list_rownames$genes
  return(counts_data)
}

individual_matrix <- function(factors_for_matrix_devision, mmatrix) {
  for (i in factors_for_matrix_devision) {
    chosen_columns <- BiocGenerics::grep(i, BiocGenerics::colnames(mmatrix))
    rmatrix <- BiocGenerics::cbind(mmatrix[, chosen_columns])
    mmatrix <- rmatrix
    rmatrix <- matrix()
  }
  return(mmatrix)
}

filtering_for_top_exprGenes <- function(counts_data, top_number_of_genes){
  var_genes <- apply(counts_data, 1, var) #estimating the variance of each gene
  # Sorting the genes by their variance and creating a new object with chosen number of most variable genes
  select_var <- names(sort(var_genes, decreasing = TRUE))[1:top_number_of_genes]
  anaMatrix <- counts_data
  highly_variable_genes <- anaMatrix[select_var,]
  dim(highly_variable_genes)

  return(highly_variable_genes)
}

scale_counts <- function(countsmatrix){
  scaled_counts =
    countsmatrix%>%
    t(.) %>% #transpose to have genes in columns
    scale() %>%  #scale (x, center = TRUE, scale = TRUE)
    t(.)  #'transpose back in original shape
  return(scaled_counts)
}

show_data_distribution <- function(scaled_counts){
  apply(scaled_counts, MARGIN = 1, mean) %>%  #calculate the mean per row
    graphics::hist(., main = "", xlab = "Z-score values", col = "dodgerblue2") #build histogram to see data distribution
}

elbow_plot <- function(seed, top_genes_matrix){
  set.seed(seed) # reproducible outcome
  wss <- function(k) {
    kmeans(top_genes_matrix, k, nstart = 10)$tot.withinss
  }
  k.values <- 1:15 # compute and plot wws for k = 1 to k = 15
  wss_values <- purrr::map_dbl(k.values, wss)
  plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")

}

color_setting <- function(colorPalette){
  my_colors = RColorBrewer::brewer.pal(n=11, name = colorPalette)
  my_colors = RColorBrewer::colorRampPalette(my_colors)(50)
  my_colors = rev(my_colors)
  my_colors
}


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


set_annotation <- function(m_top_genes_matrix, top_x_genes_cluster){
  #get numeric indices of top_x_genes_clusters
  top_x_genes_clusters <- which(BiocGenerics::rownames(m_top_genes_matrix) %in% top_x_genes_cluster)
  #set row annotation: at = numeric indices of wanted labels, labels = names of wanted labels, which = rows or columns
  labeling = BiocGenerics::rownames(m_top_genes_matrix)[top_x_genes_clusters]
  anno = ComplexHeatmap::anno_mark(at = top_x_genes_clusters, labels = labeling, which = "row")

  return(anno)
}


performing_kMeans <- function(m_top_genes_matrix, k){
  kclus <- stats::kmeans(m_top_genes_matrix, k)  #performance of kmean clustering
  kclus$cluster #check on kmean clusters
  split <- paste0("Cluster\n", kclus$cluster) #set split as kmean clusters for heatmap
  return(split)
}


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

