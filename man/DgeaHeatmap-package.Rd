\name{DgeaHeatmap-package}
\alias{DgeaHeatmap-package}
\alias{DgeaHeatmap}
\docType{package}
\title{DgeaHeatmap: Tools for easier Heatmap generation and DGEA}
\description{

Collection of tools to generate heatmaps and to analyze differentially expressed genes.
}

\section{Package options}{


Devtools uses the following functions to simplify heatmap generation and DGEA for users:
\itemize{
\item \code{Kmean_generation}: generates K-means for the columns and rows of a matrix.
\item\code{aExprsDataQC}: function for automatized quality control.
\item\code{add_demoElem}: builds a matrix from an input NanostringGeoMxSet object.
\item\code{build_matrix}: builds a matrix from an input csv file.
\item\code{color_setting}: sets a color scheme for heatmaps.
\item\code{elbow_plot}:creates an elbow plot to help choose k for clustering by k-Means.
\item\code{filtering_for_top_exprGenes}: filters a matrix to extract a certain number of most variable rows through calculation of their variance.
\item\code{function_complexHeatmap_var}: creates a heatmap with automatic annotation.
\item\code{genRawReadCountTable}: generates a raw read count table from an NanostringGeoMxSet object.
\item\code{individual_matrix}: creates a matrix containing only specific columns extracted from an orginial matrix.
\item\code{most_variable_genes}: determines the most variable genes of each cluster to enable annotation.
\item\code{performing_kMeans}: performs k-Means clustering for a matrix and sets split in order to split a heatmap into clusters.
\item\code{print_heatmap}: builds a heatmap.
\item\code{scale_counts}: scales the values of a matrix by Z-count.
\item\code{set_annotation}: sets annotation for a heatmap.
\item\code{show_data_distribution}: visualizes the data distribution within the data of a matrix.
\item\code{split_data_by_column}: splits data by a group column with feature, pheno or protocol data to then get then get the mean.
\item\code{summarise_bio_replicates}: summarizes columns biological replicates of a matrix into one.
}
}
\details{
%%  ~~ An overview of how to use the package, ~~
%%  ~~ including the most important functions ~~
}
\author{
Leonie Lancelle, Medical School OWL, Bielefeld University

Maintainer: Leonie Lancelle, Medical School OWL, Bielefeld University
}
\references{
%%  ~~ Literature or other references for background information ~~
}
\keyword{package}

%% Uncomment below to imitate parts of library(help = DgeaHeatmap)
%\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{DgeaHeatmap}}
%\section{Documentation Index}{\packageIndices{DgeaHeatmap}}

\seealso{
%%  ~~ Optional links to other man pages, e.g. ~~
%%  ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%%  ~~ Optional simple examples of the most important functions ~~
}
