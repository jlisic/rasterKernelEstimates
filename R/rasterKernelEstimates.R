#' rasterKernelEstimates:  Kernel Based Estimates on In-Memory Raster Images
#'
#'
#' rasterKernelEstimats provides kernel based estimates for in-memory raster images 
#' from the raster package.  These kernel estimates include local means
#' variances, modes, medians, and Gaussian kernel density estimates.  All
#' results are in the form of raster images, preserving original resolution
#' and projection attributes in the raster image.
#'
#' It has two main goals:
#'
#' \itemize{
#' \item Provide a method to quickly perform common actions on in-memory raster
#'   images. 
#' \item Provide fast performance for in-memory raster images through C and OpenMP.
#' }
#'
#'
#' @docType package
#' @name rasterKernelEstimates
#' @import raster
#' @importFrom raster raster
#' @importFrom raster values
NULL