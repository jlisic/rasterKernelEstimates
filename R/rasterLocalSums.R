#' Local sums for an in memory raster image.
#' 
#' \code{rasterLocalSums} finds the local sum within the weighted neighborhood of W.
#'
#' @param r An in memory raster image.
#' @param W A matrix of weights.  The sums will be applied at each centroid. 
#'   Dimensions must be non-zero and odd.  Only non-missing neighbors are used in 
#'   the sum.
#' @details A spatial neighborhood is calculated for each pixel in \code{r}.
#'   The spatial neighborhood for each pixel is defined by the weight matrix
#'   \code{W}, where the center of the odd dimensioned matrix \code{W} is identified 
#'   with the target pixel.  The target pixel value is replaced with the sum of
#'   all pixels within the neighborhood weighted by \code{W}.   Only non-missing 
#'   or neighbors with non-zero weights are used in the calculation.
#' @return An in memory raster image of local sums.
#' @examples 
#' r <- raster::raster( matrix(rnorm(36),6,6)) 
#' W <- matrix(1,3,3)
#' sumR <- rasterLocalSums(r,W)
#' @importFrom raster raster
#' @importFrom raster values
#' @useDynLib rasterKernelEstimates
#' @export
rasterLocalSums <-
function(
  r,
  W
  ) {
  
  r.values <- raster::values(r)

  # get rid of na
  r.values.na <- is.na(r.values) 

  r.values[r.values.na]  <- Inf 


  r.result <- .C("rSmoothSums",
    as.double(r.values),
    as.double(r.values),
    as.double(c(t(W))), 
    as.integer(nrow(r)),
    as.integer(ncol(r)),
    as.integer(nrow(W)),
    as.integer(ncol(W)),
    PACKAGE='rasterKernelEstimates'
  )

  r.result <<- r.result

  r.mu <- r

  r.result[[2]][r.values.na] <- NA
  raster::values(r.mu) <- r.result[[2]] 

  return( r.mu )
}
