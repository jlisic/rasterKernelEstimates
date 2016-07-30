#' Local medians for an in memory raster image.
#' 
#' \code{rasterLocalMedians} finds the median within the positive valued neighborhood 
#' of W.
#'
#' @param r An in memory raster image.
#' @param W A matrix of weights.  The median kernel will be applied at each 
#'   centroid.  Dimensions must be non-zero and odd.  Only non-missing or 
#'   neighbors with non-zero weights are used in the median calculation.
#' @return An in memory raster image of local medians.
#' @examples 
#' r <- raster::raster( matrix(rnorm(36),6,6)) 
#' W <- matrix(1,3,3)
#' medianR <- rasterLocalMedians(r,W)
#' @importFrom raster raster
#' @importFrom raster values
#' @useDynLib rasterKernelEstimates
#' @export
rasterLocalMedians <-
function(
  r,
  W
  ) {

  r.values <- raster::values(r)

  # get rid of na
  r.values.na <- is.na(r.values) 

  r.values[r.values.na]  <- Inf 

  r.result <- .C("rSmoothLocalMedian",
    as.double(r.values),
    as.double(r.values),
    as.double(c(t(W))), 
    as.integer(nrow(r)),
    as.integer(ncol(r)),
    as.integer(nrow(W)),
    as.integer(ncol(W)),
    NAOK=TRUE,
    PACKAGE='rasterKernelEstimates'
  )

  r.mu <- r

  r.result[[2]][r.values.na] <- NA
  raster::values(r.mu) <- r.result[[2]] 

  return( r.mu )
}
