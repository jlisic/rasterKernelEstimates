#' Local quantiles for an in memory raster image
#' 
#' \code{rasterLocalQuantiles} finds the quantile within the positive valued neighborhood 
#' of \code{W}.
#'
#' @param r An in memory raster image.
#' @param W A matrix of weights used to specify a local neighborhood.  The quantile 
#'   kernel will be applied to each pixel in \code{r}.  Dimensions must be non-zero 
#'  and odd.  
#' @param q A quantile.  This value is required to be in the inclusive interval from 
#'   0 to 100.
#' @details A spatial neighborhood is calculated for each pixel in \code{r}.
#'   The spatial neighborhood for each pixel is defined by the weight matrix
#'   \code{W}, where the center of the odd dimensioned matrix \code{W} is identified 
#'   with the target pixel.  The target pixel value is replaced with the
#'   quantile of the neighborhood identified by \code{W}.  Only non-missing or neighbors
#'   with non-zero weights are used in the calculation.  Quantile calculation uses
#'   the inverse empirical CDF transform, equivalent to \code{stats::quantile} type=1.
#' @return An in memory raster image of local quantiles.
#' @examples 
#' r <- raster::raster( matrix(rnorm(36),6,6)) 
#' W <- matrix(1,3,3)
#' medianR <- rasterLocalQuantiles(r,W)
#' @importFrom raster raster
#' @importFrom raster values
#' @useDynLib rasterKernelEstimates
#' @export
rasterLocalQuantiles <-
function(
  r,
  W,
  q=50
  ) {

  if( q < 0 | q > 100 )  stop("Invalid quantile value q.")
  
  r.values <- raster::values(r)

  # get rid of na
  r.values.na <- is.na(r.values) 

  r.values[r.values.na]  <- Inf 

  r.result <- .C("rSmoothLocalQuantile",
    as.double(r.values),
    as.double(r.values),
    as.double(c(t(W))),
    as.double(q/100),
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
