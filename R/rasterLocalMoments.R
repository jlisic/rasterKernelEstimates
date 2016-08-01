#' Local moments for an in memory raster image
#' 
#' \code{rasterLocalMoments} finds the local moments within the weighted neighborhood 
#' of W.
#'
#' @param r An in memory raster image.
#' @param WMu A matrix of weights.  The mean kernel will be applied to each 
#'   pixel in \code{r}.  Dimensions must be non-zero and odd.  Only non-missing 
#'   neighbors are used in the mean.  
#' @param WVar A matrix of weights.  The variance kernel will be applied at each 
#'   centroid.  Dimensions must be non-zero and odd.  Only non-missing 
#'   neighbors are used in the variance. The dimensions of WVar must match WMu.
#' @param moments The number of moments to calculate. The local spatial mean
#'   will be calculated when moments=1.  The local spatial mean and variance
#'   wil be calculated when moments=2.  Currently no higher moments are supported.  
#' @return A list of in memory raster images, one list element for each moment.
#' @examples 
#' r <- raster::raster( matrix(rnorm(36),6,6)) 
#' W <- matrix(1,3,3)
#' rLocalMoments <- rasterLocalMoments(r,W)
#' @importFrom raster raster
#' @importFrom raster values
#' @useDynLib rasterKernelEstimates
#' @export
rasterLocalMoments <-
function(
  r,
  WMu,
  WVar=WMu,
  moments=2
  ) {

  
  r.values <- raster::values(r)

  # get rid of na
  r.values.na <- is.na(r.values) 

  r.values[r.values.na]  <- Inf 

  # we only do moments > 2, and I'm not sure what moment = 0 or smaller are
  if(!(moments %in% 1:2)) stop(sprintf("moments = %d is not a valid moment", moments))

  if(nrow(WMu) != nrow(WVar)) stop("WMu and WVar must be the same size")
  if(ncol(WMu) != ncol(WVar)) stop("WMu and WVar must be the same size")

  r.result <- .C("rSmoothLocalMoments",
    as.double(r.values),
    as.double(r.values),
    as.double(r.values), 
    as.double(c(t(WMu))), 
    as.double(c(t(WVar))), 
    as.integer(nrow(r)),
    as.integer(ncol(r)),
    as.integer(nrow(WMu)),
    as.integer(ncol(WMu)),
    as.integer(moments)
    , NAOK=TRUE,
    PACKAGE='rasterKernelEstimates'
  )

  r.mu <- r
  if( moments > 1 ) r.var <- r
  

  r.result[[2]][r.values.na] <- NA
  if( moments > 1 ) r.result[[3]][r.values.na] <- NA

  raster::values(r.mu) <- r.result[[2]] 
  if( moments > 1) raster::values(r.var) <- r.result[[3]] 

  if( moments > 1) return( list( mu=r.mu, var=r.var ) )

  return( r.mu )
}
