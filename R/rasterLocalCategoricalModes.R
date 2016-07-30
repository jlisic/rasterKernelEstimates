#' Local categorical modes for an in memory raster image.
#' 
#' \code{rasterLocalCategoricalModes} finds the most popular category within the 
#' weighted neighborhood of W.
#'
#' @param r An in memory raster image. Pixels should be whole numbers or \code{NA}. 
#'   Pixels with non whole number values will be coerced to whole numbers.
#' @param W A matrix of weights.  The modal kernel will be applied at each 
#'   centroid.  Dimensions must be non-zero and odd.
#' @return An in memory raster image by most popular categories.
#' @importFrom raster raster
#' @importFrom raster values
#' @examples 
#' r <- raster::raster( matrix(runif(81),9,9)) 
#' W <- matrix(1,3,3)
#' modeR <- rasterLocalCategoricalModes(r,W)
#' @useDynLib rasterKernelEstimates
#' @export
rasterLocalCategoricalModes <-
function(
  r,
  W
  ) {

  r.values <- round(raster::values(r))

  # get rid of na
  r.values.na <- is.na(r.values)
  
  # find the min category
  min.offset <- min(r.values,na.rm=T)  

  if( min.offset > 0 ) {
    r.values[r.values.na]  <- -1 
  } else {
    r.values <- r.values + min.offset
    r.values[r.values.na]  <- -1 
  }

  r.result <- .C("rSmoothCategorical",
    as.integer(r.values),
    as.integer(r.values), 
    as.double(W),
    as.integer(nrow(r)),
    as.integer(ncol(r)),
    as.integer(nrow(W)),
    as.integer(ncol(W)),
    PACKAGE='rasterKernelEstimates'
  )

  # copy over values
  r.result <- r.result[[2]]
  
  # retain NA
  r.result[r.values.na] <- NA
  
  if( min.offset <= 0 ) r.result <- r.result + min.offset
  
  raster::values(r) <- r.result 

  return( r  )
}
