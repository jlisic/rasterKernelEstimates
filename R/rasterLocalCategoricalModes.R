#' Local categorical modes for an in memory raster image
#' 
#' \code{rasterLocalCategoricalModes} finds the most popular category within the 
#' weighted neighborhood of \code{W}.
#'
#' @param r An in memory raster image. Pixels should be whole numbers or \code{NA}. 
#'   Pixels with non-whole number values will be coerced into whole numbers.
#' @param W A matrix of weights.  The modal kernel will be applied to each 
#'   pixel in \code{r}.  Dimensions must be non-zero and odd.
#' @return An in memory raster image by most popular categories.
#' @importFrom raster raster
#' @importFrom raster values
#' @details A spatial neighborhood is calculated for each pixel in \code{r}.
#'   The spatial neighborhood for each pixel is defined by the weight matrix
#'   \code{W}, where the center of the odd dimensioned matrix \code{W} is identified 
#'   with the target pixel.  The target pixel value is replaced with the most
#'   popular value within the neighborhood weighted by \code{W}.  Ties are
#'   handled by randomly by uniformly selecting a category amongst the tied
#'   categories.  Only non-missing or neighbors with non-zero weights are used 
#'   in the calculation.
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
    r.values <- r.values - min.offset
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
  
  if( min.offset < 0 ) r.result <- r.result + min.offset
  
  raster::values(r) <- r.result 

  return( r  )
}
