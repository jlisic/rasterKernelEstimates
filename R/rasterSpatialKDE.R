#' Calculates a local Gaussian kernel estimate for an in memory raster image.
#' 
#' \code{rasterLocalMoments} calculates a kernel density estimate truncated by 
#' the weighted neighborhood of W.
#'
#' @param r An in memory raster image.
#' @param h The bandwidth parameter for the raster values. 
#' @param h.x An optional truncation parameter for the spatial support. This
#'   parameter is specified in number of horizontal pixels to be evaluated.
#'   The truncated neighborhood is centered around the current pixel, requiring
#'   this parameter to be odd.  If the parameter is not odd, it will be 
#'   increased in size by one. The default parameter is the minimum odd value
#'   of the number of columns, or the number of columns plus one.
#' @param h.y An optional truncation parameter for the Gaussian kernel. This
#'   parameter is specified in number of vertical pixels to be evaluated.
#'   The truncated neighborhood is centered around the current pixel, requiring
#'   this parameter to be odd. If the parameter is not odd it will be 
#'   increased in size by one.  The default parameter is the minimum odd value
#'   of number or rows, or the number of rows plus one.
#' @return A raster image of the kernel density estimate.
#' @examples 
#' r <- raster::raster( matrix(rnorm(36),6,6)) 
#' W <- matrix(1,3,3)
#' rLocalKDE <- rasterSpatialKDE(r,W)
#' @importFrom raster raster
#' @importFrom raster values
#' @useDynLib rasterKernelEstimates
#' @export
rasterSpatialKDE <-
function(
  r,
  h,
  h.x=ncol(r),
  h.y=nrow(r) 
  ) {

  # convert to pixels
  h.x <- as.integer(h.x)
  h.y <- as.integer(h.y)
  
  #convert to integers
  h.x <- ifelse(as.integer(h.x) %% 2 == 0, h.x, h.x+1)
  h.y <- ifelse(as.integer(h.y) %% 2 == 0, h.y, h.y+1)

  # bandwidth checks
  if( h.x < 1 ) stop( "h.x < 1")
  if( h.y < 1 ) stop( "h.y < 1")
  if( h <= 0 ) stop("bandwidth must be greater than 0")

  # get bandwidth value
  
  r.values <- raster::values(r)

  # get rid of na
  r.values.na <- is.na(r.values) 

  r.values[r.values.na]  <- -1 

  r.result <- .C("rSpatialKDE",
    as.double(r.values),
    as.double(r.values), 
    as.double(h),
    as.integer(nrow(r)),
    as.integer(ncol(r)),
    h.y,
    h.x,
    PACKAGE='rasterKernelEstimates'
  )


  r.mu <- r
  r.result[[2]][r.values.na] <- NA

  raster::values(r.mu) <- r.result[[2]] 

  return( r.mu  )
}
