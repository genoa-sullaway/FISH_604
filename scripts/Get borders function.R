## Function 'get.borders' to determine convex hull for a set of points 
## written by Margaret Short
get.borders <- function(my.lonlat, frac=0.15 )
{
  if( !is.matrix(my.lonlat) )
    stop("get.borders expects a matrix (w/lon in 1st col, lat in 2nd)")
  
  if( ncol(my.lonlat) != 2 )
    stop(paste("get.borders expects a matrix with",
               " exactly 2 columns, not ",ncol(my.lonlat),"columns"))
  
  lons <- my.lonlat[,1]
  n.distinct.lons <- length(unique(lons))
  
  lats <- my.lonlat[,2]
  n.distinct.lats <- length(unique(lats))
  
  minlon <- min(lons)
  maxlon <- max(lons)
  minlat <- min(lats)
  maxlat <- max(lats)
  
  if( n.distinct.lons == 1 ) {
    dd <- (maxlat - minlat) * frac
    my.borders <- cbind(
      lons[1] + dd*c(-1,1,1,-1),
      c(minlat-dd,minlat-dd,maxlat+dd,maxlat+dd) )
  } else if (n.distinct.lats == 1 ) {
    dd <- (maxlon - minlon) * frac
    my.borders <- cbind(
      c(minlon-dd,maxlon+dd,maxlon+dd,minlon-dd),
      lats[1] + dd*c(-1,-1,1,1))
  } else {
    
    tmp <- chull( my.lonlat ) # list of vtx numbers
    xx <- my.lonlat[tmp,1]
    yy <- my.lonlat[tmp,2]
    init.borders <- cbind( xx,yy )
    
    # find the (approximate) centroid by overlaying
    # a fairly dense grid of points; the mean of the
    # x-coordinates is approx. the x-coord of centroid,etc.
    
    # 1000 points = dlon * dlat / dxdy^2
    dxdy <- sqrt( (maxlon-minlon)*(maxlat-minlat)/1000 )
    
    grid.x <- seq( from=minlon, to=maxlon, by=dxdy )
    grid.y <- seq( from=minlat, to=maxlat, by=dxdy )
    my.grid <- as.matrix(expand.grid( x=grid.x,y=grid.y ))
    
    my.gr.in <- locations.inside(my.grid, init.borders)
    cx <- mean( my.gr.in[,1] ) # approx. x-coord of centroid
    cy <- mean( my.gr.in[,2] )
    xx <- (xx-cx)*(1.0+frac) + cx
    yy <- (yy-cy)*(1.0+frac) + cy
    my.borders <- cbind(xx,yy)
  }
  return( my.borders )
}
