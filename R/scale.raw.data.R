scale.raw.data <- function(x, y) {
  x.scale = x
  y.scale = y

  for (i in 1:nrow(x)) {
    #cat('iteration ', i, '\n', sep = '')
    x.temp <- unlist(x[i,1:ncol(x)])
    y.temp <- unlist(y[i,1:ncol(y)])
    theta.rad <- rad.theta.conv(x.temp, y.temp)
    r <- r.convert(x.temp, y.temp)

    rnorm <- scale(r, center = F, scale = T)
    x.scale[i,] <- rnorm * sin(theta.rad)
    y.scale[i,] <- rnorm * cos(theta.rad)

  }
  return(list(x.scale, y.scale))
}

rad.theta.conv <- function(x, y) {return(atan(y/x))}
r.convert <- function(x, y) {return(sqrt(y^2 + x^2))}
