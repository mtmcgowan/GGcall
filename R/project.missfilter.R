#' project.missfilter
#'
#' Checks for problematic samples and markers. Samples with a set missing rate are flagged. Markers with zero-variance (either all theta or R values are the same) are also flagged for removal.
#'
#' @param project a project object from imported GenomeStudio data
#' @param missrate The missing rate threshold where markers with a proportion of NA values above this threshold will be removed (default = 0.2 i.e. 20\% cutoff)
#' @param type Which type of data should be check. GenomeStudio theta/R values will have more missing values than the raw X/Y. (either 'raw' or 'GSnorm')
#'
#' @return A project object with filter flags applied to markers


project.missfilter <- function(project, missrate = 0.2, type = 'raw') {
  calc.missrate <- function(x)
  {
    return(sum(is.na(x)) / length(x))
  }

  if (type == 'raw')
  {
    # Calculate 'x' and 'y' missing rates for MARKERS
    x.missrate <- apply(project$data$x, 1, calc.missrate)
    y.missrate <- apply(project$data$y, 1, calc.missrate)
    x.scaled.missrate <- apply(project$data$x.scaled, 1, calc.missrate)
    y.scaled.missrate <- apply(project$data$x.scaled, 1, calc.missrate)

    # Calculate variance for 'theta' and 'r' variables for MARKERS
    x.var <- apply(project$data$x, 1, function(x) {var(x, na.rm = T)})
    y.var <- apply(project$data$y, 1, function(x) {var(x, na.rm = T)})
    x.scaled.var <- apply(project$data$x.scaled, 1, function(x) {var(x, na.rm = T)})
    y.scaled.var <- apply(project$data$y.scaled, 1, function(x) {var(x, na.rm = T)})

    # Determine which MARKERS should be removed
    rm.x <- which(x.missrate  > missrate | x.var == 0)
    rm.y <- which(y.missrate > missrate | y.var == 0)
    rm.x.scaled <- which(x.scaled.missrate > missrate | x.scaled.var == 0)
    rm.y.scaled <- which(y.scaled.missrate > missrate | y.scaled.var == 0)
    rm.all <- unique(c(rm.x, rm.y, rm.x.scaled, rm.y.scaled))

    if (length(rm.all) > 0)
    {
      for (i in rm.all)
      {
        project$markers[[i]]$filt.flag <- 1
      }
    }
  }

  if (type == 'GSnorm')
  {
    # Calculate 'x' and 'y' missing rates for MARKERS
    theta.missrate <- apply(project$data$theta, 1, calc.missrate)
    r.missrate <- apply(project$data$r, 1, calc.missrate)
    theta.trans.missrate <- apply(project$data$theta.trans, 1, calc.missrate)

    # Calculate variance for 'theta' and 'r' variables for MARKERS
    theta.var <- apply(project$data$theta, 1, function(x) {var(x, na.rm = T)})
    r.var <- apply(project$data$r, 1, function(x) {var(x, na.rm = T)})
    theta.trans.var <- apply(project$data$theta.trans, 1, function(x) {var(x, na.rm = T)})

    # Determine which MARKERS should be removed
    rm.theta <- which(theta.missrate  > missrate | theta.var == 0)
    rm.r <- which(r.missrate > missrate | r.var == 0)
    rm.theta.trans <- which(theta.trans.missrate > missrate | theta.trans.var == 0)
    rm.all <- unique(c(rm.theta, rm.r, rm.theta.trans))

    if (length(rm.all) > 0)
    {
      for (i in rm.all)
      {
        project$markers[[i]]$data$filt.flag <- 1
      }
    }
  }
  return(project)
}
