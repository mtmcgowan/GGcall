#' marker.create
#'
#' Creates an empty Marker object
#'
#' @param name the marker ID
#' @param x a vector of x.values
#' @param y a vector of y.values
#' @param x.scaled a vector of scaled x.values
#' @param y.scaled a vector of scaled y.values
#' @param R a vector of R values (either from GenomeStudio or calculated from raw values)
#' @param theta a vector of theta values (either from GenomeStudio or calculated from raw values)
#' @param theta.trans a vector of arcsin transformed values
#' @param mixmodout cluster results from mixmodCluster
#' @param miss.ind a vecor of missing value indices
#' @param clust.time timing results for the clustering
#' @param data.type Which data points were used for clustering
#'
#' @return An object storing all data for a single marker

marker.create <- function(name = 'NA',
                              x = vector(),
                              y = vector(),
                              x.scaled = vector(),
                              y.scaled = vector(),
                              R = vector(),
                              theta = vector(),
                              theta.trans = vector(),
                              mixmodout = NULL,
                              miss.ind = NULL,
                              clust.time = NULL,
                              data.type = NULL)
{
  marker = list(name = name,
                data = list(x = x,
                            y = y,
                            x.scaled = x.scaled,
                            y.scaled = y.scaled,
                            theta = theta,
                            theta.trans = theta.trans,
                            r = R,
                            miss.ind = miss.ind,
                            filt.flag = 0),
                model = list(mixmodout = mixmodout, data.type = NA, clust.time = clust.time),
                genotypes = NULL
  )

  return (marker)
}
