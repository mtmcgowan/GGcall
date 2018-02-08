#' marker.cluster
#'
#' Clusters a single markers with Rmixmod
#'
#' USes the 'Rmixmod' package to fit a Gaussian mixture model using a supplied model list
#'
#' @param marker A project object, must be formatted correctly according to the project.create function
#' @param data.type A string variable indicating which data should be used for clustering. ('raw', 'GSnorm','GStrans', 'thetaonly')
#' @param iteration The number of times to re-cluster (can result in more stable results, but slows clustering time)
#' @param model.list A list of Rmixmod models to use for cluster (Refer to Rmixmod documentation for more info)
#' @param clust.num The maximum number of clusters to test for (lowering this speeds clustering time)
#'
#' @return A marker object updated with the best clustering result out of a defined number of iterations
#'
#' This function relies on the following packages:
#' Rmixmod

marker.cluster <- function(marker,
                           data.type = 'raw',
                           iteration = 5,
                           model.list = c("Gaussian_pk_Lk_Ck", "Gaussian_pk_Lk_Bk"),
                           clust.num = 8)
{
  # Check if the marker should be clustered (bad markers are flagged)
  if (marker$data$filt.flag == 1)
  {
    return(marker)
  }

  # Initialize the clustering conditions
  clust.strategy <- mixmodStrategy(algo = c('EM', 'CEM'), nbTry = iteration,
                                   initMethod = "CEM",
                                   nbTryInInit = 1000, nbIterationInInit = 5,
                                   nbIterationInAlgo = 500, epsilonInInit = 0.001,
                                   epsilonInAlgo = 0.001, seed = NULL, parameter=NA,
                                   labels=NA)


  ############# SETTING UP THE DATA #############
  if (data.type == 'raw')
  {
    # Use raw x and y values as clustering dimensions
    x = marker$data$x
    y = marker$data$y
  }

  else if (data.type == 'raw.scaled')
  {
    # use scaled x and y values as clustering dimensions
    x = marker$data$x.scaled
    y = marker$data$y.scaled
  }
  else if (data.type == 'GS.norm')
  {
    x = marker$data$theta
    y = marker$data$r
  }
  else if (data.type == 'GS.trans')
  {
    x = marker$data$theta.trans
    y = marker$data$r
  }
  else
  {
    print("Invalid Data Type")
    return(marker)
  }

  # Determine the indices that are NA values
  x.na = which(is.na(x))
  y.na = which(is.na(y))
  all.na = unique(c(x.na, y.na))

  # Add NA indices to the marker object (to remember which positions to insert for classifications..)
  marker$data$miss.ind = all.na

  # Subset only non-NA values
  if (length(all.na) > 0)
  {
    x = x[-all.na]
    y = y[-all.na]
  }


  # Combine vectors into a data.frame
  data = data.frame(x, y)


  ############# CLUSTERING THE DATA ###############

  # Assign the data used for clustering
  marker$model$data.type = data.type

  # Initialize time variable
  ptm <- proc.time()

  # Perform clustering
  clust_results <- mixmodCluster(data,
                                 nbCluster = 1:clust.num,
                                 models = mixmodGaussianModel(listModels = model.list),
                                 criterion = 'ICL',
                                 strategy = clust.strategy)

  # Record clustering time
  marker$model$clust.time <- proc.time() - ptm

  # Assign the clustering results to the marker
  marker$model$mixmodout <- clust_results['bestResult']

  # Assign the data type used for clustering
  marker$model$data.type <- data.type

    # Return the marker with updated results
    return (marker)
}
