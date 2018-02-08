#' project.cluster
#'
#' Cluster markers in Parallel with Rmixmod
#'
#' Sets up a parallel backend and clusters each marker using mixture of gaussian models provided by the 'Rmixmod' package.
#'
#' @param project A project object, must be formatted correctly according to the project.create function
#' @param data.type A string variable indicating which data should be used for clustering. ('raw', 'GSnorm','GStrans', 'thetaonly')
#' @param iteration The number of times to re-cluster (can result in more stable results, but increases clustering time)
#' @param model.list A list of Rmixmod models to use for cluster (Refer to Rmixmod documentation for more info)
#' @param clust.num The maximum number of clusters to test for (lowering this speeds clustering time)
#' @return A project object updated with clustering models
#'
#' This function relies on the following packages:
# parallel, pbapply

project.cluster <- function(project,
                            marker_indices = NULL,
                            data.type = 'raw.scaled',
                            iteration = 5,
                            model.list = c("Gaussian_pk_L_Ck", "Gaussian_pk_L_Bk"),
                            clust.num = 8)
{
  if (is.null(marker_indices))
  {
    marker_indices = 1:length(project$markers)
  }
  ##### Set up the parallel environment #####
  print('Setting up the parallel environment:')
  # Initialize parallel cores
  no.cores <- detectCores() - 1
  cl <- makeCluster(no.cores, outfile = "")

  print('Done')
  print('Exporting data and functions to clusters:')

  # Export required objects to each logic
  clusterExport(cl,
                list("marker.cluster",
                     "mixmodCluster",
                     'mixmodGaussianModel',
                     'mixmodStrategy')
  )
  print('Done')

  # Set up the progress bar
  do.pb <- dopb()
  split <- splitpb(length(marker_indices), length(cl), nout = 100)
  B <- length(split)
  if (do.pb) {
    pb <- startpb(0, B)
    on.exit(closepb(pb), add = TRUE)
  }

  # Cluster the samples
  print('Clustering samples:')

  # Perform parallel clustering in blocks, updating the progress bar with each block
  par.out <- vector("list", B)
  for (i in seq_len(B)) {
    if(.Platform$OS.type == "windows")
    {
      par.out[i] <- list(parallel::parLapply(cl,
                                        project$markers[split[[i]]],
                                        function(x) marker.cluster(x,
                                                                   data.type,
                                                                   iteration,
                                                                   model.list,
                                                                   clust.num)))
    }
    else if (.Platform$OS.type == "unix")
    {
      par.out[i] <- list(parallel::mclapply(project$markers[split[[i]]],
                                            function(x) marker.cluster(x,
                                                                       data.type,
                                                                       iteration,
                                                                       model.list,
                                                                       clust.num),
                                            mc.cores = length(cl)))
    }

    if (do.pb)
      setpb(pb, i)
  }

  # Add each result to the project object
  counter = 1
  for (i in 1:length(par.out))
  {
    for (n in 1:length(par.out[[i]]))
    {
      project$markers[[counter]] <- par.out[[i]][[n]]
      counter = counter+1
    }
  }

  # Halt parallel environment
  stopCluster(cl)

  return(project)
}

