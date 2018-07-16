#' project.call
#'
#' Calls genotypes for all markers in a project
#'
#' USes the 'Rmixmod' package to fit a Gaussian mixture model using a supplied model list
#'
#' @param project A project object, must be formatted correctly according to the project.create function
#' @param call.type A string variable indicating how to call clusters ('cat.all' = All markers w/ indicator expansion)
#' @param min.prop The minimum proportion of the population for a cluster to be called
#'
#' @return A marker object updated with the best clustering result out of a defined number of iterations
#'
#' This function relies on the following packages:
#' Rmixmod

project.call <-
  function(project,
           call.type = 'cat.all',
           min.prop = 0.01)
  {
    if (call.type == 'cat.all')
    {
      genotype.matrix = matrix(
        NA,
        ncol = length(project$marker.names),
        nrow = length(project$sample.names),
        dimnames = (list(
          project$sample.names, project$marker.names
        ))
      )

      # Add genotype classifications as categorical variables
      for (i in 1:length(project$markers))
      {
        if (!is.null(project$markers[[i]]$model$mixmodout))
          # Check that the marker has been clustered
        {
          cluster.list = project$markers[[i]]$model$mixmodout['partition'] # Extract the partition
          miss.len = length(project$markers[[i]]$data$miss.ind)

          # Reinsert any NA values from the original data (before clustering)
          if (miss.len > 0)
          {
            for (n in 1:miss.len)
            {
              cluster.list <-
                append(cluster.list,
                       NA,
                       after = project$markers[[i]]$data$miss.ind[n] - 1)
            }
          }

          # Remove clusters falling below the proportion threshold
          clust.id <- unique(cluster.list)
          for (x in clust.id)
          {
            prop <- sum(cluster.list == x, na.rm = T) / length(cluster.list)
            if (prop < min.prop)
            {
              cluster.list[which(cluster.list == x)] <- NA
            }
          }
          # Add cluster partition to the genotype frame
          genotype.matrix[, i] <- cluster.list
        }
      }

      # Consolidating the data into a frame
      genotype.frame <- data.frame(row.names(test.matrix), test.matrix)
      col.names(genotype.frame)[1] <- "taxa"

      # Re-cast the matrix with indicator expansion (every cluster is a column coded as 0 and 1)
      genotype.dummy <- dummy.data.frame(genotype.frame, sep = '_', names = project$marker.names)

      # Determine the NA indicator columns and remove them
      NA_cols <- grep("_NA", colnames(genotype.dummy))
      genotype.dummy2 <- genotype.dummy[,-NA_cols]
    }
  }

