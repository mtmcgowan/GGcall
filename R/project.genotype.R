#' project.genotype
#'
#' Uses Rmixmod clusters from a GGcall project to call genotypes
#'
#' @param project A project object, must contain cluster information from the project.cluster function
#' @param format If clusters were called using raw intensities or scaled raw intensities, should they be transformed to theta an R? (can be 'none' or 'transformed')
#'
#' @return A list of genotype information: 3 matrices corresponding to the cluster, the x/theta cluster center, the y/R cluster center, the specified format, and the data type used for clustering

project.genotype <- function(project, clusters = F, format = "none")
{
    # Get the type of model (either raw, raw.scaled, GS.norm, or GS.trans)
    model_type = project$markers[[1]]$model$data.type

    # Initialize the genotype matrix
    clust_gen_mat <-  matrix(nrow = length(project$sample.names),
                             ncol = length(project$marker.names),
                             dimnames = list(project$sample.names, project$marker.names)
    )

    x_gen_mat <- matrix(nrow = length(project$sample.names),
                        ncol = length(project$marker.names),
                        dimnames = list(project$sample.names, project$marker.names)
    )

    y_gen_mat <- matrix(nrow = length(project$sample.names),
                        ncol = length(project$marker.names),
                        dimnames = list(project$sample.names, project$marker.names)
    )

    # Iterate through each marker assigning either clusters or cluster centers as the genotype
    for (i in project$markers)
    {
        clusters <- i$model[['mixmodout']]['partition']
        x_means <- i$model$mixmodout['parameters']['mean'][,1]
        y_means <- i$model$mixmodout['parameters']['mean'][,2]

        x_genotype <- rep(NA, length(clusters))
        y_genotype <- rep(NA, length(clusters))

        for (n in 1:i$model$mixmodout['nbCluster'])
        {
            x_genotype[clusters == n] <- x_means[n]
            y_genotype[clusters == n] <- y_means[n]
        }

        # Check for missing values and update genotype matrix
        if (i$model$data.type == "raw.scaled")
        {
            x_na <- which(is.na(i$data$x.scaled))
            y_na <- which(is.na(i$data$y.scaled))

            all_na <- unique(c(x_na, y_na))
            for (na_index in all_na)
            {
                x_genotype <- append(x_genotype, NA, na_index-1)
                y_genotype <- append(y_genotype, NA, na_index-1)
                clusters <- append(clusters, NA, na_index-1)
            }
        } else if (i$model$data.type == "raw")
        {
            x_na <- which(is.na(i$data$x))
            y_na <- which(is.na(i$data$y))

            all_na <- unique(c(x_na, y_na))
            for (na_index in all_na)
            {
                x_genotype <- append(x_genotype, NA, na_index-1)
                y_genotype <- append(y_genotype, NA, na_index-1)
                clusters <- append(clusters, NA, na_index-1)
            }

        } else if (i$model$data.type == "GS.norm")
        {
            x_na <- which(is.na(i$data$theta))
            y_na <- which(is.na(i$data$r))

            all_na <- unique(c(x_na, y_na))
            for (na_index in all_na)
            {
                x_genotype <- append(x_genotype, NA, na_index-1)
                y_genotype <- append(y_genotype, NA, na_index-1)
                clusters <- append(clusters, NA, na_index-1)
            }
        } else if (i$model$data.type == "GS.trans")
        {
            x_na <- which(is.na(i$data$theta.trans))
            y_na <- which(is.na(i$data$r))

            all_na <- unique(c(x_na, y_na))
            for (na_index in all_na)
            {
                x_genotype <- append(x_genotype, NA, na_index-1)
                y_genotype <- append(y_genotype, NA, na_index-1)
                clusters <- append(clusters, NA, na_index-1)
            }
        }

        # Check and transpose to theta and R if necessary
        if (format == "transformed" && (model_type != "GS.norm" && model_type != "GS.trans"))
        {
            theta_genotype <- rad.theta.conv(x_genotype, y_genotype)
            r_genotype <- r.convert(x_genotype, y_genotype)

            x_genotype = theta_genotype
            y_genotype = r_genotype
        }

        # Update the genotype matrix
        clust_gen_mat[,i$name] = clusters
        x_gen_mat[,i$name] = x_genotype
        y_gen_mat[,i$name] = y_genotype
    }

    genotypes = list(clust_gen_mat, x_gen_mat, y_gen_mat, model_type, format)
    names(genotypes) = c("clusters", "x_genotype", "y_genotype", "model_type", "format")
    return(genotypes)
}
