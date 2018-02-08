#' project.import
#'
#' Imports a Genome Studio data export file to only import useful information for clustering (theta, R, raw x, and raw y)
#'
#' @param exportpath the path to the GS data export
#' @param data.type Whether to import raw values or Theta/R values
#'
#' @return A project object containing a vector of marker names and a list of marker objects

# This function relies on the following packages:
# data.table

project.import <- function(exportpath = NULL,
                           raw = T,
                           GSnorm = T)
{
  # Set the data type variable
  if (raw == T)
  {
    data.type = "Raw"
  }  else
  {
    data.type = "GS.transformed"
  }

  # hardcoded transformation method
  transform = "asin"

  # Initialize the project
  project <- list(
    type = data.type,
    marker.names = NULL,
    sample.names = NULL,
    data = list(
      x = NULL,
      y = NULL,
      x.scaled = NULL,
      y.scaled = NULL,
      r = NULL,
      theta  = NULL,
      theta.trans = NULL
    ),
    markers = NULL,
    import.date = NULL
  )



  # Read in the header
  print("Reading header")
  GSheader <-
    as.character(fread(
      exportpath,
      header = F,
      nrows = 1,
      stringsAsFactors = F
    ))

  # Read in marker names
  marker.names <-
    fread(
      exportpath,
      header = T,
      nrows = -1,
      stringsAsFactors = F,
      select = c(2)
    )
  project$marker.names = as.vector(t(marker.names[, 1]))
  rm(marker.names)

  # Process GS normalized data
  if (GSnorm == T) {
    # Read in theta values
    print("Reading theta values")
    theta.ind <- grep('\\.Theta', GSheader)
    r.ind <- grep('\\.R', GSheader)
    theta.table <-
      fread(
        exportpath,
        header = T,
        nrows = -1,
        stringsAsFactors = F,
        select = c(theta.ind)
      )
    theta.table = as.matrix(theta.table)

    # Read in R values
    print("Reading intensity values")
    r.table <-
      fread(
        exportpath,
        header = T,
        nrows = -1,
        stringsAsFactors = F,
        select = c(r.ind)
      )
    r.table = as.matrix(r.table)

    # Removing dimnames from the matrices
    dimnames(theta.table) <- NULL
    dimnames(r.table) <- NULL

    # Adding a slight amount of noise to the data, but making sure to keep theta between 0 and 1
    print("Adding noise to address truncated values")
    theta.table.trans <-
      apply(theta.table[, 1:ncol(theta.table)], c(1, 2), function(x) {
        if (is.na(x)) {
          x
        } else if (x < 0.5) {
          x + abs(rnorm(1, mean = 0, sd = 0.01))
        } else if (x == 0.5) {
          x + rnorm(1, mean = 0, sd = 0.01)
        } else if (x > 0.5) {
          x - abs(rnorm(1, mean = 0, sd = 0.01))
        } else {
          x
        }
      })

    # Transformation functions
    print("Performing transformations")
    asinTransform <- function(p) {
      asin(sqrt(p))
    }
    logitTransform <- function(p) {
      log(p / (1 - p))
    }

    # Perform theta transformation if specified
    if (transform == 'asin') {
      theta.table.trans[, 1:ncol(theta.table.trans)] <-
        asinTransform(theta.table.trans[, 1:ncol(theta.table.trans)])
    } else if (transform == 'logit') {
      theta.table.trans[, 2:ncol(theta.table.trans)] <-
        logitTransform(theta.table.trans[, 2:ncol(theta.table.trans)])
    }

    # Populate the project object
    project$data$r = r.table
    project$data$theta = theta.table
    project$data$theta.trans = theta.table.trans

    # Free up memory
    rm(r.table,
       theta.table,
       theta.table.trans
    )
  }

  # Process raw data
  if (raw == T)
  {
    print("Reading X values")
    # Read in x and y intensity values
    x.index <- grep('\\.X Raw', GSheader)
    y.index <- grep('\\.Y Raw', GSheader)

    x.table <-
      fread(
        exportpath,
        header = T,
        nrows = -1,
        stringsAsFactors = F,
        select = c(x.index)
      )

    x.table <- as.matrix(x.table)

    print("Reading Y values")
    y.table <-
      fread(
        exportpath,
        header = T,
        nrows = -1,
        stringsAsFactors = F,
        select = c(y.index)
      )

    y.table <- as.matrix(y.table)

    # Removing dimnames from the matrices
    dimnames(x.table) <- NULL
    dimnames(y.table) <- NULL

    # Scale the data
    print("Scaling X and Y values")
    scale.out <- scale.raw.data(x.table, y.table)

    # Populate the project object
    project$data$x <- x.table
    project$data$y <- y.table
    project$data$x.scaled <- scale.out[[1]]
    project$data$y.scaled <- scale.out[[2]]

    # Free up memory
    rm(x.table,
       y.table,
       scale.out
    )
  }

  # Extract sample names
  if (GSnorm == T)
  {
    sample.names <- GSheader[theta.ind]
    sample.names <- sub(".Theta", "", sample.names)
  } else if (raw == T)
  {
    sample.names <- GSheader[x.index]
    sample.names <- sub(".X Raw", "", sample.names)
  }

  # Add sample names to the project
  project$sample.names <- sample.names

  # Add import date to the data
  project$import.date <- Sys.time()

  # Create marker objects for each marker
  project$markers <- rep(list(marker.create()), length(project$marker.names))
  for (i in 1:length(project$marker.names))
  {
    name <- project$marker.names[i]

    if (raw == T)
    {
      x <- project$data$x[i,]
      y <- project$data$y[i,]
      x.scaled <- project$data$x.scaled[i,]
      y.scaled <- project$data$y.scaled[i,]
    } else
    {
      x <- NULL
      y <- NULL
      x.scaled <- NULL
      y.scaled <- NULL
    }

    if (GSnorm == T)
    {
      r <- project$data$r[i,]
      theta <- project$data$theta[i,]
      theta.trans <- project$data$theta.trans[i,]
    } else
    {
      r <- NULL
      theta <- NULL
      theta.trans <- NULL
    }

    project$markers[[i]] <- marker.create(name, x, y, x.scaled, y.scaled, r, theta, theta.trans)
  }

  # Update marker object names
  names(project$markers) <- project$marker.names

  return (project)
}
