#' project.plot
#'
#' Plots the clustering results for a marker object. Includes elipses to indicate distribution variance/rotation.
#'
#' @param project a project object, markers should already be clustered
#' @param indices A numerical index of the markers to plot
#' @param save.loc A string path to the folder to save the images
#' @param recast Logical, should raw data be recast to polar theta and r coordinates?
#'
#' @return Nothing
#'
#' This function relies on the following packages:
#' ggplot2, car, RColorBrewer, stringr

project.plot <- function(project, indices = NULL, save.loc = getwd(), recast = F)
{
  current_dir <- getwd() # Save the current directory
  setwd(save.loc) # Navigate to the save location
  cat('Creating plot folder')
  if (!dir.exists('marker.plots')) {dir.create('marker.plots')} # Make a directory for storing the plots

  cat(paste('\n', 'Plotting markers', sep = ''))

  # Determine how many digits to pad the file names so they sort correctly
  pad.n <- nchar(length(project$markers))

  # Plotting the markers
  for (n in indices) {
    marker <- project$markers[[n]]
    file.name <- paste(str_pad(n, pad.n, pad = "0"), '_', marker$name, '.png', sep = "")

    clustplot <- marker.plot(marker, recast = recast)

    ggsave(paste('marker.plots/', file.name, '.png', sep = ""), plot = clustplot, device = 'png')
  }
}
