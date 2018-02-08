#' marker.plot
#'
#' Plots the clustering results for a marker object. Includes elipses to indicate distribution variance/rotation.
#'
#' @param marker A marker object, must have a valid 'mixmodout' cluster model#'
#' @param recast Should
#'
#' @return An object storing all data for a single marker
#'
#' This function relies on the following packages:
# ggplot2, car, RColorBrewer

marker.plot <- function(marker,
                        recast = F)
{
  if (is.null(marker$model$data.type))
  {
    print(paste("Marker", marker$name, "has not been clustered yet!", sep = " "))
    return(NULL)
  }

  data.type = marker$model$data.type

  if (marker$model$data.type == "raw")
  {
    data = data.frame(marker$data$x, marker$data$y)
    names(data) = c("x", "y")
    axes = c("X", "Y")
  }
  else if (marker$model$data.type == "raw.scaled")
  {
    data = data.frame(marker$data$x.scaled, marker$data$y.scaled)
    names(data) = c("x", "y")
    axes = c("X (scaled)", "Y (scaled)")
  }
  else if (marker$model$data.type == "GSnorm")
  {
    data = data.frame(marker$data$theta, marker$data$r)
    names(data) = c("x", "y")
    axes = c("Theta", "Intensity")
  }
  else if (marker$model$data.type == "GS.trans")
  {
    data = data.frame(marker$data$theta.trans, marker$data$r)
    names(data) = c("x", "y")
    axes = c("Theta (Asin transformed)", "Intensity")
  }

  # Remove any NA rows
  if (length(marker$data$miss.ind > 0))
  {
    data = data[-marker$data$miss.ind, ]
  }

  # Add cluster assignments
  data$cluster <- marker$model$mixmodout["partition"]

  # The scatter plot base
  marker.plot <-
    ggplot(data, aes(
      x = x,
      y = y,
      color = factor(cluster)
    ))

  # Calculate max and min values observed
  X.min <- min(data$x, na.rm = T)
  X.max <- max(data$x, na.rm = T)
  Y.min <- min(data$y, na.rm = T)
  Y.max <- max(data$y, na.rm = T)

  # Extract best model
  clust <- marker$model$mixmodout

  # Calculating elipse data
  ellipse.data <- matrix(nrow = 0, ncol = 3)
  center.data <- matrix(nrow = 0, ncol = 3)
  for (i in 1:clust['nbCluster']) {
    ctr <- clust['parameters'][2][i, ]
    A <- clust['parameters'][3][[i]]
    ell.points <-
      data.frame(car::ellipse(
        ctr,
        shape = A,
        radius = 0.98,
        col = "red",
        lty = 2,
        draw = F
      ))
    ell.points$cluster <- i
    center <- data.frame(t(ctr))
    center$cluster <- i
    ellipse.data <- rbind(ellipse.data, ell.points)
    center.data <- rbind(center.data, center)
  }

  # Consolidating ellipse info into a better format
  ellipse.data <- data.frame(ellipse.data)
  center.data <- data.frame(center.data)
  names(ellipse.data) <- c('x', 'y', 'clust')
  names(center.data) <- c('x', 'y', 'clust')

  # Setting up R colors
  myColors <- brewer.pal(8, "Set2")
  names(myColors) <- levels(1:8)
  colScale <- scale_colour_manual(name = "grp", values = myColors)

  # Plotting the marker
  plot <- marker.plot + geom_point(size = 1) +
    geom_path(
      data = ellipse.data,
      aes(x = x, y = y, group = clust),
      color = 'red',
      size = 1,
      inherit.aes = F
    ) +
    geom_point(
      data = center.data,
      aes(x = x, y = y),
      shape = c(center.data[, 3] + 48),
      color = 'black',
      size = 3,
      stroke = 5,
      inherit.aes = F
    ) +
    colScale +
    coord_cartesian(xlim = c((0), X.max + 0.05),
                    ylim = c(0, Y.max + 0.05)) +
    labs(x = axes[1],
         y = axes[2],
         title = marker$name)


  # Transforming the data and converting to polar if data is raw or scaled raw
  if ((data.type == "raw" || data.type == "raw.scaled") && recast == T)
  {
    data2 <- data
    data2$theta <- rad.theta.conv(data2$x, data2$y)
    data2$r <- r.convert(data2$x, data2$y)
    marker.plot2 <-
      ggplot(data2, aes(
        x = theta,
        y = r,
        color = factor(cluster)
      ))

    ellipse.data2 <- ellipse.data
    ellipse.data2$theta <-
      rad.theta.conv(ellipse.data2$x, ellipse.data2$y)
    ellipse.data2$r <- r.convert(ellipse.data2$x, ellipse.data2$y)

    center.data2 <- center.data
    center.data2$theta <-
      rad.theta.conv(center.data2$x, center.data2$y)
    center.data2$r <- r.convert(center.data2$x, center.data2$y)

    # Calculate max and min values observed
    X.min <- min(data2$theta, na.rm = T)
    X.max <- max(data2$theta, na.rm = T)
    Y.min <- min(data2$r, na.rm = T)
    Y.max <- max(data2$r, na.rm = T)

    # Plotting the marker
    plot <- marker.plot2 + geom_point(size = 1) +
      geom_path(
        data = ellipse.data2,
        aes(x = theta, y = r, group = clust),
        color = 'red',
        size = 1,
        inherit.aes = F
      ) +
      geom_point(
        data = center.data2,
        aes(x = theta, y = r),
        shape = c(center.data[, 3] + 48),
        color = 'red',
        size = 3,
        stroke = 3,
        inherit.aes = F
      ) +
      colScale +
      coord_cartesian(xlim = c((0), X.max+0.05), ylim = c(0, Y.max + 0.05)) +
      labs(x = "Theta",
           y = "Intensity",
           title = marker$name)
  }
  return (plot)
}

