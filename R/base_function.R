#' Normalize data using Z-score
#'
#' This function normalizes the data to a Z-score.
#'
#' @param data A numeric vector.
#' @return A numeric vector of Z-scores.
#' @export
zscore <- function(data){
  z_scores <- (data - mean(data, na.rm = TRUE)) / sd(data, na.rm = TRUE)
  return(z_scores)
}

#' Calculate annual mean
#'
#' This function calculates the annual mean from monthly data.
#'
#' @param data A numeric vector of monthly data.
#' @return A numeric vector of annual means.
#' @export
an.mean <- function(data){
  nyear <- length(data) / 12
  data <- matrix(data, nyear, 12, byrow = TRUE)
  an <- apply(data, 1, mean, na.rm = TRUE)
  return(an)
}

#' Calculate annual sum
#'
#' This function calculates the annual sum from monthly data.
#'
#' @param data A numeric vector of monthly data.
#' @return A numeric vector of annual sums.
#' @export
an.sum <- function(data){
  year <- length(data) / 12
  data <- matrix(data, year, 12, byrow = TRUE)
  an <- apply(data, 1, sum, na.rm = TRUE)
  return(an)
}

#' Calculate annual maximum
#'
#' This function calculates the annual maximum from monthly data.
#'
#' @param data A numeric vector of monthly data.
#' @return A numeric vector of annual maximums.
#' @export
an.max <- function(data){
  year <- length(data) / 12
  data <- matrix(data, year, 12, byrow = TRUE)
  an <- apply(data, 1, max, na.rm = TRUE)
  return(an)
}

#' Deseasonalize data by subtracting monthly mean
#'
#' This function removes seasonal effects by subtracting the monthly mean.
#'
#' @param data A numeric vector of monthly data.
#' @return A deseasonalized numeric vector.
#' @export
deseason_add <- function(data){
  nyear <- length(data) / 12
  mat <- matrix(data, nrow = nyear, ncol = 12, byrow = TRUE)
  monthly.mean <- colMeans(mat, na.rm = TRUE)
  monthly.mean <- rep(monthly.mean, nyear)
  results <- data - monthly.mean
  return(results)
}

#' Deseasonalize data by dividing by monthly mean
#'
#' This function removes seasonal effects by dividing by the monthly mean.
#'
#' @param data A numeric vector of monthly data.
#' @return A deseasonalized numeric vector.
#' @export
deseason_mul <- function(data){
  nyear <- length(data) / 12
  mat <- matrix(data, nrow = nyear, ncol = 12, byrow = TRUE)
  monthly.mean <- colMeans(mat, na.rm = TRUE)
  monthly.mean <- rep(monthly.mean, nyear)
  results <- data / monthly.mean
  return(results)
}

#' Convert raster data to array
#'
#' This function converts a raster object to a matrix format.
#'
#' @param raster A raster object.
#' @param lonlat A file path to a CSV with coordinates, or NULL.
#' @param nrows Number of rows in the matrix (default is 360).
#' @param ncols Number of columns in the matrix (default is 720).
#' @return A matrix.
#' @export
as_array <- function(raster, lonlat = NULL, nrows = 360, ncols = 720){
  if (!is.null(lonlat)) {
    lonlat <- read.csv(lonlat)[, 2]
  } else {
    lonlat <- 1:(nrows * ncols)
  }
  matrix <- matrix(NA, nrows * ncols, nlayers(raster))
  for(i in 1:nlayers(raster)){
    matrix[, i] <- as.numeric(as.matrix(subset(raster, i)))
  }
  rownames(matrix) <- lonlat
  colnames(matrix) <- rep(NA, nlayers(raster))
  return(matrix)
}

#' Convert vector to raster
#'
#' This function converts a vector to a raster format.
#'
#' @param x A numeric vector.
#' @param nrows Number of rows in the raster (default is 360).
#' @param ncols Number of columns in the raster (default is 720).
#' @return A raster object.
#' @export
as_raster <- function(x, nrows = 360, ncols = 720){
  mat <- matrix(x, nrows, ncols)
  r <- raster(mat, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(r)
}

#' Convert vector to raster stack
#'
#' This function converts a vector to a raster stack format.
#'
#' @param x A numeric vector.
#' @param nrows Number of rows in the raster (default is 360).
#' @param ncols Number of columns in the raster (default is 720).
#' @return A raster stack object.
#' @export
as_stack <- function(x, nrows = 360, ncols = 720){
  array_data <- array(x, dim = c(nrows, ncols, ncol(x)))
  r_stack <- brick(array_data, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(r_stack)
}

#' Remove outliers from data frame column
#'
#' This function removes outliers from a specified column in a data frame.
#'
#' @param data A data frame.
#' @param column The column name as a string.
#' @return A data frame with outliers removed.
#' @export
remove_outliers_df <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  cleaned_data <- data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ]
  return(cleaned_data)
}

#' Remove outliers from numeric vector
#'
#' This function replaces outliers in a numeric vector with NA.
#'
#' @param data A numeric vector.
#' @return A numeric vector with outliers replaced by NA.
#' @export
remove_outliers <- function(data) {
  if (!is.numeric(data)) {
    stop("Please provide a numeric vector.")
  }
  Q1 <- quantile(data, 0.25, na.rm = TRUE)
  Q3 <- quantile(data, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data[data < lower_bound | data > upper_bound] <- NA
  return(data)
}

#' Remove all objects except specified ones
#'
#' This function removes all objects from the global environment except those specified.
#'
#' @param keep_vars A character vector of variable names to keep.
#' @export
remove_except <- function(keep_vars) {
  all_vars <- ls(envir = .GlobalEnv)
  vars_to_remove <- setdiff(all_vars, keep_vars)
  rm(list = vars_to_remove, envir = .GlobalEnv)
  cat("Kept variables:", keep_vars, "\n")
}
