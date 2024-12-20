% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_function.R
\name{as_array}
\alias{as_array}
\title{Convert raster data to array}
\usage{
as_array(raster, lonlat = NULL, nrows = 360, ncols = 720)
}
\arguments{
\item{raster}{A raster object.}

\item{lonlat}{A file path to a CSV with coordinates, or NULL.}

\item{nrows}{Number of rows in the matrix (default is 360).}

\item{ncols}{Number of columns in the matrix (default is 720).}
}
\value{
A matrix.
}
\description{
This function converts a raster object to a matrix format.
}
