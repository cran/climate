% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearest_stations_imgw.R
\name{nearest_stations_imgw}
\alias{nearest_stations_imgw}
\title{List of nearby meteorological or hydrological IMGW-PIB stations in Poland}
\usage{
nearest_stations_imgw(
  type = "meteo",
  rank = "synop",
  year = 2018,
  add_map = TRUE,
  point = NULL,
  no_of_stations = 50,
  ...
)
}
\arguments{
\item{type}{data name; "meteo" (default), "hydro"}

\item{rank}{rank of the stations: "synop" (default), "climate", or "precip"; Only valid if type = "meteo"}

\item{year}{select year for searching nearest station}

\item{add_map}{logical - whether to draw a map for a returned data frame (requires maps/mapdata packages)}

\item{point}{a vector of two coordinates (longitude, latitude) for a point we want to find nearest stations to (e.g. c(15, 53)); If not provided calculated as a mean longitude and latitude for the entire dataset}

\item{no_of_stations}{how many nearest stations will be returned from the given geographical coordinates. 50 used by default}

\item{...}{extra arguments to be provided to the \code{\link[graphics:plot.default]{graphics::plot()}} function (only if add_map = TRUE)}
}
\value{
A data.frame with a list of nearest stations. Each row represents metadata for station which collected measurements in a given year. Particular columns contain stations metadata (e.g. station ID, geographical coordinates, official name, distance in kilometers from a given coordinates).
}
\description{
Returns a data frame of meteorological or hydrological stations with their coordinates in particular year.
The returned object is valid only for a given year and type of stations (e.g. "synop", "climate" or "precip"). If \code{add_map = TRUE} additional map of downloaded data is added.
}
\examples{
\donttest{
  df =  nearest_stations_imgw(type = "meteo", 
  rank = "synop",
  year = 2018,
  point = c(17, 52),
  add_map = TRUE, 
  no_of_stations = 4)
}

}
