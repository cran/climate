% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearest_stations_noaa.R
\name{nearest_stations_nooa}
\alias{nearest_stations_nooa}
\title{List of nearby SYNOP stations for a defined geographical location}
\usage{
nearest_stations_nooa(
  country,
  date = Sys.Date(),
  add_map = TRUE,
  point = NULL,
  no_of_stations = 10,
  ...
)
}
\arguments{
\item{country}{country name; use CAPITAL LETTERS (e.g., "SRI LANKA"), if not used function will found selected
number of nearest stations without country classification}

\item{date}{optionally, a day when measurements were done in all available locations; current Sys.Date used by default}

\item{add_map}{logical - whether to draw a map for a returned data frame (requires maps/mapdata packages)}

\item{point}{a vector of two coordinates (longitude, latitude) for a point we want to find nearest stations to (e.g. c(80, 6))}

\item{no_of_stations}{how many nearest stations will be returned from the given geographical coordinates}

\item{...}{extra arguments to be provided to the \code{\link[graphics:plot.default]{graphics::plot()}} function (only if add_map = TRUE)}
}
\value{
A data.frame with number of nearest station according to given point columns describing stations parameters
(e.g.  ID station, distance from point, geographic coordinates, etc.) where each row represent a measurement,
each station which has a measurements on selected date. If \code{add_map = TRUE} additional map of downloaded data is added.
}
\description{
Returns a data frame of meteorological stations with their coordinates and distance from a given location based on the noaa website.
The returned list is valid only for a given day.
}
\examples{
\donttest{
  nearest_stations_nooa(country = "SRI LANKA", 
  point = c(80, 6),
  add_map = TRUE, 
  no_of_stations = 10)
}

}
