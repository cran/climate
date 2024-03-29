% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteo_imgw_daily.R
\name{meteo_imgw_daily}
\alias{meteo_imgw_daily}
\title{Daily IMGW meteorological data}
\usage{
meteo_imgw_daily(
  rank = "synop",
  year,
  status = FALSE,
  coords = FALSE,
  station = NULL,
  col_names = "short",
  ...
)
}
\arguments{
\item{rank}{rank of the stations: "synop" (default), "climate", or "precip"}

\item{year}{vector of years (e.g., 1966:2000)}

\item{status}{leave the columns with measurement and observation statuses
(default status = FALSE - i.e. the status columns are deleted)}

\item{coords}{add coordinates of the station (logical value TRUE or FALSE)}

\item{station}{name of meteorological station(s).
It accepts names (characters in CAPITAL LETTERS); stations' IDs (numeric) are
no longer valid}

\item{col_names}{three types of column names possible:
"short" - default, values with shorten names,
"full" - full English description,
"polish" - original names in the dataset}

\item{...}{other parameters that may be passed to the 'shortening' function that
shortens column names}
}
\description{
Downloading daily (meteorological) data from the SYNOP / CLIMATE / PRECIP stations
available in the danepubliczne.imgw.pl collection
}
\examples{
\donttest{
  daily = meteo_imgw_daily(rank = "climate", year = 2000)
  head(daily)
}

}
