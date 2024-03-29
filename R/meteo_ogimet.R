#' Scrapping meteorological (Synop) data from the Ogimet webpage
#'
#' Downloading hourly or daily (meteorological) data from the Synop stations available at https://www.ogimet.com/
#'
#' @param interval 'daily' or 'hourly' dataset to retrieve - given as character
#' @param date start and finish date (e.g., date = c("2018-05-01", "2018-07-01")) - character or Date class object
#' @param coords add geographical coordinates of the station (logical value TRUE or FALSE)
#' @param station WMO ID of meteorological station(s). Character or numeric vector
#' @param precip_split whether to split precipitation fields into 6/12/24h
#' numeric fields (logical value TRUE (default) or FALSE); valid only for hourly time step
#' @importFrom XML readHTMLTable
#' 
#' @export
#' @return A data.frame of measured values with columns describing the meteorological parameters (e.g. air temperature, wind speed, cloudines). 
#' Depending on the interval, at a given hour or day. Different parameters are returned for daily and hourly datasets.
#' \enumerate{
#'  \item station_ID - WMO station identifier
#'  \item Lon - longitude
#'  \item Lat - latitude
#'  \item Date - date (and time) of observations
#'  \item TC - air temperature at 2 metres  above ground level. Values given in Celsius degrees
#'  \item TdC - dew point temperature at 2 metres  above ground level. Values given in Celsius degrees
#'  \item TmaxC - maximum air temperature at 2 metres  above ground level. Values given in Celsius degrees
#'  \item TminC - minimum air temperature at 2 metres  above ground level. Values given in Celsius degrees
#'  \item ddd - wind direction
#'  \item ffkmh - wind speed in km/h
#'  \item Gustkmh - wind gust in km/h
#'  \item P0hpa - air pressure at elevation of the station in hPa
#'  \item PseahPa - sea level pressure in hPa
#'  \item PTnd - pressure tendency in hPa
#'  \item Nt - total cloud cover
#'  \item Nh - cloud cover by high-level cloud fraction
#'  \item HKm - height of cloud base
#'  \item InsoD1 - insolation in hours
#'  \item Viskm - visibility in kilometres
#'  \item Snowcm - depth of snow cover in centimetres
#'  \item pr6 - precicipitation totals in 6 hours
#'  \item pr12 - precicipitation totals in 12 hours
#'  \item pr24 - precicipitation totals in 24 hours
#'  \item TemperatureCAvg - average air temperature at 2 metres above ground level. Values given in Celsius degrees
#'  \item TemperatureCMax - maximum air temperature at 2 metres above ground level. Values given in Celsius degrees
#'  \item TemperatureCMin - minimum air temperature at 2 metres above ground level. Values given in Celsius degrees
#'  \item TdAvgC - average dew point temperature at 2 metres above ground level. Values given in Celsius degrees
#'  \item HrAvg - average relative humidity. Values given in %
#'  \item WindkmhDir - wind direction
#'  \item WindkmhInt - wind speed in km/h
#'  \item WindkmhGust - wind gust in km/h
#'  \item PresslevHp - Sea level pressure in hPa
#'  \item Precmm - precipitation totals in mm
#'  \item TotClOct - total cloudiness in octants
#'  \item lowClOct - cloudiness by low level clouds in octants
#'  \item SunD1h - sunshine duration in hours
#'  \item PreselevHp - atmospheric pressure measured at altitude of station in hPa
#'  \item SnowDepcm - depth of snow cover in centimetres
#'  }


#' @examples 
#' \donttest{
#'   # downloading data for Poznan-Lawica
#'   # poznan = meteo_ogimet(interval = "daily", 
#'   #                      date = c(Sys.Date()-30, Sys.Date()),
#'   #                      station = 12330, 
#'   #                      coords = TRUE)
#'   # head(poznan)
#' }
#'
meteo_ogimet = function(interval, date, coords = FALSE, station, precip_split = TRUE) {
  if (interval == "daily") {
    # daily
    if (!precip_split) {
      warning("The `precip_split` argument is only valid for hourly time step", call. = FALSE)
    }
    all_data = ogimet_daily(date = date,  coords = coords, station = station)
  } else if (interval == "hourly") {
    #hourly
    all_data = ogimet_hourly(date = date,  coords = coords, station = station,
                              precip_split = precip_split)
  } else{
    stop("Wrong `interval` value. It should be either 'hourly' or 'daily'")
  }
  return(all_data)
}
