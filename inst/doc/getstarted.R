## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
old <- options(scipen = 999)

## ----stations , eval=T, fig.width=7,fig.height=7, fig.fullwidth=TRUE----------
library(climate)
ns = nearest_stations_ogimet(country = c("United Kingdom", "France"),
                             point = c(-3, 50),
                             no_of_stations = 50, 
                             add_map = TRUE)

## ----stations-2, eval=T-------------------------------------------------------
if (is.data.frame(ns)) {
 knitr::kable(head(ns, 15))
}

## ----stations-3, eval=T, fig.width=7, fig.height=7, fig.fullwidth=T-----------
library(climate)
PL = stations_ogimet(country = "Poland", add_map = TRUE)

if (is.data.frame(PL)) {
    knitr::kable(head(PL))
}

## ----noaa_svalbard, include=FALSE---------------------------------------------
df = readRDS(system.file("extdata/vignettes/svalbard_noaa.rds", package = "climate"))

## ----windrose,eval=F----------------------------------------------------------
# # downloading data with NOAA service:
# df = meteo_noaa_hourly(station = "010080-99999", year = 2016)
# 
# # You can also download the same (but more granular) data with Ogimet.com (example for year 2016):
# # df = meteo_ogimet(interval = "hourly",
# #                   date = c("2016-01-01", "2016-12-31"),
# #                   station = c("01008"))

## ----noaa-kable,eval=T--------------------------------------------------------
knitr::kable(head(df))

## ----sonda-read, eval=T, include=F, echo=F------------------------------------
library(climate)
data("profile_demo")
df2 = profile_demo[[1]] 
colnames(df2)[c(1, 3:4)] = c("PRESS", "TEMP", "DEWPT") # changing column names

## ----sonda, eval=F, include=T-------------------------------------------------
# profile_demo <- sounding_wyoming(wmo_id = 12120,
#                                  yy = 2000,
#                                  mm = 3,
#                                  dd = 23,
#                                  hh = 0)
# df2 = profile_demo[[1]]
# colnames(df2)[c(1, 3:4)] = c("PRESS", "TEMP", "DEWPT") # changing column names

## ----sonda2, echo=FALSE-------------------------------------------------------
knitr::kable(head(df2, 10), caption = "Exemplary data frame of sounding preprocessing")

## ----imgw_meteo, include=FALSE------------------------------------------------
df = readRDS(system.file("extdata/vignettes/leba_monthly.rds", package = "climate"))

## ----imgw_meteo-2, eval=FALSE, include=TRUE-----------------------------------
# library(climate)
# df = meteo_imgw(interval = "monthly", rank = "synop", year = 1991:2000, station = "ŁEBA")
# # please note that sometimes 2 names are used for the same station in different years

## ----imgw_meteo-3, fig.width=7, fig.height=7, fig.fullwidth=TRUE, error=TRUE, eval=TRUE, include=TRUE----
try({
suppressMessages(library(dplyr))
df2 = dplyr::select(df, station:t2m_mean_mon, rr_monthly)

monthly_summary = df2 %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarise(tmax = mean(tmax_abs, na.rm = TRUE), 
                   tmin = mean(tmin_abs, na.rm = TRUE),
                   tavg = mean(t2m_mean_mon, na.rm = TRUE), 
                   precip = sum(rr_monthly) / n_distinct(yy))            

monthly_summary = as.data.frame(t(monthly_summary[, c(5, 2, 3, 4)])) 
monthly_summary = round(monthly_summary, 1)
colnames(monthly_summary) = month.abb
})

## ----imgw_meteo2, echo=FALSE, error=TRUE--------------------------------------
try({
knitr::kable(head(monthly_summary), 
             caption = "Exemplary data frame of meteorological preprocessing.")
})

## ----data, eval=TRUE, include=FALSE, echo=FALSE-------------------------------
h = readRDS(system.file("extdata/vignettes/hydro_monthly.rds", package = "climate"))

## ----data-2, eval=FALSE, include=TRUE-----------------------------------------
# library(climate)
# library(dplyr)
# library(tidyr)
# h = hydro_imgw(interval = "monthly", year = 2001:2002, coords = TRUE)

## ----data-3, eval=TRUE, include=TRUE, echo=TRUE-------------------------------
knitr::kable(head(h))

## ----filtering, eval=TRUE, include=TRUE---------------------------------------
h2 = h %>%
  dplyr::filter(idex == 3) %>%
  dplyr::select(id, station, X, Y, hyy, Q) %>%
  dplyr::group_by(hyy, id, station, X, Y) %>%
  dplyr::summarise(annual_mean_Q = round(mean(Q, na.rm = TRUE), 1)) %>% 
  tidyr::pivot_wider(names_from = hyy, values_from = annual_mean_Q)

knitr::kable(head(h2))

## ----filtering2, echo=FALSE, eval=FALSE---------------------------------------
# 
# knitr::kable(head(h2),
#              caption = "Exemplary data frame of hydrological preprocesssing.")

## ----setup_restore, include = FALSE-------------------------------------------
options(old)

