#' Monthly hydrological data
#'
#' Downloading monthly hydrological data from the danepubliczne.imgw.pl collection
#'
#' @param year vector of years (e.g., 1966:2000)
#' @param coords add coordinates of the stations (logical value TRUE or FALSE)
#' @param station name or ID of hydrological station(s).
#' It accepts names (characters in CAPITAL LETTERS) or stations' IDs (numeric)
#' @param col_names three types of column names possible:
#' "short" - default, values with shorten names,
#' "full" - full English description,
#' "polish" - original names in the dataset
#' @param ... other parameters that may be passed to the 'shortening' function that shortens column names
#' @importFrom XML readHTMLTable
#' @importFrom utils download.file unzip read.csv
#' @importFrom data.table fread
#' @export
#'
#' @examples \donttest{
#'   monthly = hydro_imgw_monthly(year = 2000)
#'   head(monthly)
#' }
#'
hydro_imgw_monthly = function(year, coords = FALSE, station = NULL, col_names= "short", ...) {

  translit = check_locale()

  base_url = "https://danepubliczne.imgw.pl/data/dane_pomiarowo_obserwacyjne/dane_hydrologiczne/"
  interval = "monthly"
  interval_pl = "miesieczne"

  temp = tempfile()
  test_url(link = paste0(base_url, interval_pl, "/"), output = temp)
  a = readLines(temp, warn = FALSE)

  ind = grep(readHTMLTable(a)[[1]]$Name, pattern = "/")
  catalogs = as.character(readHTMLTable(a)[[1]]$Name[ind])
  catalogs = gsub(x = catalogs, pattern = "/", replacement = "")
  catalogs = catalogs[catalogs  %in% as.character(year)]
  if (length(catalogs) == 0) {
    stop("Selected year(s) is not available in the database.", call. = FALSE)
  }
  meta = hydro_metadata_imgw(interval)

  all_data = vector("list", length = length(catalogs))
  for (i in seq_along(catalogs)) {
    catalog = catalogs[i]

    adres = paste0(base_url, interval_pl, "/", catalog, "/mies_", catalog, ".zip")

    temp = tempfile()
    temp2 = tempfile()
    test_url(adres, temp)
    unzip(zipfile = temp, exdir = temp2)
    file1 = paste(temp2, dir(temp2), sep = "/")[1]

    if (translit) {
      data1 = as.data.frame(data.table::fread(cmd = paste("iconv -f CP1250 -t ASCII//TRANSLIT", file1)))
    } else {
      data1 = read.csv(file1, header = FALSE, stringsAsFactors = FALSE, fileEncoding = "CP1250")
    }

    colnames(data1) = meta[[1]][, 1]
    all_data[[i]] = data1
  }
  all_data = do.call(rbind, all_data)

  all_data[all_data == 9999] = NA
  all_data[all_data == 99999.999] = NA
  all_data[all_data == 99.9] = NA
  colnames(all_data) = meta[[1]][, 1]
  # coords
  if (coords) {
    all_data = merge(climate::imgw_hydro_stations, all_data, by.x = "id", by.y = "Kod stacji", all.y = TRUE)
  }
  #station selection
  if (!is.null(station)) {
    if (is.character(station)) {
      all_data = all_data[substr(all_data$`Nazwa stacji`, 1, nchar(station)) == station, ]
      if (nrow(all_data) == 0) {
        stop("Selected station(s) is not available in the database.", call. = FALSE)
      }
    } else if (is.numeric(station)) {
      all_data = all_data[all_data$`Kod stacji` %in% station, ]
      if (nrow(all_data) == 0) {
        stop("Selected station(s) is not available in the database.", call. = FALSE)
      }
    } else {
      stop("Selected station(s) are not in the proper format.", call. = FALSE)
    }
  }

  all_data = all_data[order(all_data$`Nazwa stacji`,
                            all_data$`Rok hydrologiczny`,
                            all_data$`Wskaznik miesiaca w roku hydrologicznym`), ]
  all_data = hydro_shortening_imgw(all_data, col_names = col_names, ...)

  return(all_data)
}
