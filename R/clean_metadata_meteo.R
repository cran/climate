#' Meteorological metadata cleaning
#'
#' Internal function for meteorological metadata cleaning
#' @param address URL address of the metadata file
#' @param rank stations' rank
#' @param interval temporal interval
#' @importFrom utils read.fwf
#' @importFrom stats na.omit
#' @keywords internal
#'
#' @examples
#' \donttest{
#'   my_add = paste0("https://danepubliczne.imgw.pl/data/dane_pomiarowo_obserwacyjne/",
#'                   "dane_meteorologiczne/dobowe/synop/s_d_format.txt")
#'   climate:::clean_metadata_meteo(address = my_add, rank = "synop", interval = "hourly")
#' }
#'

clean_metadata_meteo = function(address, rank = "synop", interval = "hourly") {

  temp = tempfile()
  test_url(link = address, output = temp)
  a = readLines(temp, warn = FALSE)

  a = iconv(a, from = "cp1250", to = "ASCII//TRANSLIT") # usuwamy polskie znaki, bo to robi spore "kuku"
  a = gsub(a, pattern = "\\?", replacement = "") # usuwamy znaki zapytania powstale po konwersji

  # additional workarounds for mac os but not only...
  a = gsub(x = a, pattern = "'", replacement = "")
  a = gsub(x = a, pattern = "\\^0", replacement = "")
  a = data.frame(V1 = a[nchar(a) > 0], stringsAsFactors = FALSE)
  # this one does not work on windows
  # a = suppressWarnings(na.omit(read.fwf(address, widths = c(1000),
  #                                        fileEncoding = "CP1250", stringsAsFactors = FALSE)))
  length_char = max(nchar(a$V1), na.rm = TRUE)

  if (rank == "precip" && interval == "hourly") length_char = 40 # wyjatek dla precipow
  if (rank == "precip" && interval == "daily") length_char = 40 # wyjatek dla precipow dobowych
  if (rank == "synop" && interval == "hourly") length_char = 60 # wyjatek dla synopow terminowych
  if (rank == "climate" && interval == "monthly") length_char = 52 # wyjatek dla synopow terminowych

  field = substr(a$V1, length_char - 3, length_char)

  if (rank == "synop" && interval == "monthly") {
    length_char = as.numeric(names(sort(table(nchar(a$V1)), decreasing = TRUE)[1])) + 2
    field = substr(a$V1, length_char - 3, length_char + 2)
  }

  a$field1 = suppressWarnings(as.numeric(unlist(lapply(strsplit(field, "/"), function(x) x[1]))))
  a$field2 = suppressWarnings(as.numeric(unlist(lapply(strsplit(field, "/"), function(x) x[2]))))

  a$V1 = trimws(substr(a$V1, 1, nchar(a$V1) - 3))

  strsplit(x = a$V1, split = "/")
  #a = a[nchar(a$V1)>2,] # remove empty or almost empty rows
  a = a[!(is.na(a$field1) & is.na(a$field2)), ] # remove info about status
  colnames(a)[1] = "parameters"
  return(a)
}
