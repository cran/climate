context("meteo_imgw")

test_that("sounding_wyoming works!", {
  profile <- sounding_wyoming(wmo_id = 12120,
                              yy = 2019,
                              mm = 4,
                              dd = 4,
                              hh = 0)

  # expected error
  testthat::expect_message(sounding_wyoming(wmo_id = 12220,
                                            yy = 2019,
                                            mm = 4,
                                            dd = 4,
                                            hh = 0))

  # expected error
  testthat::expect_error(sounding_wyoming(wmo_id = c(12220, 12375),
                                          yy = 2019,
                                          mm = 4,
                                          dd = 4,
                                          hh = 0))
  
  # expected error for bufr
  testthat::expect_error(sounding_wyoming(wmo_id = 12375,
                                          yy = 2019,
                                          mm = 4,
                                          dd = 4,
                                          hh = 0,
                                          bufr = TRUE))
})
