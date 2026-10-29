test_that("seasonalSolarFunctions initializes with the documented defaults", {
  solar <- seasonalSolarFunctions$new()

  expect_s3_class(solar, "seasonalSolarFunctions")
  expect_identical(solar$update_method(), "spencer")
  expect_true(solar$legal_hour)
  expect_equal(solar$Gsc, 1367)
})

test_that("seasonalSolarFunctions accepts and updates supported methods", {
  solar <- seasonalSolarFunctions$new(method = "cooper", legal_hour = FALSE)

  expect_identical(solar$update_method(), "cooper")
  expect_false(solar$legal_hour)
  expect_identical(solar$update_method("spencer"), "spencer")
  expect_identical(solar$update_method(), "spencer")
})

test_that("seasonalSolarFunctions rejects unsupported methods", {
  expect_error(seasonalSolarFunctions$new(method = "unknown"))

  solar <- seasonalSolarFunctions$new()

  expect_error(solar$update_method("unknown"))
})

test_that("B spans one full cycle over 365 days", {
  solar <- seasonalSolarFunctions$new()

  result <- solar$B(c(0, 365))

  expect_equal(result, c(0, 2 * pi), tolerance = 1e-12)
})

test_that("angle conversions are inverse and vectorised", {
  solar <- seasonalSolarFunctions$new()
  degrees <- c(-180, -90, 0, 45, 180)

  radians <- solar$radiant(degrees)

  expect_length(radians, length(degrees))
  expect_equal(solar$degree(radians), degrees, tolerance = 1e-12)
})

test_that("angle conversions preserve empty and missing values", {
  solar <- seasonalSolarFunctions$new()

  expect_length(solar$radiant(numeric()), 0)
  expect_length(solar$degree(numeric()), 0)
  expect_true(is.na(solar$radiant(NA_real_)))
  expect_true(is.na(solar$degree(NA_real_)))
})

test_that("equation of time returns the accepted January reference value", {
  solar <- seasonalSolarFunctions$new()

  result <- as.numeric(solar$E(as.Date("2022-01-01")), units = "mins")

  expect_equal(result, -2.90, tolerance = 0.02)
})

test_that("elevation is zero at sea level and increases with altitude", {
  solar <- seasonalSolarFunctions$new()

  result <- solar$elevation(c(0, 1000, 2000))

  expect_equal(result[[1]], 0)
  expect_true(all(diff(result) > 0))
  expect_true(all(result >= 0))
})

test_that("solar_time preserves instants and applies the longitude correction", {
  solar <- seasonalSolarFunctions$new()
  clock_time <- as.POSIXct("2022-01-01 12:00:00", tz = "UTC")

  at_meridian <- solar$solar_time(clock_time, lon = 15, lon_st = 15, tz = "UTC")
  west_of_meridian <- solar$solar_time(clock_time, lon = 14, lon_st = 15, tz = "UTC")

  expect_s3_class(at_meridian, "POSIXt")
  expect_equal(
    as.numeric(difftime(west_of_meridian, at_meridian, units = "mins")),
    -4,
    tolerance = 1e-10
  )
})

test_that("solar_hour includes minutes and seconds", {
  solar <- seasonalSolarFunctions$new()
  solar_time <- as.POSIXct("2022-01-01 12:30:30", tz = "UTC")

  result <- solar$solar_hour(solar_time)

  expect_equal(result, 12 + 30 / 60 + 30 / 3600, tolerance = 1e-12)
})

test_that("hour_angle is zero at solar noon and changes 15 degrees per hour", {
  solar <- seasonalSolarFunctions$new()
  solar_time <- as.POSIXct(
    c("2022-01-01 11:00:00", "2022-01-01 12:00:00", "2022-01-01 13:00:00"),
    tz = "UTC"
  )

  result <- solar$hour_angle(solar_time)

  expect_equal(result, c(-15, 0, 15), tolerance = 1e-12)
})

test_that("incidence_angle returns bounded cosine values", {
  solar <- seasonalSolarFunctions$new()
  solar_time <- as.POSIXct(
    c("2022-06-21 00:00:00", "2022-06-21 12:00:00"),
    tz = "UTC"
  )

  result <- solar$incidence_angle(solar_time, lat = 43, alt = 0)

  expect_length(result, 2)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(result[[1]], 0)
  expect_gt(result[[2]], 0.9)
})

test_that("azimut_angle returns finite angles after the undefined-symbol regression", {
  solar <- seasonalSolarFunctions$new()
  solar_time <- as.POSIXct(
    c("2022-06-21 10:00:00", "2022-06-21 14:00:00"),
    tz = "UTC"
  )

  result <- solar$azimut_angle(solar_time, lat = 43, alt = 0)

  expect_length(result, 2)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= -180 & result <= 180))
  expect_lt(result[[1]], 0)
  expect_gt(result[[2]], 0)
})

test_that("Gon stays close to the solar constant for both supported methods", {
  dates <- as.Date(c("2022-01-03", "2022-07-04"))

  spencer <- seasonalSolarFunctions$new(method = "spencer")$Gon(dates)
  cooper <- seasonalSolarFunctions$new(method = "cooper")$Gon(dates)

  expect_length(spencer, 2)
  expect_length(cooper, 2)
  expect_true(all(spencer > 1300 & spencer < 1420))
  expect_true(all(cooper > 1300 & cooper < 1420))
  expect_gt(spencer[[1]], spencer[[2]])
  expect_gt(cooper[[1]], cooper[[2]])
})

test_that("Gon derivative agrees with a centred finite difference", {
  solar <- seasonalSolarFunctions$new(method = "spencer")
  day <- 172
  step <- 1e-3

  numerical <- (solar$Gon(day + step) - solar$Gon(day - step)) / (2 * step)
  result <- solar$Gon(day, deriv = TRUE)

  expect_equal(result, numerical, tolerance = 1e-6)
})

test_that("declination follows equinox and solstice reference values", {
  solar <- seasonalSolarFunctions$new(method = "spencer")
  dates <- as.Date(c("2022-03-20", "2022-06-21", "2022-12-21"))

  result <- solar$declination(dates)

  expect_equal(result[[1]], 0, tolerance = 0.6)
  expect_equal(result[[2]], 23.44, tolerance = 0.1)
  expect_equal(result[[3]], -23.44, tolerance = 0.1)
})

test_that("declination derivative agrees with a centred finite difference", {
  solar <- seasonalSolarFunctions$new(method = "cooper")
  day <- 172
  step <- 1e-3

  numerical <-
    (solar$declination(day + step) - solar$declination(day - step)) / (2 * step)
  result <- solar$declination(day, deriv = TRUE)

  expect_equal(result, numerical, tolerance = 1e-6)
})

test_that("sunset angle and daylight duration are physically consistent", {
  solar <- seasonalSolarFunctions$new()
  dates <- as.Date(c("2022-06-21", "2022-12-21"))

  sunset <- solar$sunset_hour_angle(dates, lat = 43, alt = 0)
  daylight <- as.numeric(solar$sun_hours(dates, lat = 43, alt = 0), units = "hours")

  expect_true(all(sunset >= 0 & sunset <= 180))
  expect_equal(daylight, 2 * sunset / 15, tolerance = 1e-10)
  expect_gt(daylight[[1]], daylight[[2]])
})

test_that("sunset angle derivative agrees with a centred finite difference", {
  solar <- seasonalSolarFunctions$new()
  day <- 172
  step <- 1e-3

  numerical <- (
    solar$sunset_hour_angle(day + step, lat = 43, alt = 0) -
      solar$sunset_hour_angle(day - step, lat = 43, alt = 0)
  ) / (2 * step)
  result <- solar$sunset_hour_angle(day, lat = 43, alt = 0, deriv = TRUE)

  expect_equal(result, numerical, tolerance = 1e-6)
})

test_that("Hon is positive and larger in summer at mid-latitudes", {
  solar <- seasonalSolarFunctions$new()
  dates <- as.Date(c("2022-06-21", "2022-12-21"))

  result <- solar$Hon(dates, lat = 43, alt = 0)

  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))
  expect_gt(result[[1]], result[[2]])
})

test_that("Hon derivative agrees with a centred finite difference", {
  solar <- seasonalSolarFunctions$new()
  day <- 172
  step <- 1e-3

  numerical <- (
    solar$Hon(day + step, lat = 43, alt = 0) -
      solar$Hon(day - step, lat = 43, alt = 0)
  ) / (2 * step)
  result <- solar$Hon(day, lat = 43, alt = 0, deriv = TRUE)

  expect_equal(result, numerical, tolerance = 1e-6)
})

test_that("solar_altitude matches seasonal noon reference values", {
  solar <- seasonalSolarFunctions$new()
  dates <- as.Date(c("2022-06-21", "2022-12-21"))

  result <- solar$solar_altitude(dates, lat = 43)

  expect_equal(result[[1]], 70.44, tolerance = 0.1)
  expect_equal(result[[2]], 23.56, tolerance = 0.1)
})

test_that("solar_angles returns a complete deterministic daily grid", {
  solar <- seasonalSolarFunctions$new()
  instant <- as.POSIXct("2022-01-15 12:00:00", tz = "UTC")

  result <- solar$solar_angles(
    instant,
    lat = 43,
    lon = 12,
    alt = 0,
    lon_st = 15,
    by = "1 hour",
    tz = "UTC"
  )

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c(
      "date", "clocktime", "solartime", "lat", "lon", "omega",
      "declination", "omega_min", "omega_max", "sun_hours", "theta", "E"
    )
  )
  expect_equal(nrow(result), 25)
  expect_s3_class(result$clocktime, "POSIXct")
  expect_s3_class(result$solartime, "POSIXt")
  expect_false(anyNA(result$solartime))
  expect_equal(result$omega_min, -result$omega_max)
  expect_true(all(result$theta >= 0 & result$theta <= 1))
})

test_that("clearsky returns bounded irradiance for supported climates", {
  solar <- seasonalSolarFunctions$new()
  climates <- c("No Correction", "Summer", "Winter", "Subartic Summer", "Tropical")

  result <- vapply(
    climates,
    function(climate) solar$clearsky(cosZ = 0.8, G0 = 1000, alt = 0, clime = climate),
    numeric(1)
  )

  expect_named(result, climates)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1000))
})

test_that("clearsky rejects unsupported climate corrections", {
  solar <- seasonalSolarFunctions$new()

  expect_error(
    solar$clearsky(cosZ = 0.8, G0 = 1000, alt = 0, clime = "Monsoon")
  )
})
