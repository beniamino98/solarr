#' Bologna model specification
#'
#' A saved `solarModel_spec` object for Bologna. The data file
#' `Bologna.RData` loads an object named `spec`.
#'
#' @format An R6 object of class `solarModel_spec`.
#' @details
#' The object is included as an example package specification. It has
#' `place` set to `"Bologna"` and `target` set to `"GHI"`.
#'
#' @usage data(Bologna)
#' @docType data
#' @keywords datasets
#' @name Bologna
NULL

#' CAMS data by location
#'
#' A list of daily solar radiation data sets for selected locations.
#'
#' @format A named list. Each element is a tibble with 7,465 rows and 8
#' columns:
#' \describe{
#'   \item{date}{Date.}
#'   \item{n}{Numeric day index.}
#'   \item{Year}{Numeric year.}
#'   \item{Month}{Numeric month.}
#'   \item{Day}{Integer day of month.}
#'   \item{GHI}{Numeric global horizontal irradiance value.}
#'   \item{clearsky}{Numeric clear-sky radiation value.}
#'   \item{H0}{Numeric extraterrestrial radiation value.}
#' }
#'
#' @details
#' The list contains named elements for the locations available to
#' `solarModel_spec$specification()`. Individual elements include `place` and
#' `coords` attributes.
#'
#' @usage data(CAMS_data)
#' @docType data
#' @keywords datasets
#' @name CAMS_data
NULL

#' Solar spectral data
#'
#' A tabular data set of spectral irradiance values.
#'
#' @format A tibble with 2,002 rows and 4 numeric columns:
#' \describe{
#'   \item{Lambda}{Wavelength values.}
#'   \item{Am0}{Numeric spectral values.}
#'   \item{Am1}{Numeric spectral values.}
#'   \item{Am2}{Numeric spectral values.}
#' }
#'
#' @usage data(solarSpectra)
#' @docType data
#' @keywords datasets
#' @name solarSpectra
NULL
