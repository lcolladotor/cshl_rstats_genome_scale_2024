## ----"initial_weekday_praise"----------------------------------------------
weekday_praise <- function(date = Sys.Date()) {
    date <- as.Date(date)
    date_weekday <- weekdays(date)
    paste0(date_weekday, ": ", praise::praise())
}
weekday_praise()
weekday_praise("2024-06-09")


## ----"weekday_praise_full_function"----------------------------------------
#' Praise a weekday
#'
#' Given a date, figure out which weekday it was, then write a positive
#' message.
#'
#' @param date A `base::Date` object or a `character()` in a format that can be
#' converted to a `base::Date` object with `base::as.Date()`.
#'
#' @importFrom praise praise
#' @export
#' @examples
#'
#' ## Praise the current weekday
#' weekday_praise()
#'
#' ## Praise the date we started teaching
#' weekday_praise("2024-06-09")
#'
#' ## Praise the current weekday in a reproducible way
#' set.seed(20240610)
#' weekday_praise()
#'
#' ## Verify that it's reproducible
#' set.seed(20240610)
#' weekday_praise()
weekday_praise <- function(date = Sys.Date()) {
    date <- as.Date(date)
    date_weekday <- weekdays(date)
    paste0(date_weekday, ": ", praise::praise())
}


## ----"weekday_praise_tests"------------------------------------------------
library("testthat")

## Verify that we get the result we wanted
set.seed(20240610)
expect_equal(weekday_praise("2024-06-09"), "Sunday: You are wondrous!")

## Verify that we get an error if the input is not correct
expect_error(weekday_praise("240609"))

## Should work for a vector input
expect_equal(length(weekday_praise(c("2024-06-09", "2024-06-10"))), 2L)

