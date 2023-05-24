##' We need to map "dates" onto [`dust::dust`]'s concept of model
##' "step" and we do this by mapping a date such as `2020-01-01` into
##' the number of days after 2019-12-31 (1 here, with the 1st January 2020
##' being day 1). We call this integer number a "covid_multi_strain date".
##'
##' There are several related functions here
##'
##' * `covid_multi_strain_date` converts its argument into an R `Date` object,
##'   then applies this tranformation. If the argument is not a `Date`
##'   object or a string representing one, an error will be thrown.
##'
##' * `covid_multi_strain_date_to_date` does the reverse conversion to
##'   `covid_multi_strain_date`, converting an integer covid_multi_strain date into an R
##'   `Date`
##'
##' * `as_covid_multi_strain_date` does the same conversion as `covid_multi_strain_date`
##'   but will assume that an integer *already* represents a covid_multi_strain
##'   date and will return it unmodified rather than erroring.
##'
##' * `as_date` does a string to date conversion, using [as.Date()]
##'   but requiring the dates are in ISO 8601 (YYYY-MM-DD) format (it
##'   is a helper that avoids conversion to `NA`, instead throwing an
##'   error)
##'
##' @title Date handling for covid_multi_strain
##'
##' @param date A Date object, or something that can be converted to
##'   one, or a "covid_multi_strain date"; see Details
##'
##' @return An integer, being the number of days into 2020
##' @export
##' @examples
##' # Convert dates into covid_multi_strain dates:
##' covid_multi_strain::covid_multi_strain_date("2020-07-20")
##' covid_multi_strain::covid_multi_strain_date(c("2020-07-20", "2020-10-01"))
##'
##' # Reverse the conversion:
##' covid_multi_strain::covid_multi_strain_date_as_date(1)
##' covid_multi_strain::covid_multi_strain_date_as_date(c(61, 275))
##'
##' # Double conversion not possible with covid_multi_strain_date...
##' try(covid_multi_strain::covid_multi_strain_date(61))
##' # ...but allowed with as_covid_multi_strain_date
##' covid_multi_strain::as_covid_multi_strain_date(61)
##'
##' # Strict date conversion with as_date
##' covid_multi_strain::as_date("2020-03-01")
##' try(covid_multi_strain::as_date("03-01-2020"))
covid_multi_strain_date <- function(date) {
    days_into_2020 <- as.numeric(as_date(date) - as_date("2019-12-31"))
    if (any(days_into_2020 < 0)) {
        stop("Negative dates, covid_multi_strain_date likely applied twice")
    }
    days_into_2020
}


##' @export
##' @rdname covid_multi_strain_date
covid_multi_strain_date_as_date <- function(date) {
    assert_covid_multi_strain_date(date)
    as_date("2019-12-31") + date
}


##' @export
##' @rdname covid_multi_strain_date
as_covid_multi_strain_date <- function(date) {
    if (is.character(date) || is_date(date)) {
        covid_multi_strain_date(as_date(date))
    } else {
        assert_covid_multi_strain_date(date)
    }
}

##' @export
##' @rdname covid_multi_strain_date
as_date <- function(date) {
    if (is_date(date)) {
        return(date)
    }
    if (!all(grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", date))) {
        stop("Expected ISO dates or R dates - please convert")
    }
    as.Date(date)
}


assert_covid_multi_strain_date <- function(date) {
    if (!is.numeric(date)) {
        stop("'date' must be numeric - did you forget covid_multi_strain_date()?")
    }
    if (any(date < 0)) {
        stop("Negative dates, covid_multi_strain_date likely applied twice")
    }
    date
}

is_date <- function(x) {
    inherits(x, "Date")
}