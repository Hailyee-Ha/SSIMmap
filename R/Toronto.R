#' Sample polygon data of Toronto
#'
#' A sf(simple feature) containing geometric boundaries of Toronto DAs(Dissemination Area) with their codes.
#'
#' @format A \code{sf} object with 58 rows and 2 variables
#' \describe{
#'   \item{DAUID}{Dissemination Area ID}
#'   \item{CIMD_SDD}{Factor score of CIMD(The Canadian Index of Multiple Deprivation) social deprivation dimension}
#'   \item{PP_SDD}{Principal score of Pampalon social deprivation dimension}
#'   \item{P_commute}{Percentage of households who commute within census subdivision (CSD) of residence }
#'   \item{geometry}{the geometry column for counties(CRS: NAD83)}
#' }
"Toronto"
