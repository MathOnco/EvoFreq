#' Clone information dataset in Long Format
#'
#' A dataset from a simulation taken from "Niche engineering drives early passage through an immune bottleneck in progression to colorectal cancer"
#'
#' @format A data frame with 1181 rows and 5 variables:
#' \describe{
#'   \item{parents}{Parent clone}
#'   \item{fitness}{Clone metadata}
#'   \item{drivers}{Driver mutations of clone}
#'   \item{passengers}{Passenger mutations of clone}
#'   \item{clone_id}{ID of the clone}

#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/623959v2}
#' @examples
#' if (require("EvoFreq")) {
#' example.long.edges
#' }
"example.long.edges"