#' Clone information dataset in Wide Format
#'
#' A dataset from a simulation taken from "Niche engineering drives early passage through an immune bottleneck in progression to colorectal cancer"
#'
#' @format A data frame with 6591 rows and 43 variables:
#' \describe{
#'   \item{master_random_seed}{Simulation metadata}
#'   \item{macrophage_suppression}{Simulation metadata}
#'   \item{macrophage_growth_benefit}{Simulation metadata}
#'   \item{pdl1_protection}{Simulation metadata}
#'   \item{tumor_random_seed}{Simulation metadata}
#'   \item{tumor_lineage}{Simulation metadata}
#'   \item{tumor_status}{Simulation metadata}
#'   \item{onset}{Simulation metadata}
#'   \item{parent}{Simulation metadata}
#'   \item{og_parent}{Simulation metadata}
#'   \item{clone}{Simulation metadata}
#'   \item{type}{Simulation metadata}
#'   \item{division_rate}{Simulation metadata}
#'   \item{n_driver_mutations}{Simulation metadata}
#'   \item{new_antigenicity}{Simulation metadata}
#'   \item{og_antigenicity}{Simulation metadata}
#'   \item{final_antigenicity}{Simulation metadata}
#'   \item{macrophage}{Simulation metadata}
#'   \item{pdl1}{Simulation metadata}
#'   \item{immune_susceptibility}{Simulation metadata}
#'   \item{detected}{Simulation metadata}
#'   \item{size_at_time_of_detection}{Simulation metadata}
#'   \item{X2227}{Timepoint}
#'   \item{X2257}{Timepoint}
#'   \item{X2287}{Timepoint}
#'   \item{X2317}{Timepoint}
#'   \item{X2347}{Timepoint}
#'   \item{X2377}{Timepoint}
#'   \item{X2407}{Timepoint}
#'   \item{X2437}{Timepoint}
#'   \item{X2467}{Timepoint}
#'   \item{X2497}{Timepoint}
#'   \item{X2527}{Timepoint}
#'   \item{X2557}{Timepoint}
#'   \item{X2587}{Timepoint}
#'   \item{X2617}{Timepoint}
#'   \item{X2647}{Timepoint}
#'   \item{X2677}{Timepoint}
#'   \item{X2707}{Timepoint}
#'   \item{X2737}{Timepoint}
#'   \item{X2767}{Timepoint}
#'   \item{X2797}{Timepoint}
#'   \item{X2808}{Timepoint}
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/623959v2}
#' @examples
#' if (require("EvoFreq")) {
#' example.wide
#' }
"example.wide"