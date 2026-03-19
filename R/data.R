#' Microglia Development Index
#'
#' A gene signature for scoring microglial developmental maturation state.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{direction}{Direction of regulation: "UP" or "DOWN"}
#'   \item{valence}{Magnitude of change}
#' }
"mdi"

#' LPS Inflammatory Microglia Index
#'
#' A gene signature for scoring LPS-induced microglial inflammatory activation.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{direction}{Direction of regulation: "UP" or "DOWN"}
#'   \item{valence}{Magnitude of change}
#' }
"lps_index"

#' Fast-Spiking Interneuron Index
#'
#' A gene signature for scoring fast-spiking interneuron transcriptional identity.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Gene}{Gene symbol}
#'   \item{direction}{Direction of regulation: "UP" or "DOWN"}
#'   \item{valence}{Magnitude of change}
#' }
"fsi_index"
