#' Body Fat of Little Brown Bat
#'
#' A dataset containing the percent body fat of little brown bats sampled in
#' Aeolus Cave, located in East Dorset, Vermont, in the North Eastern United States.
#'
#' @usage data(bodyfat_Aeolus)
#' @format A data frame with 159 rows and 4 variables:
#' \describe{
#'   \item{sex}{sex of the sampled bat, \code{M} for masculine and \code{F} for
#'        feminine.}
#'   \item{percentfat}{percent body fat.}
#'   \item{year}{year the bat was sampled.}
#'   \item{days}{hibernation time, defined as days since the fall equinox.}
#' }
#' @details The complete dataset was collected by Cheng et al. (2019) in five
#'     different regions of the United States. \code{bodyfat_Aeolus} is the subset
#'     of the data collected in Aeolus Cave, located in East Dorset, Vermont,
#'     in the North Eastern United States.
#' @keywords datasets
#' @source \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.sh487nh}
#' @references Cheng TL, Gerson A, Moore MS, et al. (2019) Higher fat stores
#'  contribute to persistence of little brown bat populations with white-nose
#'  syndrome. \emph{J Anim Ecol}. 88:591-600.
#'  https://doi.org/10.1111/1365-2656.12954
#' @examples
#' data("bodyfat_Aeolus")
#' # Model with zeta = 2
#' fit <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
#'              family = "PE", zeta = 2)
#' summary(fit)
"bodyfat_Aeolus"
