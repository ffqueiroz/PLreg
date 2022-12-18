#' Peru Blank Votes
#'
#' A dataset on the blank votes in the 2006 Peruvian general election.
#'
#' @usage data(PeruVotes)
#' @format A data frame with 194 rows and 2 variables:
#' \describe{
#'   \item{votes}{proportion of blank votes in the 2006 Peruvian general election.}
#'   \item{HDI}{Human Development Index.}
#' }
#' @details The dataset was collected by Bayes et al. (2012) and the response variable
#'     is \code{votes}, proportion of blank votes in the 2006 Peruvian general election.
#'     Bayes et al. (2012) analyzed the influence of the Human Development Index (HDI) on
#'     the proportion of blank votes using a beta rectangular regression model. Lemonte and
#'     Bazan (2015) also analyzed the data using GJS regression models.
#' @keywords datasets
#' @source \url{http://www.undp.org/es/peru}
#' @references Bayes, C., Bazan, J. L. and Garcia, C. (2012). A new robust regression model for proportions.
#'     \emph{Bayesian Analysis}. 7:771-796 \cr \cr
#'     Lemonte, A. J. and Bazan, J. L. (2015). New class of Johnson SB distributions
#'     and its associated regression model for rates and proportions. \emph{Biometrical Journal}. 58:727-746.
#' @examples data("PeruVotes")
#' fitPL <- PLreg(votes ~ HDI | HDI,
#'                data = PeruVotes,
#'                family = "TF",
#'                zeta = 5,
#'                control = PLreg.control(lambda = 1))
#' summary(fitPL)
#' plot(fitPL, type = "standardized")
"PeruVotes"
