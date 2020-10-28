#' Generalized pivotal quantities
#'
#' Simulates from the generalized pivotal quantities.
#'
#' @param fit an \code{\link{aov1r}} object
#' @param n number of simulations
#'
#' @return The simulations in a list.
#'
#' @importFrom stats rchisq
#' @export
#'
#' @examples
#' dat <- simAOV1R(I=20, J=5, mu=10, sigmab=1, sigmaw=1)
#' fit <- aov1r(y ~ group, data=dat)
#' pivsims <- pivotal(fit)
#' pivsims$G_sigma2tot <- pivsims$G_sigma2b + pivsims$G_sigma2w
#' # Generalized confidence intervals:
#' lapply(pivsims, quantile, probs = c(0.025, 0.975))
#' # compare with the frequentist confidence intervals:
#' confint(fit, SDs = FALSE)
#' # Generalized prediction interval:
#' with(
#'   pivsims,
#'   quantile(rnorm(length(G_mu), G_mu, sqrt(G_sigma2tot)),
#'            probs = c(0.025, 0.975))
#' )
#' # compare with the frequentist prediction interval:
#' predict(fit)
pivotal <- function(fit, n=10000){
  I <- fit[["Design"]][["I"]]
  J <- fit[["Design"]][["Jh"]]
  N <- fit[["Design"]][["N"]]
  ssb <- fit[["Sums of squares"]][["ssb"]]
  ssw <- fit[["Sums of squares"]][["ssw"]]
  Z <- rnorm(n)
  U2b <- rchisq(n, I-1)
  U2w <- rchisq(n, N-I)
  list(
    G_mu = fit[["grandMean"]] - Z/sqrt(U2b)*sqrt(ssb/I/J),
    G_sigma2b = 1/J*(ssb/U2b - ssw/U2w),
    G_sigma2w = ssw/U2w
  )
}

# pivotal0 <- function(fit, Z, U2b, U2w){
#   I <- fit[["Design"]][["I"]]
#   J <- fit[["Design"]][["Jh"]]
#   ssb <- fit[["Sums of squares"]][["ssb"]]
#   ssw <- fit[["Sums of squares"]][["ssw"]]
#   list(
#     G_mu = fit[["grandmean"]] - Z/sqrt(U2b)*sqrt(ssb/I/J),
#     G_sigma2b = 1/J*(ssb/U2b - ssw/U2w),
#     G_sigma2w = ssw/U2w
#   )
# }
