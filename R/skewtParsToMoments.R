#' @title Conversion between parametrizations of a skew-t distribution
#'
#' @description Compute the mean, standard deviation, skewness
#' and excess kurtosis of the skew-t distribution from the parameters
#'
#' @param dp a skew-t direct parameter
#'
#' @return A vector of the mean, standard deviation, skewness
#' and excess kurtosis
#'
#' @author Chindhanai Uthaisaad
#'
#' @examples
#' require("sn")
#' data("Dreturns")
#' dp <- st.mple(y = Dreturns, penalty = "Qpenalty")$dp
#' skewtParsToMoments(dp)
#'
#' @export

skewtParsToMoments <- function(dp) {
  t <- numeric(length(dp))
  if (dp[4] == Inf)
    return(sn.cumulants(dp[1], dp[2], dp[3]))

  delta <- dp[3]/sqrt(1+dp[3]^2)

  if (dp[4] > 343){
    b_nu <- sqrt(2/pi)
  } else if (dp[4] > 1) {
    b_nu <- gamma((dp[4]-1)/2)/sqrt(pi)/gamma(dp[4]/2)*sqrt(dp[4])
  } else {
    stop("Your degree of freedom is less than 1: the moments do not edpist")
  }

  t[1] <- dp[1] + dp[2] * delta * b_nu

  if (dp[4] > 2) {
    sigma_z <- sqrt(dp[4]/(dp[4]-2) - (delta*b_nu)^2)
    t[2] <- dp[2]*sigma_z
  } else {
    sigma_z <- t[2] <- Inf
  }

  if (dp[4] > 3) {
    t[3] <- (delta*b_nu)/(sigma_z^1.5) * (dp[4]*(3-delta^2)/(dp[4]-3) - (3*dp[4]/(dp[4]-2)) + 2*(delta*b_nu)^2)
  } else {
    t[3] <- Inf
  }

  if (dp[4] > 4) {
    t[4] <- (3*dp[4]^2/((dp[4]-2)*(dp[4]-4)) - 4*(delta*b_nu)^2*dp[4]*(3-delta^2)/(dp[4]-3) + 6*(delta*b_nu)^2*dp[4]/(dp[4]-2) - 3*(delta*b_nu)^4)/sigma_z^4-3
  } else {
    t[4] <- Inf
  }

  t
}
