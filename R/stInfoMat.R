#'
#' @title Skew-t Information Matrix Computation
#'
#' @description Analytically compute skew-t observed information matrix, and
#' compute expected information matrix numerically.
#'
#' @param y a vector of skew-t random samples. It is required when
#' computing the observed information matrix.
#' @param dp a vector of 4 elements of the direct skew-t parameters.
#' @param type a character string with the type of information matrix
#' to be computed
#'
#' @return
#' A list of 4 elements:
#' \item{dp}{the direct skew-t parameters used in the computation}
#' \item{stInfoMat}{the expected information matrix when type = "expected",
#' and the observed information matrix when type = "observed"}
#' \item{SEMat}{the asymptotic standard errors of the skew-t parameters
#' when type = "expected", and the element-wise standard error of the
#' observed information matrix in the case oftype = "observed".}
#' \item{type}{the type of information matrix from the computation}
#'
#' @author Chindhanai Uthaisaad
#'
#' @examples
#' data("cusumData")
#' results = cusumActMgr(portfolioName = "Parvest", benchmarkName = "RUS2500", data = cusumData)
#' chartCusum(results, which = 1)
#' chartCusum(results, which = c(1,2))
#' @export



stInfoMat <- function(y=NULL, dp, type = c("expected", "observed")) {

  if (!(type %in% c("observed", "expected"))) {
    stop("Invalid arg: type should be either 'observed' or 'expected'.")
  }

  if (type == "observed" && is.null(y)) {
    stop("Invalid arg: the realization vector y is required when type = 'observed'.")
  }

  if (!is.vector(dp) || !is.numeric(dp)) {
    stop("Invalid arg: the vector of skew-t parameters must be a numerical vector")
  }

  if (length(dp) < 4) {
    stop("Invalid arg: the vector of skew-t parameters must be of length 4")
  }

  type = type[1]
  # Define the parameters
  xi=dp[1]; omega=dp[2]; alpha=dp[3]; nu=dp[4]

  if (type == "observed") {
    # Apply the second derivative formulas to the sample vector y
    xixi <- sapply(y, S_xixi,xi, omega, alpha, nu)
    xiomega <- sapply(y, S_xiomega,xi, omega, alpha, nu)
    xialpha <- sapply(y, S_xialpha,xi, omega, alpha, nu)
    xinu <- sapply(y, S_xinu,xi, omega, alpha, nu)
    omegaomega <- sapply(y, S_omegaomega,xi, omega, alpha, nu)
    omegaalpha <- sapply(y, S_omegaalpha,xi, omega, alpha, nu)
    omeganu <- sapply(y, S_omeganu,xi, omega, alpha, nu)
    alphaalpha <- sapply(y, S_alphaalpha,xi, omega, alpha, nu)
    alphanu <- sapply(y, S_alphanu,xi, omega, alpha, nu)
    nunu <- sapply(y, S_nunu,xi, omega, alpha, nu)

    # Standard Error Matrix
    SE11 <- sd(xixi)
    SE12 <- sd(xiomega)
    SE13 <- sd(xialpha)
    SE14 <- sd(xinu)
    SE22 <- sd(omegaomega)
    SE23 <- sd(omegaalpha)
    SE24 <- sd(omeganu)
    SE33 <- sd(alphaalpha)
    SE34 <- sd(alphanu)
    SE44 <- sd(nunu)
    SEMat <- (matrix(c(SE11, SE12, SE13, SE14,
                        SE12, SE22, SE23, SE24,
                        SE13, SE23, SE33, SE34,
                        SE14, SE24, SE34, SE44),
                      nrow = 4, ncol = 4))

    # Observed Information Matrix
    S11 <- mean(xixi)
    S12 <- mean(xiomega)
    S13 <- mean(xialpha)
    S14 <- mean(xinu)
    S22 <- mean(omegaomega)
    S23 <- mean(omegaalpha)
    S24 <- mean(omeganu)
    S33 <- mean(alphaalpha)
    S34 <- mean(alphanu)
    S44 <- mean(nunu)

    stInfoMat <- (-matrix(c(S11, S12, S13, S14,
                             S12, S22, S23, S24,
                             S13, S23, S33, S34,
                             S14, S24, S34, S44),
                           nrow = 4, ncol = 4))

  } else {

    I11 <- I_xixi(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I12 <- I_xiomega(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I13 <- I_xialpha(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I14 <- I_xinu(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I22 <- I_omegaomega(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I23 <- I_omegaalpha(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I24 <- I_omeganu(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I33 <- I_alphaalpha(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I34 <- I_alphanu(xi = xi, omega = omega, alpha = alpha, nu = nu)
    I44 <- I_nunu(xi = xi, omega = omega, alpha = alpha, nu = nu)

    stInfoMat <- matrix(c(I11, I12, I13, I14,
                          I12, I22, I23, I24,
                          I13, I23, I33, I34,
                          I14, I24, I34, I44),
                        nrow = 4, ncol = 4)
    SEMat <- sqrt(diag(solve(stInfoMat)))

  }

  return(list(dp = dp,
              stInfoMat = stInfoMat,
              SEMat = SEMat,
              type = type))
}

# Observed Information Matrix functions -----------------------
gamma_integrand <- function(u,nu){
  dt(u,df=nu+1)*(((nu+2)*u^2)/((nu+1)*(nu+1+u^2))-log(1+u^2/(nu+1)))
}

Gamma <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  value <- integrate(gamma_integrand,-Inf, varxi, nu = nu, rel.tol = 1e-10)$val
  return(value)
}

beta_integrand <- function(u,nu){
  dt(u,df=nu+1)*(((nu+2)*u^2)/((nu+1)*(nu+1+u^2))-log(1+u^2/(nu+1)))^2
}

Beta <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  value <- integrate(beta_integrand,-Inf, varxi, nu = nu, rel.tol = 1e-10)$val
  return(value)
}

delta_integrand <- function(u, nu){
  dt(u,df=nu+1)*((nu*u^2-2*nu-2)*u^2)/((nu+1)*(nu+1+u^2))^2
}

Delta <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  value <- integrate(delta_integrand, -Inf, varxi, nu=nu, rel.tol=1e-10)$val
  return(value)
}


S_z <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  value <- (alpha*tau*nu*w)/(nu+z^2) - tau^2*z
  return(value)
}


S_zz <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  w_z <- -((nu*(nu+2)*alpha^2*z*w)/((nu+z^2+alpha^2*z^2)*(nu+z^2)))-((nu*alpha*tau*w^2)/(nu+z^2))
  value <- (2*tau^2*z^2)/(nu+z^2)-tau^2-(3*alpha*tau*nu*z*w)/(nu+z^2)^2+(alpha*tau*nu*w_z)/(nu+z^2)
  return(value)
}


S_zalpha <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  w_alpha <- -(nu+2)*alpha*z^2*w/(nu+z^2+alpha^2*z^2)-z*tau*(w^2)
  value <- nu*tau*(w+alpha*w_alpha)/(nu+z^2)
  return(value)
}


S_znu <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  w_nu <-   w / 2 * (((nu+2) * alpha^2 * z^2)/((nu + z^2 + alpha^2*z^2)*(nu + z^2)) - log(1 + alpha^2 * z^2 / (nu + z^2)) - Gamma(y, xi, omega, alpha, nu)/pt(varxi, df = nu + 1)) + (alpha * z * (1 - z^2)* w^2)/(2 * tau * (nu + z^2)^2)

  value <- (z*(1-z^2))/(nu+z^2)^2+(alpha*(nu*(3*z^2-1)+2*z^2)*w)/(2*tau*(nu+z^2)^3)+(alpha*tau*nu*w_nu)/(nu+z^2)
  return(value)
}

S_xixi <- function(y, xi, omega, alpha, nu){
  value <- S_zz(y, xi, omega, alpha, nu)/omega^2
  return(value)
}

S_xiomega <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  value <- z/omega^2*S_zz(y, xi, omega, alpha, nu) + S_z(y, xi, omega, alpha, nu)/omega^2
  return(value)
}

S_xialpha <- function(y, xi, omega, alpha, nu){
  value <- -S_zalpha(y, xi, omega, alpha, nu)/omega
  return(value)
}

S_omegaomega <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  value <- 1/omega^2+z^2/omega^2*S_zz(y, xi, omega, alpha, nu)+2*z/omega^2*S_z(y, xi, omega, alpha, nu)
  return(value)
}

S_omegaalpha <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  value <- -z/omega * S_zalpha(y, xi, omega, alpha, nu)
  return(value)
}

S_omeganu <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  value <- -z/omega * S_znu(y, xi, omega, alpha, nu)
  return(value)
}

S_alphaalpha <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  w_alpha <- -((nu+2)*alpha*z^2*w)/(nu+z^2+alpha^2*z^2)-z*tau*w^2
  value <- z*tau*w_alpha
  return(value)
}

S_nunu <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  w_nu <- w / 2 * (((nu+2) * alpha^2 * z^2)/((nu + z^2 + alpha^2*z^2)*(nu + z^2)) - log(1 + alpha^2 * z^2 / (nu + z^2)) - Gamma(y, xi, omega, alpha, nu)/pt(varxi, df = nu + 1)) + (alpha * z * (1 - z^2)* w^2)/(2 * tau * (nu + z^2)^2)

  return((trigamma(nu/2 + 1) - trigamma(nu/2)) / 4 + (2 * nu^2 + 2*nu + 1)/(2 * nu^2 * (nu + 1)^2) + z^2/(2 * nu * (nu + z^2))
         -(z^2 * (nu^2 + 2*nu + z^2))/(2 * nu^2 * (nu+z^2)^2) - (alpha * z * (z^2 -1) * (z^2 + 4*nu + 3) * w)/(4 * tau * (nu + 1) * (nu + z^2)^3)
         + (alpha * z * (1 - tau^2) * w_nu)/(2 * tau * (nu + z^2)) - Gamma(y, xi, omega, alpha, nu)^2 / (4 * pt(varxi, df = nu+1)^2)
         - (alpha * z * (z^2 -1) * Gamma(y, xi=xi, omega, alpha, nu) * w)/(4 * pt(varxi, df = nu+1) * tau * (nu + z^2)^2)
         + (2 * Delta(y, xi, omega, alpha, nu) + Beta(y, xi, omega, alpha, nu))/(4 * pt(varxi, df = nu+1))
         + (((nu+2) * alpha^2 * z^2)/((nu+1) * (nu + z^2 + alpha^2*z^2)) - log(1 + alpha^2 * z^2 / (nu + z^2))) * (alpha * z * (z^2 - 1) * w)/(4 * tau * (nu + z^2)^2))

}

S_xinu <- function(y, xi, omega, alpha, nu){
  value <- -S_znu(y, xi, omega, alpha, nu)/omega
  return(value)
}

S_alphanu <- function(y, xi, omega, alpha, nu){
  z <- (y-xi)/omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha*z*tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)
  w_nu <- w / 2 * (((nu+2) * alpha^2 * z^2)/((nu + z^2 + alpha^2*z^2)*(nu + z^2)) - log(1 + alpha^2 * z^2 / (nu + z^2)) - Gamma(y, xi, omega, alpha, nu)/pt(varxi, df = nu + 1)) + (alpha * z * (1 - z^2)* w^2)/(2 * tau * (nu + z^2)^2)

  return((z*(z^2-1)*w)/(2*tau*(nu+z^2)^2)+z*tau*w_nu)
}

# Expected Information Matrix Functions------------------
g_xixi <- function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (z * tau ^ 2 / omega - alpha * tau * nu * w / (omega * (nu + z ^ 2))) ^ 2 * dst(y, xi = xi, omega = omega, alpha = alpha, nu = nu)
}

I_xixi <- function(xi, omega, alpha, nu) {
  integrate(g_xixi, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_omegaomega <- function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (-1 / omega + z^2 * tau^2 / omega  + alpha * w * (tau * z^3/ (omega * (nu+z^2)) - z * tau/omega)) ^ 2 * dst(y, xi = xi, omega = omega, alpha = alpha, nu = nu)
}

I_omegaomega <- function(xi, omega, alpha, nu) {
  integrate(g_omegaomega, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_alphaalpha <- function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (z * tau * w)^2 * dst(y, xi = xi, omega = omega, alpha = alpha, nu = nu)
}

I_alphaalpha <- function(xi, omega, alpha, nu) {
  integrate(g_alphaalpha, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_xiomega <- function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (z * tau ^ 2 / omega - alpha * tau * nu * w / (omega * (nu + z ^ 2))) * (-1 / omega + z^2 * tau^2 / omega  + alpha * w * (tau * z^3/ (omega * (nu+z^2)) - z * tau/omega)) * dst(y, xi = xi, omega = omega, alpha = alpha, nu = nu)
}

I_xiomega <- function(xi, omega, alpha, nu) {
  integrate(g_xiomega, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_xialpha <- function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (z * tau ^ 2 / omega - alpha * tau * nu * w / (omega * (nu + z ^ 2))) * (z * tau * w) * dst(y, xi, omega, alpha, nu)
}

I_xialpha <- function(xi, omega, alpha, nu) {
  integrate(g_xialpha, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_omegaalpha <- function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (-1 / omega + z^2 * tau^2 / omega  + alpha * w * (tau * z^3/ (omega * (nu+z^2)) - z * tau/omega)) * (z * tau * w) * dst(y, xi, omega, alpha, nu)
}

I_omegaalpha <- function(xi, omega, alpha, nu) {
  integrate(g_omegaalpha, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

gamma_integrand <- function(u, nu){
  dt(u,df = nu + 1) * (((nu+2) * u^2 )/((nu+1)*(nu + 1 + u^2))-log(1 + u^2/(nu + 1)))
}

g_nunu <- Vectorize(function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (0.5 * (digamma(1+nu/2) - digamma(nu/2) - (2*nu + 1)/(nu * (nu + 1))
          - log(1 + z^2 / nu) + z^2 * tau^2 / nu + (alpha * z * (z^2 - 1) * w)/((nu + z^2)^2 * tau)
          + sapply(y, function(y) { # Gamma
            integrate(function(u, nu) gamma_integrand(u, nu), lower = -Inf, upper = varxi, nu = nu, rel.tol = 1e-12)$value
          }) / pt(varxi, nu + 1))) ^ 2 * dst(y, xi, omega, alpha, nu)
})

I_nunu <- function(xi, omega, alpha, nu) {
  integrate(g_nunu, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_xinu <- Vectorize(function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (z * tau ^ 2 / omega - alpha * tau * nu * w / (omega * (nu + z ^ 2))) * (0.5 * (digamma(1+nu/2) - digamma(nu/2) - (2*nu + 1)/(nu * (nu + 1))
                                                                                  - log(1 + z^2 / nu) + z^2 * tau^2 / nu + (alpha * z * (z^2 - 1) * w)/((nu + z^2)^2 * tau)
                                                                                  + sapply(y, function(y) { # Gamma
                                                                                    integrate(function(u, nu) gamma_integrand(u, nu), lower = -Inf, upper = varxi, nu = nu, rel.tol = 1e-12)$value
                                                                                  }) / pt(varxi, nu + 1))) * dst(y, xi, omega, alpha, nu)
})

I_xinu <- function(xi, omega, alpha, nu) {
  integrate(g_xinu, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_omeganu <- Vectorize(function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (-1 / omega + z^2 * tau^2 / omega  + alpha * w * (tau * z^3/ (omega * (nu+z^2)) - z * tau/omega)) * (0.5 * (digamma(1+nu/2) - digamma(nu/2) - (2*nu + 1)/(nu * (nu + 1))
                                                                                                              - log(1 + z^2 / nu) + z^2 * tau^2 / nu + (alpha * z * (z^2 - 1) * w)/((nu + z^2)^2 * tau)
                                                                                                              + sapply(y, function(y) { # Gamma
                                                                                                                integrate(function(u, nu) gamma_integrand(u, nu), lower = -Inf, upper = varxi, nu = nu, rel.tol = 1e-12)$value
                                                                                                              }) / pt(varxi, nu + 1))) * dst(y, xi, omega, alpha, nu)
})

I_omeganu <- function(xi, omega, alpha, nu) {
  integrate(g_omeganu, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}

g_alphanu <- Vectorize(function(y, xi, omega, alpha, nu) {
  z <- (y - xi) / omega
  tau <- sqrt((nu+1)/(z^2+nu))
  varxi <- alpha * z * tau
  w <- dt(varxi, nu+1)/pt(varxi, nu+1)

  (z * tau * w) * (0.5 * (digamma(1+nu/2) - digamma(nu/2) - (2*nu + 1)/(nu * (nu + 1))
                          - log(1 + z^2 / nu) + z^2 * tau^2 / nu + (alpha * z * (z^2 - 1) * w)/((nu + z^2)^2 * tau)
                          + sapply(y, function(y) { # Gamma
                            integrate(function(u, nu) gamma_integrand(u, nu), lower = -Inf, upper = varxi, nu = nu, rel.tol = 1e-12)$value
                          }) / pt(varxi, nu + 1))) * dst(y, xi, omega, alpha, nu)
})

I_alphanu <- function(xi, omega, alpha, nu) {
  integrate(g_alphanu, lower = -Inf, upper = Inf,
            xi=xi, omega=omega, alpha=alpha, nu=nu, rel.tol = 1e-10)$val
}
