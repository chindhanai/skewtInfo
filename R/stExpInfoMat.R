# Numerical integration for information matrix of Azzalini's skew-t distribution

# S_xi <-  z * tau ^ 2 / omega - alpha * tau * nu * w / (omega * (nu + z ^ 2))
# S_omega <- -1 / omega + z^2 * tau^2 / omega  + alpha * w * (tau * z^3/ (omega * (nu+z^2)) - z * tau/omega)
# S_alpha <- z * tau * w

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

# S_nu <- 0.5 * (digamma(1+nu/2) - digamma(nu/2) - (2*nu + 1)/(nu * (nu + 1))
#         - log(1 + z^2 / nu) + z^2 * tau^2 / nu + (alpha * z * (z^2 - 1) * w)/((nu + z^2)^2 * tau)
#         + sapply(y, function(y) { # Gamma
# integrate(function(u, nu) gamma_integrand(u, nu), lower = -Inf, upper = varxi, nu = nu, rel.tol = 1e-12)$value
# }) / pt(varxi, nu + 1))

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

# S_xi <-  z * tau ^ 2 / omega - alpha * tau * nu * w / (omega * (nu + z ^ 2))
# S_omega <- -1 / omega + z^2 * tau^2 / omega - alpha * z * tau * nu * w / (w * (nu + z ^ 2))
# S_alpha <- z * tau * w

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

#Define the function to compute the information matrix
stExpInfoMat <- function(xi, omega, alpha, nu){
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

  matrix(c(I11, I12, I13, I14,
           I12, I22, I23, I24,
           I13, I23, I33, I34,
           I14, I24, I34, I44), nrow = 4, ncol = 4)
}


