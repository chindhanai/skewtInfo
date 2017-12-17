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


stObsInfoMat <- function(n, xi, omega, alpha, nu) {

  r <- rst(n, xi, omega, alpha, nu)
  S11 <- mean(sapply(r, S_xixi,xi, omega, alpha, nu))
  S12 <- mean(sapply(r, S_xiomega,xi, omega, alpha, nu))
  S13 <- mean(sapply(r, S_xialpha,xi, omega, alpha, nu))
  S14 <- mean(sapply(r, S_xinu,xi, omega, alpha, nu))
  S22 <- mean(sapply(r, S_omegaomega,xi, omega, alpha, nu))
  S23 <- mean(sapply(r, S_omegaalpha,xi, omega, alpha, nu))
  S24 <- mean(sapply(r, S_omeganu,xi, omega, alpha, nu))
  S33 <- mean(sapply(r, S_alphaalpha,xi, omega, alpha, nu))
  S34 <- mean(sapply(r, S_alphanu,xi, omega, alpha, nu))
  S44 <- mean(sapply(r, S_nunu,xi, omega, alpha, nu))

  return(matrix(c(S11, S12, S13, S14,
                  S12, S22, S23, S24,
                  S13, S23, S33, S34,
                  S14, S24, S34, S44),
                nrow = 4, ncol = 4))
}
