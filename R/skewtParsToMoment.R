skewtTransform <- function(x) {
  t <- numeric(length(x))
  if (nu == Inf) 
    return(sn.cumulants(x[1], x[2], x[3]))
  
  delta <- x[3]/sqrt(1+x[3]^2)
  
  if (x[4] > 343){
    b_nu <- sqrt(2/pi)
  } else if (x[4] > 1) {
    b_nu <- gamma((x[4]-1)/2)/sqrt(pi)/gamma(x[4]/2)*sqrt(x[4])
  } else {
    stop("Your degree of freedom is less than 1: the moments do not exist")
  }
  
  t[1] <- x[1] + x[2] * delta * b_nu
  
  if (x[4] > 2) {
    sigma_z <- sqrt(x[4]/(x[4]-2) - (delta*b_nu)^2)
    t[2] <- x[2]*sigma_z
  } else {
    sigma_z <- t[2] <- Inf
  } 
  
  if (x[4] > 3) {
    t[3] <- (delta*b_nu)/(sigma_z^1.5) * (x[4]*(3-delta^2)/(x[4]-3) - (3*x[4]/(x[4]-2)) + 2*(delta*b_nu)^2)
  } else {
    t[3] <- Inf
  } 
  
  if (x[4] > 4) {
    t[4] <- (3*x[4]^2/((x[4]-2)*(x[4]-4)) - 4*(delta*b_nu)^2*x[4]*(3-delta^2)/(x[4]-3) + 6*(delta*b_nu)^2*x[4]/(x[4]-2) - 3*(delta*b_nu)^4)/sigma_z^4-3
  } else {
    t[4] <- Inf
  } 
  
  t
}