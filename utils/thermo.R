library(rootSolve)

fHelmholtz <- function(x, beta, n)
{ #' Helmholtz free energy and its derivative
  p_fac <- pi * d^3 / 6.0
  eta <- p_fac * x
  a_const <- (-4.374863)

  if(n == 0)
  {
    f_val <- x*log(x) - x - x*log(1-eta) +
      1.5*x*eta*(2-eta)*(1-eta)^(-2.0) +
       (beta * a_const*x*x)
  }
  if (n == 1)
  {
    f_val <- log(x)-log(1-eta) +
      0.5 * eta * ( 5 * eta*eta -13*eta +14 )*(1-eta)^(-3.0) +
      (2*beta*a_const*x);
  }

  return(f_val)
}

mu <-function(x, beta)
{ #' Computes the chemical potential needed in a grand-canonical setting

  return(fHelmholtz(x = x, beta = beta, n = 1))
}

omega <- function(x, mu, beta)
{ #' Computes Landau's potential density

  omega <- fHelmholtz(x = x, beta = beta, n = 0) - mu * x
  return(omega)
}

p <- function(x, beta) {
  #' Computes the value of the (negative) pressure (delta-omega)

  mu_val <- mu(x = x, beta = beta)
  return(-omega(x = x, mu = mu_val, beta = beta))
}


eos <- function(x, beta) {
  y <- numeric(2)

  x_v <- x[1]
  x_l <- x[2]

  mu_v <- mu(x = x_v, beta = beta)
  mu_l <- mu(x = x_l, beta = beta)

  p_v <- p(x = x_v, beta = beta)
  p_l <- p(x = x_l, beta = beta)

  y[1] <- p_v - p_l
  y[2] <- mu_v - mu_l

  return(y)
}

find_liq <- function(x, x_v, beta) {

  mu_v <- mu(x = x_v, beta = beta)
  mu_l <- mu(x = x, beta = beta)

  f1 <- mu_l - mu_v

  return(f1)
}


V <- function(r) {
  #' Computes the volume of a cluster:
  return(4 * pi/3 * r^3) 
}

S <- function(r) {
  #' Computes the surface of a cluster:
  return(4 * pi * r^2)
}

r_path <- function(x, rho_v, rho_l) {
  #' Computes the radius of a cluster:

  y <- (1/(x + rho_v))^0.7 + (-1/(x - rho_l))^0.8
  y <- y - min(y, na.rm = T) + r_min

  return(y)
}

W <- function(x, rho_v, rho_l, beta, sfe = 0.5 * 2.16574 ) {
  #' Computes the work of cluster formation:

  r <- r_path(x = x, rho_v = rho_v, rho_l = rho_l)

  mu_vap <- mu(rho_v, beta = beta)
  d_rho <- (x - rho_v)

  O <- (omega(x, mu_vap, beta) - omega(rho_v, mu_vap, beta)) * V(r) +
    beta * sfe * d_rho * d_rho * S(r)

  return(O)
}

f_att <- function(x, rho_v, rho_l) {
  #' Computes the attachment rate (metric):
  r <- r_path(x+0.1, rho_v, rho_l)
  
  g_xx <- (4*pi/45)*r^5/x #+ (4*pi/45)*(r^5 *(1-r)^2/(rho_0-rho * r^3)/(1+r+r^2))*r^3
  g_xr <- (4*pi/30)*(((x-rho_v)/(rho_v-x*r^3))*r^7) #*(1-r)/((1+r+r^2)^2))*r^3
  g_rr <- (4*pi/5)*(((x-rho_v)/(rho_v-x*r^3))*r^6) # /((1+r+r^2)^3))*r^3
  return(g_xx^2 + 2*g_xr*c(diff(r), 0) + g_rr*c(diff(r), 0)^2)
}  

W_eff <- function(x, rho_v, rho_l, beta) {
  w <- W(x = x, rho_v = rho_v, rho_l = rho_l, beta = beta)
  w_k <- log(f_att(x = x, rho_v = rho_v, rho_l = rho_l))
  return( w - 0.5 * w_k )
}

supersat_graph <- function(rho_v, rho_l, rho_c, beta, func, n_points = 500)
{
  x <- seq(rho_v, rho_l, length=n_points)
  x_axis <- x / rho_v
  y_axis <- func(x = x, rho_v = rho_v, rho_l = rho_l, beta = beta)
  
  y_max_idx <- which.max(y_axis) 
  y_max <- y_axis[y_max_idx]
  x_max <- x_axis[y_max_idx]
  
  print(y_max)
  print(x_max)
  
  plot(x_axis,
       y_axis,
       type="o",
       xlim=c(1, x_max*1.05),
       ylim=c(-5, y_max*1.05)
  )
  abline(h=0, lty=2)
  
  for (s in seq(rho_v/rho_c, 5, length=20)) {
    rho_v <- s * rho_c

    print(rho_v)
    print(rho_l)

    coex_liq <- function(x) {return(find_liq(x, x_v = rho_v, beta =  b_beta))}
    root <- multiroot(f = coex_liq, start = 1.01 * x_l.coex)
    rho_l <- root$root
    
    x <- seq(rho_v, rho_l, length=n_points)
    x_axis <- x / rho_v

    func(x = x, rho_v = rho_v, rho_l = rho_l, beta = beta)

    lines(x_axis,
          y_axis,
          lty=2,
          col=s
    )
  }
}




