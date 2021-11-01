# Creation of our bulk:
bulk <- new.env()

bulk$rho_vap <- 0.025
bulk$rho_liq <- 0.51
bulk$d <- 1.1244132

bulk$eta <- function(x) 
{
  #' Computes the packing fraction
  return(pi * x * bulk$d^3 / 6)
}

bulk$mu <-function(x) {
  #' Computes the chemical potential needed in a grand-canonical setting
  eta <- bulk$eta(x)

  b_mu <- log(x) - log(1-eta) 
  b_mu <- b_mu + 0.5* eta * ( 5*eta^2 - 13*eta + 14)/((1-eta)^3) 
  b_mu <- b_mu - 8.749726*x
  
  return(b_mu)
}

bulk$p <- function(x, mu) {
  #' Computes the value of the (negative) pressure (delta-omega)
  eta <- bulk$eta(x)

  omega <- x*log(x) - x - x*log(1-eta) 
  omega <- omega + (3/2) * x * (eta * (2-eta) )/( (1-eta)^2 )
  omega <- omega - mu * x - 4.374863 *x^2
 
   return(omega)
}

# Creation of our cluster object:
cluster <- new.env()

#' Defines the minimum size:
cluster$r_min <- 0.95

#' Defines the range of densities to be studied:
cluster$rho <- seq(bulk$rho_vap, bulk$rho_liq, length=100)

cluster$V <- function(r) {
  #' Computes the volume of a cluster:
  return(4 * pi/3 * r^3) 
}

cluster$S <- function(r) {
  #' Computes the surface of a cluster:
  return(4 * pi * r^2)
}

cluster$radius <- function(bulk) {
  #' Computes the radius of a cluster:

  rho <- cluster$rho
  rho_v <- bulk$rho_vap
  rho_l <- bulk$rho_liq
  
  y <- (1/(rho + rho_v))^0.7 + (-1/(rho - rho_l))^0.8
  y <- y - min(y, na.rm = T) + cluster$r_min
  return(y)
}

cluster$W <- function(sfe) {
  #' Computes the work of cluster formation:

  rho_vap <- bulk$rho_vap
  mu_vap <- bulk$mu(rho_vap)
  rho <- cluster$rho
  
  p <- bulk$p
  V <- cluster$V
  S <- cluster$S
  
  radii <- cluster$radius(bulk)

  return(V(radii) * (p(rho, mu_vap)-p(rho_vap, mu_vap)) + sfe * S(radii) * (rho - rho_vap)^2)
}

cluster$f_att <- function(bulk) {
  #' Computes the attachment rate (metric):

  rho <- cluster$rho
  
  rho_v <- bulk$rho_vap
  rho_l <- bulk$rho_liq
  R <- cluster$radius
  
  r <- R(rho+0.1, rho_v, rho_l)
  
  g_xx <- (4*pi/45)*r^5/rho #+ (4*pi/45)*(r^5 *(1-r)^2/(rho_0-rho * r^3)/(1+r+r^2))*r^3
  g_xr <- (4*pi/30)*(((rho-rho_v)/(rho_v-rho*r^3))*r^7) #*(1-r)/((1+r+r^2)^2))*r^3
  g_rr <- (4*pi/5)*(((rho-rho_v)/(rho_v-rho*r^3))*r^6) # /((1+r+r^2)^3))*r^3
  return(g_xx^2 + 2*g_xr*c(diff(r), 0) + g_rr*c(diff(r), 0)^2)
}


# The set of radii associated to the range of cluster densities under study
radii <- cluster$radius(bulk = bulk)

plot(cluster$rho, radii, "o", ylim=c(0,8))
abline(h=cluster$r_min, lty=2)
abline(v=bulk$rho_liq, lty=2)

# Plot omega:

mu_vap <- bulk$mu(bulk$rho_vap)
plot(x = cluster$rho, 
     y = cluster$W(sfe = 0.000001), 
     xlim=c(bulk$rho_vap, bulk$rho_liq),
     ylim=c(-1,10),
     type="o")


# abline(v=rho_0, lty=2)
# abline(v=rho_l, lty=2)
# 
# # Plot g:
# plot(rho, log(f_att(rho, rho_0, rho_l)), type="o")
# lines(rho, radius(rho, rho_0, rho_l), type="o", ylim=c(0,15))
# 
# # Plot omega - (1/2)*log(g)
# plot(rho, omega(rho, rho_0, rho_l, sfe = 16) - 0.5*log(f_att(rho, rho_0, rho_l)), type="o", xlim=c(0, 0.4), ylim=c(-40, 10))
# lines(rho, 0.5*log(f_att(rho, rho_0, rho_l,diff = 10)), type="o")
# 
# #install.packages("plotly")
# library(plotly)
# 
# fig <- plot_ly(#
#   x = rho, #
#   y = radius(rho, rho_0, rho_l), #
#   z = omega(rho, rho_0, rho_l, sfe = 15), #
#   type = 'scatter3d', #
#   mode = 'lines+markers', #
#   line = list(width = 3, color = 2),
#   marker = list(size = 2, color = 2)
# )
# 
# fig <- fig %>% add_markers()
# fig <- fig %>% layout(#
#   scene = list(xaxis = list(title = 'rho', range=c(rho_0,rho_l), autorange="reversed"),
#                yaxis = list(title = 'radius', range=c(0,10)),
#                zaxis = list(title = 'Omega', range=c(0,20)))
#   )
# 
# fig
#
# plot(rho, f_att(rho, rho_0, rho_l), type="o", xlim=c(rho_0,0.4), log="y")
# 
# plot(rho, omega(rho, rho_0, rho_l, sfe = 15), type="o", xlim=c(0, 0.5), ylim=c(-10,20))
# abline(h=0)
# lines(rho,-log(f_att(rho, rho_0, rho_l)))
# 
# plot(rho, omega(rho, rho_0, rho_l, sfe = 25) - log(f_att(rho, rho_0, rho_l)), type="o", xlim=c(0, 0.5), ylim=c(-10, 10))
# lines(rho,-log(f_att(rho, rho_0, rho_l)))
# abline(h=0)
# 
# plot(
#   rho, 
#   exp(-(omega(rho, rho_0, rho_l, sfe = 15) - 0.5*log(f_att(rho, rho_0, rho_l)))), 
#   xlim=c(rho_0,0.4), 
#   type='o',
#   log = "y",
#   ylim=c(1e-2,1e3),
# )