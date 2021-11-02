plot_radius <- function(bulk) {
  radii <- r_path(x = bulk$rho, rho_v = bulk$rho_vap, rho_l = bulk$rho_liq)

  plot(bulk$rho, radii, "o", ylim=c(0,8))
  abline(h=bulk$r_min, lty=2)
  abline(v=bulk$rho_liq, lty=2)
  abline(v=bulk$rho_vap, lty=2)
}

plot_omega <- function(bulk)
{
  Omega <- W(x = bulk$rho,
             rho_v = bulk$rho_vap,
             rho_l = bulk$rho_liq,
             beta = bulk$beta)

  cc_idx <- which.max(Omega)
  o_cc <- Omega[cc_idx]
  r_cc <- rho[cc_idx]

  plot(x = bulk$rho,
       y = Omega,
       xlim=c(bulk$rho_vap, 1.1 * r_cc),
       ylim=c(-1, 1.1 * o_cc ),
       type="o")

  points(r_cc, o_cc, col="red", pch=16)
  abline(h = o_cc, lty = 2, col="grey")
  abline(h = 0, lty = 2)
}

plot_g <- function(bulk) {
  g_m <- log(f_att( x = bulk$rho, rho_v = bulk$rho_vap, rho_l = bulk$rho_liq))
  plot(x = bulk$rho,
       y = g_m,
       xlim=c(0,0.6*bulk$rho_liq),
       ylim=c(0,10),
       type="o")
  abline(h=0, lty=2)
}


plot_potentials <- function(bulk) {
  par(mfrow=c(1,2))
  Omega <- W(x = bulk$rho,
             rho_v = bulk$rho_vap,
             rho_l = bulk$rho_liq,
             beta = bulk$beta)

  cc_idx <- which.max(Omega)
  o_cc <- Omega[cc_idx]
  r_cc <- rho[cc_idx]

  plot(x = bulk$rho,
       y = Omega,
       xlim=c(0, 1.1*r_cc),
       ylim=c(-10, 1.1*o_cc),
       type="l",
       lwd = 2
  )
  abline(h=0)
  abline(h = o_cc, lty = 2, col="grey")
  points(x = r_cc, y = o_cc, col = 2, pch = 16)

  Omega_eff <- W_eff(x = bulk$rho,
                     rho_v = bulk$rho_vap,
                     rho_l = bulk$rho_liq,
                     beta = bulk$beta)

  cc_eff_idx <- which.max(Omega_eff)
  r_cc_eff <- rho[cc_eff_idx]
  o_cc_eff <- Omega_eff[cc_eff_idx]

  plot(x = bulk$rho,
       y = Omega_eff,
       xlim=c(0, 1.1*r_cc_eff),
       ylim=c(-10, 1.1*o_cc_eff),
       type="l",
       lwd = 2
  )
  cc_eff_idx <- which.max(Omega_eff)
  abline(h=0)
  abline(h = o_cc_eff, lty = 2, col="grey")
  points(x = r_cc_eff, y = o_cc_eff, col = 2, pch = 16)
}
