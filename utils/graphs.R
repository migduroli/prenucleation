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
  r_cc <- bulk$rho[cc_idx]

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


compute_omega <- function(bulk) {
  Omega <- W(x = bulk$rho,
             rho_v = bulk$rho_vap,
             rho_l = bulk$rho_liq,
             beta = bulk$beta)

  cc_idx <- which.max(Omega)
  o_cc <- Omega[cc_idx]
  r_cc <- bulk$rho[cc_idx]

  return(list(
    omega=Omega,
    w_c=o_cc,
    r_c=r_cc
  ))
}

compute_omega_eff <- function(bulk) {
  Omega_eff <- W_eff(x = bulk$rho,
                     rho_v = bulk$rho_vap,
                     rho_l = bulk$rho_liq,
                     beta = bulk$beta)

  cc_eff_idx <- which.max(Omega_eff)
  r_cc_eff <- bulk$rho[cc_eff_idx]
  o_cc_eff <- Omega_eff[cc_eff_idx]

  return(list(
    omega=Omega_eff,
    w_c=o_cc_eff,
    r_c=r_cc_eff
  ))
}

get_limits <- function(is_first, omega_obj, limits_obj, n) {
  if (is_first) {
    x_lims <- c(0, 1.1*omega_obj$r_c)
    y_lims <- c(-10, 1.1*omega_obj$w_c)
  } else {
    if (n == 1) {
      x_lims <- limits_obj$xlim
      y_lims <- limits_obj$ylim
    } else {
      x_lims <- limits_obj$xlim_eff
      y_lims <- limits_obj$ylim_eff
    }
  }
  return(list(xlim=x_lims, ylim=y_lims))
}

plot_panel <- function(n, bulk_obj, omega_obj, lims_obj, is_first) {

  if (is_first) {
    xlab <- expression(rho)

    if (n == 1) {
      ylab <- expression(paste(beta, Omega))
    }
    else {
      ylab <- expression(paste(beta, Omega[eff]))
    }

    plot(
      x = bulk_obj$rho,
      y = omega_obj$omega,
      xlim = lims_obj$xlim,
      ylim = lims_obj$ylim,
      xlab = expression(rho),
      ylab = ylab,
      type = "l",
      lwd = 2,
      cex.lab=1.25
    )
    abline(h=0)

  } else {
    par(mfg = c(1, n))
    plot(
      NULL,
      xlim = lims_obj$xlim,
      ylim = lims_obj$ylim,
      xlab = "",
      ylab = "",
      xaxt="n", yaxt="n"
    )
    lines(
      x = bulk_obj$rho,
      y = omega_obj$omega,
      lty = 2,
    )
  }

  abline(h = omega_obj$w_c, lty = 2, col="grey")
  points(x = omega_obj$r_c, y = omega_obj$w_c, col = 2, pch = 16)
}

compute_potentials <- function(bulk) {
  o_obj <- compute_omega(bulk = bulk)
  o_eff_obj <- compute_omega_eff(bulk = bulk)

  return(list(w=o_obj, w_eff=o_eff_obj))
}

plot_potentials <- function(bulk, limits) {

  is_first <- is.list(limits) & length(limits) == 0

  if (is_first) { par(mfrow=c(1,2)) }

  w_objs <- compute_potentials(bulk = bulk)

  o_obj <- w_objs$w
  lims <- get_limits(
    is_first = is_first,
    omega_obj = o_obj,
    limits_obj = limits,
    n = 1
  )

  o_eff_obj <- w_objs$w_eff
  lims_eff <- get_limits(
    is_first = is_first,
    omega_obj = o_eff_obj,
    limits_obj = limits,
    n = 2
  )

  plot_panel(
    n=1,
    bulk_obj = bulk,
    omega_obj = o_obj,
    lims_obj = lims,
    is_first = is_first
  )

  plot_panel(
    n=2,
    bulk_obj = bulk,
    omega_obj = o_eff_obj,
    lims_obj = lims_eff,
    is_first = is_first
  )

  return(list(
    xlim=lims$xlim,
    ylim=lims$ylim,
    xlim_eff=lims_eff$xlim,
    ylim_eff=lims_eff$ylim))
}
