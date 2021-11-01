initialize <- function()
{
  coexistence_func <- function(x)
  {#' The coexistence condition for a given temperature:
    return(eos(x = x, beta = b_beta))
  }

  # Find the coexistence conditions:
  xstart <- c(0.01, 0.7)
  roots <- multiroot(f = coexistence_func, start = xstart, atol = 1e-12)

  x_v.coex <- roots$root[1]
  x_l.coex <- roots$root[2]

  print(paste("Conditions: kT = ", kT, " beta = ", b_beta))
  print(paste("rho_v[coex] = ",x_v.coex, " rho_l[coex] = ", x_l.coex))

  # Set up the supersaturated gas:
  rho_vap <- supsat * x_v.coex

  coex_liq <- function(x)
  {#' Finds the coexistent liquid for a vapor:
    return(find_liq(x, x_v = rho_vap, beta =  b_beta))
  }

  # Find coexisting liquid with metastabel vapor:
  root <- multiroot(f = coex_liq, start = c(1.01 * x_l.coex))
  rho_liq <- root$root

  bulk <- new.env()

  bulk$supsat <- supsat
  bulk$beta <- b_beta
  bulk$x_v.coex <- x_v.coex
  bulk$x_l.coex <- x_l.coex
  bulk$rho_vap <- rho_vap
  bulk$rho_liq <- rho_liq
  bulk$r_min <- r_min
  bulk$rho <- seq(rho_vap, rho_liq, length=n_points)

  return(bulk)
}

