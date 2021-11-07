require(tidyverse)

new.bulk <- function(beta, supsat, x_v.coex, x_l.coex, rho_vap, rho_liq, r_min, n_points) {
  #' Generates a bulk object with all the properties initialised:
  bulk <- new.env()
  bulk$supsat <- supsat
  bulk$beta <- beta
  bulk$kT <- 1/beta
  bulk$x_v.coex <- x_v.coex
  bulk$x_l.coex <- x_l.coex
  bulk$rho_vap <- rho_vap
  bulk$rho_liq <- rho_liq
  bulk$r_min <- r_min
  bulk$rho <- seq(rho_vap, rho_liq, length=n_points)
  return(bulk)
}

print_bulk <- function(m_bulk) {
  writeLines(str_interp(
    paste("Bulk conditions:\n",
          "|-> kT = $[.4f]{kT}\n",
          "|-> beta = $[.4f]{beta}\n",
          "|-> Supersat = $[.4f]{supsat}\n",
          "|-> rho_v[coex] = $[.4f]{x_v}\n",
          "|-> rho_l[coex] = $[.4f]{x_l}\n",
          "|-> rho_vap = $[.4f]{rho_v}\n",
          "+-> rho_liq = $[.4f]{rho_l}\n"),
          list(
            kT=m_bulk$kT,
            beta=m_bulk$beta,
            x_v=m_bulk$x_v.coex,
            x_l=m_bulk$x_l.coex,
            supsat=m_bulk$supsat,
            rho_v=m_bulk$rho_vap,
            rho_l=m_bulk$rho_liq
          )
  ))
}

initialize <- function(beta, supsat, verbose = FALSE) {
  #' Initialises a bulk object for a given temperature (beta) and supersaturation:

  coexistence_func <- function(x)
  {#' The coexistence condition for a given temperature:
    return(eos(x = x, beta = beta))
  }

  # Find the coexistence conditions:
  xstart <- c(0.01, 0.7)
  roots <- multiroot(f = coexistence_func, start = xstart, atol = 1e-12)

  x_v.coex <- roots$root[1]
  x_l.coex <- roots$root[2]

  # Set up the supersaturated gas:
  rho_vap <- supsat * x_v.coex

  coex_liq <- function(x)
  {#' Finds the coexistent liquid for a vapor:
    return(find_liq(x, x_v = rho_vap, beta =  beta))
  }

  # Find coexisting liquid with metastabel vapor:
  root <- multiroot(f = coex_liq, start = c(1.01 * x_l.coex))
  rho_liq <- root$root

  bulk <- new.bulk(
    beta = beta,
    supsat = supsat,
    x_v.coex = x_v.coex,
    x_l.coex = x_l.coex,
    rho_vap = rho_vap,
    rho_liq = rho_liq,
    r_min = r_min,
    n_points = n_points
  )

  if (verbose) { print_bulk(m_bulk = bulk) }

  return(bulk)
}

