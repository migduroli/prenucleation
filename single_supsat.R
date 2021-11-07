source("utils/params.R")
source("utils/thermo.R")
source("utils/bulk.R")
source("utils/graphs.R")


m_bulk <- initialize(beta = b_beta, supsat = b_supsat)

## The following assignments are not needed
## Bulk vars:
# supsat <- m_bulk$supsat
# rho_vap <- m_bulk$rho_vap
# rho_liq <- m_bulk$rho_liq
# r_min <- m_bulk$r_min
#
## Density grid:
# rho <- m_bulk$rho

# The set of radii associated to the range of cluster densities under study
plot_radius(m_bulk)

# Plot omega:
plot_omega(m_bulk)

# Plot of g:
plot_g(m_bulk)


# Plot omega vs omega_eff:
plot_potentials(m_bulk, limits = list())

