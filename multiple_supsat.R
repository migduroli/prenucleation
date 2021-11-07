source("utils/params.R")
source("utils/thermo.R")
source("utils/bulk.R")
source("utils/graphs.R")

supsat_values <- seq(1.75, 3.25, length = 10)

limits <- list()
for (s in supsat_values)
{
  m_bulk <- initialize(beta = b_beta, supsat = s, verbose = T)
  limits <- plot_potentials(bulk = m_bulk, limits = limits)
}
