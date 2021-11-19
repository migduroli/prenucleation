source("utils/params.R")
source("utils/thermo.R")
source("utils/bulk.R")
source("utils/graphs.R")



pdf(file = "figures/potentials.Srange.pdf", width = 6, height = 6)
supsat_values <- seq(1.75, 3.25, length = 10)
limits <- list()
for (s in supsat_values)
{
  m_bulk <- initialize(beta = b_beta, supsat = s, verbose = F)
  limits <- plot_potentials(bulk = m_bulk, limits = limits)
  # limits <- plot_probability(bulk = m_bulk, limits = limits)
}
dev.off()

pdf(file = "figures/probabilities.S.range.pdf", width = 6, height = 6)
limits <- list()
for (s in supsat_values)
{
  m_bulk <- initialize(beta = b_beta, supsat = s, verbose = F)
  limits <- plot_probability(bulk = m_bulk, limits = limits, arrow = TRUE, critical = TRUE)
}
dev.off()

supsat_values <- seq(1.01, 1.75, length = 10)
pdf(file = "figures/probabilities.S.low.pdf", width = 6, height = 6)
limits <- list()
for (s in supsat_values)
{
  m_bulk <- initialize(beta = b_beta, supsat = s, verbose = F)
  limits <- plot_probability(bulk = m_bulk, limits = limits, arrow = FALSE, critical = TRUE)
}
dev.off()


supsat_values <- seq(0.75, 1.0, length = 10)
pdf(file = "figures/probabilities.S.subsat.pdf", width = 6, height = 6)
limits <- list()
for (s in supsat_values)
{
  m_bulk <- initialize(beta = b_beta, supsat = s, verbose = F)
  limits <- plot_probability(bulk = m_bulk, limits = limits, arrow = FALSE, critical = FALSE)
}
dev.off()

