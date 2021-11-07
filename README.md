# Prenucleation 💥

Welcome to prenucleation! 

This is nothing but a scientific repository 
which aims at helping with the research on the
occurrence and effects of intermediate (metastable) phases
in nucleation.

The code is written in `R` as this is the most convenient
language for the different parties of the project.

So far the project consists in:

- [bulk](utils/bulk.R): Module where the initialisation of the bulk phase is done
- [thermo](utils/thermo.R): Module with the main functions about the bulk and cluster thermodynamic properties
- [params](utils/params.R): Module where the configuration parameters (e.g. `kT` or `hard_sphere`) are set
- [main](single_supsat.R): This is the entry point to generate the figures required

