#!/bin/bash
## for mencia
## THIS NEEDS TO BE RUN FROM /soe/gdinolov/PDE-solvers/

rho_data=0.95

cd ~/PDE-solvers
./src/finite-element-igraph/likelihood-profiles.sh 0.60 0.10 0.10 16 16 ${rho_data}

cd ~/PDE-solvers/src/finite-element-igraph/
save_file_prefix=~/PDE-solvers/src/kernel-expansion/documentation/chapter-3/figures/matched-rho-${rho_data}-
Rscript plot-likelihood-profiles-with-matching.R ${save_file_prefix}
