#!/bin/bash
## for mencia
## THIS NEEDS TO BE RUN FROM /soe/gdinolov/PDE-solvers

cd /soe/gdinolov/PDE-solvers

echo $PWD
rhos_basis=(0.70)
sigmas_x_basis=(0.20)
sigma_y_basis=0.20

# data_sets=($(ls ./src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-8))
path_to_dir=/src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-0.95-n-8/

for rho_basis_current in ${rhos_basis[@]};
do
    for sigma_x_basis in ${sigmas_x_basis[@]}; 
    do
	echo $rho_basis_current $sigma_x_basis
	results_files=($(ls .${path_to_dir}TEST-dx-analytic-deriv*rho_basis-${rho_basis_current}-sigma_x_basis-${sigma_x_basis}-sigma_y_basis-${sigma_y_basis}-mle*.csv))
	
	output_name=${PWD}${path_to_dir}dx-analytic-deriv-rho_basis-${rho_basis_current}-sigma_x_basis-${sigma_x_basis}-sigma_y_basis-${sigma_y_basis}-mle-results-ALL.csv

	echo sigma_x,sigma_y,rho > ${output_name}

	for i in ${!results_files[@]}; do
	    echo $(cat ${results_files[$i]} | tail -n +2) >> ${output_name}
	done

    done 
done

Rscript ./src/kernel-expansion/documentation/chapter-2/results-rogers.R list_for_mle_tmp.txt

Rscript ./src/kernel-expansion/documentation/chapter-2/plot-results.R /soe/gdinolov/PDE-solvers${path_to_dir}dx-analytic-deriv-rho_basis-0.70-sigma_x_basis-0.20-sigma_y_basis-0.20-mle-results-ALL.csv /soe/gdinolov/PDE-solvers${path_to_dir}rogers-results.csv /soe/gdinolov/PDE-solvers${path_to_dir}classic-results.csv /soe/gdinolov/PDE-solvers${path_to_dir}