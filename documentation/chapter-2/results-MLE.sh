#!/bin/bash
## for mencia
## THIS NEEDS TO BE RUN FROM /soe/gdinolov/PDE-solvers

cd ~/PDE-solvers

echo $PWD
rhos_basis=(0.80)
sigmas_x_basis=(0.20)
sigma_y_basis=0.20
dx=250
data_size=64

data_sets=($(ls ./src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-${data_size}))
for i in ${!data_sets[@]};
do
    if [ $i -eq 0 ]; then
	echo ${PWD}/src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-${data_size}/${data_sets[$i]} > list_for_mle_results_tmp.txt
    elif [ $i -lt 50 ]; then
       	echo ${PWD}/src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-${data_size}/${data_sets[$i]} >> list_for_mle_results_tmp.txt
    fi
done

path_to_dir=/src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-0.95-n-${data_size}/
for rho_basis_current in ${rhos_basis[@]};
do
    for sigma_x_basis in ${sigmas_x_basis[@]}; 
    do
	echo $rho_basis_current $sigma_x_basis $data_size
	results_files=($(ls .${path_to_dir}TEST-dx-${dx}-analytic-deriv*rho_basis-${rho_basis_current}-sigma_x_basis-${sigma_x_basis}-sigma_y_basis-${sigma_y_basis}-mle*.csv))
	
	output_name=${PWD}${path_to_dir}dx-${dx}-analytic-deriv-rho_basis-${rho_basis_current}-sigma_x_basis-${sigma_x_basis}-sigma_y_basis-${sigma_y_basis}-mle-results-ALL.csv

	echo sigma_x,sigma_y,rho > ${output_name}

	for i in ${!results_files[@]}; do
	    echo $(cat ${results_files[$i]} | tail -n +2) >> ${output_name}
	done

    done 
done

Rscript ./src/kernel-expansion/documentation/chapter-2/results-rogers.R list_for_mle_results_tmp.txt ${data_size}

Rscript ./src/kernel-expansion/documentation/chapter-2/plot-results.R ~/PDE-solvers${path_to_dir}dx-${dx}-analytic-deriv-rho_basis-${rhos_basis[0]}-sigma_x_basis-${sigmas_x_basis[0]}-sigma_y_basis-${sigma_y_basis}-mle-results-ALL.csv ~/PDE-solvers${path_to_dir}rogers-results.csv ~/PDE-solvers${path_to_dir}classic-results.csv ~/PDE-solvers${path_to_dir}

rm list_for_mle_results_tmp.txt
