#!/bin/bash
## for mencia
## THIS NEEDS TO BE RUN FROM /soe/gdinolov/PDE-solvers

cd /soe/gdinolov/PDE-solvers

bazel --output_user_root=/home/gdinolov-tmp build //src/mle-estimator:2d-mle-finite-element-data

echo $PWD
number_threads=8
rhos_basis=(0.80)
sigmas_basis=(0.20)

data_sets=($(ls ./src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-8))
path_to_dir=./src/kernel-expansion/documentation/chapter-2/results/mle-results-rho-0.95-n-8/
for i in ${!data_sets[@]};
do
    if [ $i -eq 0 ]; then
	echo ${PWD}/src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-8/${data_sets[$i]} > list_for_mle_tmp.txt
    elif [ $i -lt 100 ]; then
     	echo ${PWD}/src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-0.95-n-8/${data_sets[$i]} >> list_for_mle_tmp.txt
    fi
done

# echo ${data_sets[1]} >> list_for_mle_tmp.txt
# echo ${data_sets[2]} >> list_for_mle_tmp.txt

for rho_basis_current in ${rhos_basis[@]};
do
    for sigma_basis in ${sigmas_basis[@]}; 
    do
	echo $rho_basis_current $sigma_basis
	if [ ! -d "${path_to_dir}" ]; then
	    mkdir ${path_to_dir}
	fi
	
	./bazel-bin/src/mle-estimator/2d-mle-finite-element-data list_for_mle_tmp.txt ${path_to_dir} 0.0001 1.0 1.0 0.0 $rho_basis_current $sigma_basis $sigma_basis TEST-dx-600-analytic-deriv-rho_basis-${rho_basis_current}-sigma_x_basis-${sigma_basis}-sigma_y_basis-${sigma_basis}- -SUFFIX-TEST ${number_threads}
    done 
done