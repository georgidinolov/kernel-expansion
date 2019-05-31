#!/bin/bash
## for mencia
## THIS NEEDS TO BE RUN FROM /soe/gdinolov/PDE-solvers

cd /soe/gdinolov/PDE-solvers

bazel --output_user_root=/home/gdinolov-tmp build //src/mle-estimator:2d-mle-finite-element-data

echo $PWD
number_threads=16
rhos_basis=(0.70)
sigmas_y_basis=(0.05)
sigmas_x_basis=(0.20)
dx=300
data_sizes=(16)
rho_data=0.95


for data_size in ${data_sizes[@]};
do 

    data_sets=($(ls ./src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-${rho_data}-n-${data_size}/data-set-*.csv))
    path_to_dir=./src/kernel-expansion/documentation/chapter-3/results/mle-results-rho-${rho_data}-n-${data_size}/

    for i in ${!data_sets[@]};
    do
	if [ $i -eq 0 ]; then
	    echo ${data_sets[$i]} > list_for_mle_tmp.txt
	elif [ $i -lt 50 ]; then
       	    echo ${data_sets[$i]} >> list_for_mle_tmp.txt
	fi
    done

    for rho_basis_current in ${rhos_basis[@]};
    do
 	for j in ${!sigmas_y_basis[@]};
	do
	    if [ ! -d "${path_to_dir}" ]; then
    		mkdir ${path_to_dir}
    	    fi

	    prefix=TEST-MATCHING-dx-${dx}-analytic-deriv-rho_basis-${rho_basis_current}-sigma_x_basis-${sigmas_x_basis[$j]}-sigma_y_basis-${sigmas_y_basis[$j]}-
	    
	    echo path_to_dir = ${path_to_dir} 
	    echo rho_data = ${rho_data} 
	    echo rho_basis_current = $rho_basis_current 
	    echo sigma_x_basis[$j] = ${sigmas_x_basis[$j]} 
	    echo sigma_y_basis[$j] = ${sigmas_y_basis[$j]} 
	    echo prefix = ${prefix} 
	    echo numbers_thread = ${number_threads} 
	    echo dx = ${dx}

	    ./bazel-bin/src/mle-estimator/2d-mle-finite-element-data list_for_mle_tmp.txt ${path_to_dir} 0.0001 1.0 1.0 ${rho_data} $rho_basis_current ${sigmas_x_basis[$j]} ${sigmas_y_basis[$j]} ${prefix} -SUFFIX-TEST ${number_threads} ${dx}
	done 
    done
done 


