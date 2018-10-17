#!/bin/bash
## for mencia
## THIS NEEDS TO BE RUN FROM /soe/gdinolov/PDE-solvers

cd /soe/gdinolov/PDE-solvers

bazel --output_user_root=/home/gdinolov-tmp build //src/brownian-motion:generate-OCHL-data-points-for-chapter-2

rhos=(0.0 0.60 0.95)
ns=(8 16)
number_data_sets_per_iteration=50

for i in ${!rhos[@]}; 
do 
    for j in ${!ns[@]};
    do

    	for k in `seq 0 $((${number_data_sets_per_iteration}-1))`; 
    	do
	    echo $i, $j, $k, rho = ${rhos[$i]}, n = ${ns[$j]}, $(((((($i+$j))*$number_data_sets_per_iteration)) + $k))
	    seed=$(((((($i+$j))*$number_data_sets_per_iteration)) + $k))
	    echo n = ${ns[$j]}, rho = ${rhos[$i]}, seed = $seed
	    
	    path_to_dir=./src/kernel-expansion/documentation/chapter-2/data/mle-data-sets-rho-${rhos[$i]}-n-${ns[$j]}
	    if [ -d "${path_to_dir}" ]; then
		echo Directory ${path_to_dir} exists
		
	    else 
		echo Directory ${path_to_dir} DOES NOT exist
		mkdir ${path_to_dir}
		
	    fi

	    ./bazel-bin/src/brownian-motion/generate-OCHL-data-points-for-chapter-2 ${ns[$j]} ${rhos[$i]} $seed ${path_to_dir}/data-set-${k}.csv
    	
    	done
    done
done
