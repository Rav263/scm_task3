#! /bin/bash
echo "./run_exp.sh <from> <to runs>"


for L in 0 1
do
    for cpu in 1 2 4 8
    do
        for grid in 128 256 512
        do
            for run in $(seq $1 1 $2)
            do
                if ((${cpu} == 8)); then
    	            echo "mpirun --bind-to core:overload-allowed --cpus-per-proc 4 -np ${cpu} omp_main ${L} ${grid} 20 > res_omp/res_${L}_${grid}_${cpu}_run_${run}"
    	            mpirun --bind-to core:overload-allowed --cpus-per-proc 4 -np ${cpu} omp_main ${L} ${grid} 20 > res_omp/res_${L}_${grid}_${cpu}_run_${run}
                else
    	            echo "mpirun --cpus-per-proc 4 -np ${cpu} omp_main ${L} ${grid} 20 > res_omp/res_${L}_${grid}_${cpu}_run_${run}"
    	            mpirun --cpus-per-proc 4 -np ${cpu} omp_main ${L} ${grid} 20 > res_omp/res_${L}_${grid}_${cpu}_run_${run}
                fi
            done
        done
    done
done
