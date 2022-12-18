#! /bin/bash
echo "./run_exp.sh <from> <to runs>"


for L in 0 1
do
    for cpu in 1 4 8 16 32
    do
        for grid in 128 256 512
        do
            for run in $(seq $1 1 $2)
            do
                if ((${cpu} == 32)); then
    	            echo "mpirun --use-hwthread-cpus -n ${cpu} mpi_main ${L} ${grid} 20 > res_mpi/res_${L}_${grid}_${cpu}_run_${run}"
    	            mpirun --use-hwthread-cpus -n ${cpu} mpi_main ${L} ${grid} 20 > res_mpi/res_${L}_${grid}_${cpu}_run_${run}
                else
    	            echo "mpirun -n ${cpu} mpi_main ${L} ${grid} 20 > res_mpi/res_${L}_${grid}_${cpu}_run_${run}"
    	            mpirun -n ${cpu} mpi_main ${L} ${grid} 20 > res_mpi/res_${L}_${grid}_${cpu}_run_${run}
                fi
            done
        done
    done
done
