#! /bin/bash

for grid in 128 256 512
do
	tt=15
	if ((grid == 512))
	then
		tt=30
	fi
	for l in 0 1
  	do
        	for run in 1
        	do
            		for cpu in 1 4 8 16
            		do
                		bsub -x -n ${cpu} -R span[hosts=1] -W 00:${tt} -o polus_mpi_re_exc/res_${grid}_${cpu}_${l}_${run} -e all_errors_re_exc/err_${grid}_${cpu}_${l}_${run} "mpiexec ./reorder_mpi_main ${l} ${grid} 20"
            		done
	    		for cpu in 32
			do
                		bsub -x -n ${cpu} -R span[ptile=16] -W 00:${tt} -o polus_mpi_re_exc/res_${grid}_${cpu}_${l}_${run} -e all_errors_re_exc/err_${grid}_${cpu}_${l}_${run} "mpiexec ./reorder_mpi_main ${l} ${grid} 20"
			done
            		for cpu in 1 2 4
            		do
                		bsub -x -n ${cpu} -R span[hosts=1] -R "affinity[core(4)]" -W 00:${tt} -o polus_omp_re_exc/res_${grid}_${cpu}_${l}_${run} -e all_errors_re_exc/err_${grid}_${cpu}_${l}_${run} "OMP_NUM_THREADS=4; mpiexec ./reorder_omp_main ${l} ${grid} 20"
            		done
            		for cpu in 8
			do
                		bsub -x -n ${cpu} -R span[ptile=4] -R "affinity[core(4)]" -W 00:${tt} -o polus_omp_re_exc/res_${grid}_${cpu}_${l}_${run} -e all_errors_re_exc/err_${grid}_${cpu}_${l}_${run} "OMP_NUM_THREADS=4; mpiexec ./reorder_omp_main ${l} ${grid} 20"
			done	
        	done
	done
done
