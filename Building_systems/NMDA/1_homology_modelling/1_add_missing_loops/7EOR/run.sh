remodel=/wrk2/pirnessa/rosetta_bin_linux_2018.48.60516_bundle/main/source/bin/remodel.mpi.linuxgccrelease


mpirun -np 94 $remodel @flag_missing_loops_A
mpirun -np 94 $remodel @flag_missing_loops_B
mpirun -np 94 $remodel @flag_missing_loops_C
mpirun -np 94 $remodel @flag_missing_loops_D

