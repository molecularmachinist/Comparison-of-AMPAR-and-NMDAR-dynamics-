# Read https://www.rosettacommons.org/demos/latest/tutorials/loop_modeling/loop_modeling for more information.
# The .remodel file is created by editing file that can be optained by uncommenting the following line
#/home/fabiolol/rosetta_bin_linux_2018.48.60516_bundle/tools/remodel/getBluePrintFromCoords.pl -pdbfile ../0_initial_structure/5weo_cleaned.pdb -chain A > 5weo_chain_temp.remodel


mpirun -np 94 /home/fabiolol/rosetta_bin_linux_2018.48.60516_bundle/main/source/bin/remodel.mpi.linuxgccrelease @flag_missing_loops_A

mpirun -np 94 /home/fabiolol/rosetta_bin_linux_2018.48.60516_bundle/main/source/bin/remodel.mpi.linuxgccrelease @flag_missing_loops_B

mpirun -np 94 /home/fabiolol/rosetta_bin_linux_2018.48.60516_bundle/main/source/bin/remodel.mpi.linuxgccrelease @flag_missing_loops_C

mpirun -np 94 /home/fabiolol/rosetta_bin_linux_2018.48.60516_bundle/main/source/bin/remodel.mpi.linuxgccrelease @flag_missing_loops_D

