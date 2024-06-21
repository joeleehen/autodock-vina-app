idev:
	idev -p development -m 120 -n 32 -N 1 -A TACC-SCI

pull:
	apptainer pull docker://wjallen/autodock_vina:1.2.3
	
test:
	MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=0 ibrun -np 32 apptainer exec autodock_vina_1.2.3.sif python3 autodock.py -r fullHMGB1_human.pdb --center="-19.4,4.4,14" -s "20,20,20" -f vina -d rigid -l /scratch/projects/docking/Test-set-compressed -n 15 -x 'EMPTY'

test-big:
	MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=0 ibrun -np 32 apptainer exec autodock_vina_1.2.3.sif python3 autodock.py -r fullHMGB1_human.pdb --center="-19.4,4.4,14" -s "20,20,20" -f vina -d rigid -l /scratch/projects/docking/ZINC-in-trials-compressed -n 50 -x 'EMPTY'

clean:
	rm -rf configs/ output/ results* merged_results.txt autodock.log results.tar.gz


