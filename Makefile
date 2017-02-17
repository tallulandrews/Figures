
*Simulation*.png : Figure_Simulations.R
	export R_LIBS_USER=/nfs/users/nfs_t/ta6/R-modules-3.3
	bsub -R"select[mem>25000] rusage[mem=25000]" -M25000 -o simfig2.out /software/R-3.3.0/bin/Rscript Figure_Simulations.R 

Figure1.png : Figure1S1S2.R
	export R_LIBS_USER=/nfs/users/nfs_t/ta6/R-modules-3.3
	bsub -R"select[mem>25000] rusage[mem=25000]" -M25000 -o fig1.out /software/R-3.3.0/bin/Rscript $@

Figure_Reproducibility.png : Reproducibility.R
	export R_LIBS_USER=/nfs/users/nfs_t/ta6/R-modules-3.3
	bsub -R"select[mem>25000] rusage[mem=25000]" -M25000 -o fig1.out /software/R-3.3.0/bin/Rscript $@
