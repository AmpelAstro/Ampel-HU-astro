
# directory where ampel job should be run (NEEDS .tmp FOLDER)
execute_directory=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/ligo-healpix/test2

# directory where jobfiles are found
jobfile_save_dir=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/calibrateKilonovaEval_jobfiles/

# list of names of jobfiles to run
jobfile_names=("ligo_healpix_S200115j.yaml" "ligo_healpix_S190426c.yaml")

# ampel job call without file
call_template="ampel job --schema <<JOBFILE>> --config $MASTERARBEIT/ampel/Ampel-HU-astro/ampel_conf.yaml --secrets $MASTERARBEIT/ampel/Ampel-HU-astro/vault.yaml"

cd $execute_directory
source activate ampel-hu

for file in "${jobfile_names[@]}"; do
	echo $file
	${call_template/<<JOBFILE>>/"$jobfile_save_dir$file"}
done


