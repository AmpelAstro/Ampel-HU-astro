jobfile_event_file=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/scripts/jobfileWriter/event_lists/O4_events.yaml
jobfile_template_file=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/scripts/jobfileWriter/templates/ligo_healpix_template.yml
jobfile_writer=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/scripts/jobfileWriter/jobfileWriter.py

jobfile_save_dir=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/O4_jobfiles
jobfile_prefix="O4"

#eval "$(conda shell.bash hook)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ampel-hu
python3 $jobfile_writer $jobfile_template_file $jobfile_event_file $jobfile_prefix $jobfile_save_dir

