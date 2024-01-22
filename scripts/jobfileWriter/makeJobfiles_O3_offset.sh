jobfile_event_file=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/O3_offset_jobfiles/O3_offset_events.yml

jobfile_template_file=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/scripts/jobfileWriter/ligo_healpix_template.yml

jobfile_writer=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/scripts/jobfileWriter/jobfileWriter.py

jobfile_save_dir=/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/O3_offset_jobfiles
jobfile_prefix="O3_offset"


source ~/miniconda3/etc/profile.d/conda.sh
conda activate ampel-hu
python3 $jobfile_writer $jobfile_template_file $jobfile_event_file $jobfile_prefix $jobfile_save_dir

