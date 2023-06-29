# %%
# notebook to execute a list of jobfiles without supervision for "stream not found" or bad exit codes
# repeats same file 5 times then moves on

import os
import subprocess
import shlex
import numpy as np
from select import select
from subprocess import Popen

import astropy.time as atime


# directory where ampel job should be run (NEEDS .tmp FOLDER)
execute_directory="/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/ligo-healpix/test2"

# directory where jobfiles are found
jobfile_save_dir="/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/calibrateKilonovaEval_jobfiles/"

#jobfiles_to_execute_file = os.listdir(jobfile_save_dir)
# list of names of jobfiles to run
jobfile_names = os.listdir(jobfile_save_dir)
#config file to use
conf_file = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/ampel_conf.yaml"

#vault file to use
vault_file = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/vault.yaml"

print(jobfile_names)


def execute_jobfile(job_call, execute_directory):
    command = job_call
    command_list = shlex.split(command)

    #print(command_list)

    stream_not_found_err = "Stream not found"

    os.chdir(execute_directory)

    timeout = 0.1
    retry = True
    max_retries = 5
    retry_ind = 0
    while retry:
        print("STARTING SUBPROCESS, try #", retry_ind)
        try:
            process = Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except NameError:
            print("Error")
        retry = False

        
        while process:
            #output, error = process.communicate()
            #print(output, error)

            if process.poll() is not None:
                if process.returncode != 0:
                    print("BAD EXIT")
                    print(process.stderr.read().decode("utf-8"), end="")
                    retry = True
                    break
                print(process.stdout.read().decode("utf-8"), end="")
                process.stdout.close()
                break
            
            rlist = select([process.stdout], [], [], timeout)[0]

            for f in rlist:
                line = f.readline().decode("utf-8")
                print(line, end="")
                if line.find(stream_not_found_err) != -1:
                    print("FOUND", stream_not_found_err)
                    retry = True
                    process.kill()
                #elif line.find("Error") != -1:
                    #print("Found SyntaxError")
                    #retry = True
                    #process.kill()
                #print(f"{f.readline()}")

        

        retry_ind += 1
        if retry_ind >= max_retries:
            print("MAX RETRIES REACHED, STOPPING")
            break

        print("RETRY?", retry)


if __name__=="__main__":
    for jobfile_name in jobfile_names:

        jobfile = os.path.join(jobfile_save_dir, jobfile_name)

        print(f"################################################ \n################################################ \n\n\t\tNEXT JOBFILE\n\t{jobfile_name}\n\n################################################ \n################################################ \n")

        job_call = "" #f"cd {execute_directory};"
        job_call += f"ampel job --schema {jobfile} --config {conf_file} --secrets {vault_file}"

        execute_jobfile(job_call, execute_directory)


