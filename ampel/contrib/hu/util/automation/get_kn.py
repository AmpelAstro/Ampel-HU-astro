import os
import sys
from ampel.contrib.hu.util.automation.jobfileWriter import writeJobfilesFromDict
from ampel.contrib.hu.util.automation.jobfileCaller import executeJobfile

from ligo.gracedb.rest import GraceDb
import ligo

jobfile = ""


EXEC_DIR = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/remote"
CONF_FILE = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/ampel_conf.yaml"
VAULT_FILE = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/vault.yaml"
JOBFILE_TEMP = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/remote/ligo_automated_template.yml"


########### get event info

#gw_name = "GW240109a"
gw_name = sys.argv[1]
gw_name = gw_name.replace("GW", "")
gw_name = gw_name.replace("S", "")
event_name = "S" + gw_name

client = GraceDb()

try: 
    response =  client.files(event_name)
except ligo.gracedb.exceptions.HTTPError as e:
    print("GraceDB API request exited with ", e)
    sys.exit()

#if response.status_code != 201:
#    print("GraceDB API request exited with response code ", response.status_code)
#    sys.exit()
file_contents = response.json()


fits_files = [key for key in file_contents.keys() if ".fits.gz" in key]

#print(fits_files)

map_url = ""
map_file = 'LALInference.fits.gz'
if not map_file in fits_files:
    map_file = 'bayestar.fits.gz'
if not map_file in fits_files:
    print("No working map file found")
    raise ValueError(f"No map file found for superevent {event_name}")
else: 
    map_url = file_contents[map_file]

print()
print(f"+++++++++ \tFOUND EVENT ON GRACEDB \t+++++++++ \n")
#print(map_url)

########### write the jobfile

print(f"+++++++++ \t WRITING JOBFILE\t{event_name}\t +++++++++ \n")

jobfile_list_dict = {gw_name: {"map_url": map_url,
                                #"map_name": map_name,
                                }}

jobfile = writeJobfilesFromDict(JOBFILE_TEMP, jobfile_list_dict, commonName = "automated", saveDir = EXEC_DIR)[0]




########### execute the jobfile
print(f"+++++++++\tEXECUTE JOBFILE \t{jobfile}\t +++++++++ \n")


job_call = ""
job_call += f"ampel job --schema {jobfile} --config {CONF_FILE} --secrets {VAULT_FILE}"
#print(job_call)
executeJobfile(job_call, EXEC_DIR, max_retries=2)