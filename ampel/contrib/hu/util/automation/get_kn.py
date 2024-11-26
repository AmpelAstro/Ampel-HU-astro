#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/util/automation/AmpelHealpix.py
# License:             BSD-3-Clause
# Author:              andrea ernst
# Date:                20.03.2024
# Last Modified Date:  20.03.2024
# Last Modified By:    ernstand@physik.hu-berlin.de

import sys

import ligo
from ligo.gracedb.rest import GraceDb

from ampel.contrib.hu.util.automation.jobfileCaller import executeJobfile
from ampel.contrib.hu.util.automation.jobfileWriter import writeJobfilesFromDict

jobfile = ""


EXEC_DIR = "./"  #  "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/O4"
CONF_FILE = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/ampel_conf.yaml"
VAULT_FILE = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/vault.yaml"
JOBFILE_TEMP = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/examples/remote/ligo_automated_template.yml"


########### get event info

# gw_name = "GW240109a"
gw_names = sys.argv[1:]
print(type(gw_names))

print(gw_names)
for k, gw_name in enumerate(gw_names):
    # gw_name = sys.argv[1]
    print("#########################################################")
    print(f"############### - MAP {gw_name} - {k+1}/{len(gw_names)} - ###############")
    print("#########################################################")
    event_name = "S" + (gw_name.replace("GW", "").replace("S", ""))

    client = GraceDb()

    try:
        response = client.files(event_name)
    except ligo.gracedb.exceptions.HTTPError as e:
        print("GraceDB API request exited with ", e)
        sys.exit()

    # if response.status_code != 201:
    #    print("GraceDB API request exited with response code ", response.status_code)
    #    sys.exit()
    file_contents = response.json()

    fits_files = [key for key in file_contents if ".fits.gz" in key]

    # print(fits_files)

    map_url = ""
    map_file = "LALInference.fits.gz"
    if map_file not in fits_files:
        map_file = "bayestar.fits.gz"
    if map_file not in fits_files:
        print("No working map file found")
        raise ValueError(f"No map file found for superevent {event_name}")
    else:
        map_url = file_contents[map_file]

    print()
    print("+++++++++ \tFOUND EVENT ON GRACEDB \t+++++++++ \n")
    # print(map_url)

    ########### write the jobfile

    print(f"+++++++++ \t WRITING JOBFILE\t{event_name}\t +++++++++ \n")

    jobfile_list_dict = {
        gw_name: {
            "map_url": map_url,
            # "map_name": map_name,
        }
    }

    jobfile = writeJobfilesFromDict(
        JOBFILE_TEMP, jobfile_list_dict, commonName="automated", saveDir=EXEC_DIR
    )[0]

    ########### execute the jobfile
    print(f"+++++++++\tEXECUTE JOBFILE \t{jobfile}\t +++++++++ \n")

    job_call = ""
    job_call += (
        f"ampel job --schema {jobfile} --config {CONF_FILE} --secrets {VAULT_FILE}"
    )
    # print(job_call)
    executeJobfile(job_call, EXEC_DIR, max_retries=2)
