# Author: ernstand@physik.hu-berlin.de
# Last modified date: 22.05.2023
# Last modified by: AE 

import os
import sys

import astropy.time as atime
import yaml


def writeJobfilesFromDict(jobfile_template, jobfile_list_dict, commonName = "jobfile", saveDir = "."):  
    """
    Automatically generate jobfiles based on a template, given a dictionary of values to replace in template. Template has to have <<VALUE>> placeholders:

    config:
      pvalue_limit: 0.9 
      chunk_size: 2000
      map_name: <<map_name>>
      map_url: <<map_url>>
      map_dir: './tmp'
    
    """   

    jobcall_template = "ampel job --schema <<jobfile>> --config $MASTERARBEIT/ampel/Ampel-HU-astro/ampel_conf.yaml --secrets $MASTERARBEIT/ampel/Ampel-HU-astro/vault.yaml"
    
    jobfile_list = []


    with open(jobfile_template) as template:
        for map, replaceDict in jobfile_list_dict.items():

            newJobfileName = commonName + "_" + map + ".yaml" 

            if replaceDict.get("map_name"):
                map_name = replaceDict.get("map_name")
                replaceDict["map_token"] = replaceDict["map_name"] + "_token"
            else:
                map_name = map
                if replaceDict["map_url"].find(",") != -1:
                    map_name += replaceDict["map_url"][-10:]
                else:
                    map_name += replaceDict["map_url"][-8:]
                replaceDict["map_name"] = map_name
                replaceDict["map_token"] = replaceDict["map_name"] + "_token"
            print(replaceDict)
            if replaceDict.get("trigger_jd"):
                tmp_time = atime.Time(replaceDict["trigger_jd"], format="jd")
                tmp_time.format = "iso"
                date_str = "'" + tmp_time.to_string() + "'"
                print(date_str, replaceDict["trigger_jd"])
                replaceDict["date_str"] = date_str


            print(f"Creating jobfile {newJobfileName} ...")

            with open(os.path.join(saveDir, newJobfileName), "w") as response:

                for old_line in template:
                    # go line by line and replace placeholder values 
                    newline = old_line

                    for placeHolder, value in replaceDict.items():
                        newline = newline.replace("<<" + placeHolder + ">>", str(value))

                    response.write(newline)   
            jobfile_list.append(newJobfileName)
            print(jobcall_template.replace("<<jobfile>>",  os.path.abspath(response.name)), "\n") 
            response.close()
            template.seek(0, 0)
    template.close()
    return jobfile_list
    #print(jobfile_list)
    #print("Done!")


def writeJobfilesFromYaml(jobfile_template, jobfile_list_file, commonName = "jobfile", saveDir = "."):
    """ Reads in .yaml file as dict of jobfiles to automatically generate. .yaml file must have following schema:
    
    EVENTNAME1:
      name_in_place_holder1: value1
      name_in_place_holder2: value2

    EVENTNAME2:
      name_in_place_holder1: value1
      name_in_place_holder2: value2
    
    """

    with open(jobfile_list_file) as f:
        jobfile_list_dict = yaml.safe_load(f)
    writeJobfilesFromDict(jobfile_template, jobfile_list_dict, commonName = commonName, saveDir = saveDir)

def main():
    """
    Inputs are, in order:
    - template of jobfile to copy off
    - yaml file with events and values to overwrite
    - common naming scheme for jobfiles (optional)
    - save directory (optional)
    """
    writeJobfilesFromYaml(*sys.argv[1:])

if __name__ == "__main__":
    main()
