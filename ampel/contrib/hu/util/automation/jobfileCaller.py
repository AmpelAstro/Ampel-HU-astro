#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/util/automation/AmpelHealpix.py
# License:             BSD-3-Clause
# Author:              andrea ernst
# Date:                20.03.2024
# Last Modified Date:  20.03.2024
# Last Modified By:    ernstand@physik.hu-berlin.de

import os
import shlex
import subprocess
from select import select
from subprocess import Popen


def executeJobfile(job_call, execute_directory, should_retry=True, max_retries=5):
    command = job_call
    command_list = shlex.split(command)

    # print(command)

    stream_not_found_err = "Stream not found"

    os.chdir(execute_directory)
    # print(os.getcwd())
    # print(execute_directory)

    timeout = 0.1
    retry = True
    retry_ind = 0
    while retry:
        print("STARTING SUBPROCESS, try #", retry_ind)
        try:
            process = Popen(
                command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        except NameError:
            print("Error")
        retry = False

        while process:
            # output, error = process.communicate()
            # print(output, error)

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
                # elif line.find("Error") != -1:
                # print("Found SyntaxError")
                # retry = True
                # process.kill()
                # print(f"{f.readline()}")

        retry_ind += 1
        if retry_ind >= max_retries or not should_retry:
            print("MAX RETRIES REACHED, STOPPING")
            break

        print("RETRY?", retry)
