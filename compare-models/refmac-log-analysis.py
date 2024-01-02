import glob
import os

this_weeks_refmacs_table = "this_weeks_refmacs.table"
this_weeks_refmacs_table = "15-11-refmacs.table"


def check_for_error(refmac_log_file_name_in):
    refmac_log_file_name = str(refmac_log_file_name_in)
    print('checking "{}"'.format(refmac_log_file_name)) # may not be needed
    with open(refmac_log_file_name, "r") as refmac_log_file:
        # print("analysing log", refmac_log_file_name)
        for line in refmac_log_file.readlines():
            if "Error termination." in line:
                print("This log had a error termination:", refmac_log_file_name)
            if "Error: Problem with ligand name" in line:
                print("This log had a error due to ligand name:", refmac_log_file_name)


with open(this_weeks_refmacs_table) as f:
    lines = f.readlines()
    for line in lines:
        parts = line.split()
        if parts:
            shel_script_file_name = parts[0]
            refmac_log_file_name  = parts[1]
            check_for_error(refmac_log_file_name)
