import os
import subprocess

todays_refmacs_table = "15-11-refmacs.table"

def log_name_to_code(log):
    parts_1 = log.split(".")
    pre = parts_1[0]
    parts_2 = pre.split("-")
    code = parts_2[1]
    return code

if os.path.exists("rob-proc"):
    pass
else:
    os.mkdir("rob-proc")

if os.path.exists("correlation-tables"):
    pass
else:
    os.mkdir("correlation-tables")

with open(todays_refmacs_table) as f:
    lines = f.readlines()
    for line in lines:
        parts = line.split()
        log = parts[1]
        accession_code = log_name_to_code(log)
        print(log, accession_code)

        rob_proc_log_file_name = os.path.join("rob-proc", accession_code + ".log")
        args = ["python3", "robsrepo/Ligand-Tools/compare-models/compare_models_for_ligands_of_the_week.py",
                '--code', accession_code, '--ligand-analysis-base', "."]
        comp_proc = subprocess.run(args, capture_output=True)
        with open(rob_proc_log_file_name, "bw") as f:
            f.write(comp_proc.stdout)
        if os.path.exists("results.out"):
            new_results_out_name = os.path.join("correlation-tables", "results-" + accession_code + ".table")
            os.rename("results.out", new_results_out_name)
