
from multiprocessing import Pool
import subprocess
import os
import glob

pool_size = 6

def bash_the_script(script_file_name: str):
    log_file_name = script_name_to_log_name(script_file_name)
    print("### bash", script_file_name, log_file_name)
    process = subprocess.run(["bash", script_file_name], check=True, stdout=subprocess.PIPE, text=True)
    output = process.stdout
    with open(log_file_name, "w") as f_out:
        f_out.write(output)

def script_name_to_log_name(script_file_name):
    s = script_file_name[:-3] + ".log"
    return s

if __name__ == '__main__':

    today = "22_10_24" # generate this
    # per line we have the script file name and the log file name
    this_weeks_refmacs_file_name = today + "_refmacs.table"

    ls = []
    # run refmac if there was not a log file for that script
    refmac_subdirs = glob.glob("refmac/*")
    for subdir in refmac_subdirs:
        subdir_scripts = glob.glob(subdir + "/refmac-*.sh")
        for script_file_name in subdir_scripts:
            log_file_name = script_name_to_log_name(script_file_name)
            if os.path.exists(log_file_name):
                pass
            else:
                # print("run that bad boy", script_file_name)
                ls.append(script_file_name)


    with Pool(pool_size) as pool:
        m = pool.map(bash_the_script, ls)

        with open(this_weeks_refmacs_file_name) as f_refmacs_table:
            for scr in ls:
                f_refmacs_table.write(scr)
                f_refmacs_table.write(" ")
                f_refmacs_table.writee(script_name_to_log_name(scr))
                f_refmacs_table.write("\n")

# this is the shell script:
#
# today=$(date +%y_%m_%d)
# this_weeks_refmacs=$today"_refmacs.table"
# echo "" > $this_weeks_refmacs
# for refmac_script in refmac/*/refmac*.sh ;
# do
#    stub="${refmac_script%.*}"
#    log=$stub.log
#    # force a rerun
#    log=$stub.logx
#    # echo $log
#    if [ -e $log ] ; then
#       :
#    else
#       echo run refmac script $refmac_script $log
#       # bash $refmac_script > $log 2>&1 &
#       # echo $refmac_script $log >> $this_weeks_refmacs
#    fi
# done
