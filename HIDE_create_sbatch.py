import os,sys
import argparse
import file_verification as ver
import numpy as np
###################################################################
# Check the python version and import configparser
###################################################################
if sys.version_info[0]==2:
    import ConfigParser
    config = ConfigParser.RawConfigParser()
elif sys.version_info[0]==3:
    import configparser
    config = configparser.ConfigParser()
###############################################################################
#You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument( '-s','--Hstart'       , action = 'store', dest =  'Hstart'      , default = 1                                      , help = 'First position index to be submitted. Ps.: each position index doenst correspond to a horn but a specific pointing on the sky. This is because one horn could pinpoint more than one position on the sky')
parser.add_argument( '-e','--Hend'         , action = 'store', dest =  'Hend'        , default = 140                                    , help = 'Last position index to be submitted')
parser.add_argument( '-q','--SD_queue'     , action = 'store', dest =  'SD_queue'    , default = 'sequana_cpu_dev'                      , help = 'SDumont row')
parser.add_argument( '-po','--path_dir_out', action = 'store', dest =  'path_dir_out', default = "/scratch/bingo/thiago.pena"           , help = 'directory where we will (make a new directory  sbatch_files and) save new sbatch files')
parser.add_argument( '-ph','--path_hide'   , action = 'store', dest =  'path_hide'   , default = "/scratch/bingo/thiago.pena/extra/beam_test/hide", help = 'path to HIDE directory')

parser.add_argument( '-pho','--path_hide_output'   , action = 'store', dest =  'path_hide_out'   , default = "/scratch/bingo/thiago.pena/extra/beam_test/240412/seedtest2700", help = 'path to HIDE output directory')


###############################################################################
#Variables
###############################################################################
arguments = parser.parse_args()
Hstart       = int(arguments.Hstart)
Hend         = int(arguments.Hend)
SD_queue     = str(arguments.SD_queue)
path_dir_out = str(arguments.path_dir_out)
path_hide    = str(arguments.path_hide)

path_hide_out = str(arguments.path_hide_out)
###############################################################################
###############################################################################
#path_dir = "/scratch/bingo/thiago.pena" #os.getcwd()
ver.file_verification(path_dir_out," ","sbatch_files")
path = os.path.join(path_dir_out,"sbatch_files/")
Hid  = np.arange(Hstart-1,Hend,1)
print("Saving \".sh\" files for horn between {} and {} position index".format(Hstart,Hend))
print("Directory path: {}".format(path))
for i in Hid:
	f = open("".join((path,"sbatch_HIDE_position_",str(i+1),".sh")),"w+")
	f.write("#!/bin/bash -f\n\n")
	f.write("".join(("#SBATCH --job-name=ID_",str(i+1),"\n")))
	f.write("".join(("#SBATCH -N1 -n1\n")))
	f.write("#SBATCH --mem=10G\n")
	f.write("#SBATCH -t 00:20:00\n")
	f.write("#SBATCH -p {}\n".format(SD_queue))
	f.write("#SBATCH -o {}/horn_{}.out\n".format(path_hide_out,i+1))
	f.write("#SBATCH -e {}/horn_{}.err\n\n".format(path_hide_out,i+1))
	f.write("module load openmpi/gnu/4.0.1_gcc-7.4\n")
	f.write("source /scratch/bingo/thiago.pena/extra/venv/bin/activate\n\n")
	f.write("echo Time is `date`\n")
	f.write("echo Running on host `hostname`\n")
	f.write("echo Directory is `pwd`\n")
	f.write("echo This jobs runs on the following processors:\n")
	f.write("echo $SLURM_JOB_NODELIST\n\n")
	f.write("mpirun python3 {}/run_hide.py \"bingo.py\" \"bingo_horn\" {} {}\n\n".format(path_hide,i,i+1))
	f.write("cp {}/hide/config/bingo_horn_{}.py {}/maps\n".format(path_hide, i, path_hide_out))
	f.write("echo Time is `date`\n")
	f.close()
print("Files saved.")
