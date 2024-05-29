##############################
# Adaptaded from Gabriel Hoerning's code to execute a sequential file submissions on SDumont cluster
# Last update:    05/16/2023
# Responsible(s): Alessandro Marins and Thiago Pena
##############################

import subprocess
import time
import os
from configparser import ConfigParser
import random
import numpy as np

def main(user_SD_name=None, file_pathdir=None,
	first_part_name     = None, file_extension="sh", horn_pos_first=None, horn_pos_last=None, 
	job_first_part_name = None, steptime=300):
    
	Pos_index       = np.arange(horn_pos_first,horn_pos_last+1,1)
	vfilenames      = ["{}_{}.{}".format(first_part_name, i, file_extension) for i in Pos_index] #Vector of sbatch filenames. For instance: [sbatch_HIDE_position_1.sh, sbatch_HIDE_position_2.sh,..., sbatch_HIDE_position_140.sh]
	for i,ifilename in enumerate(vfilenames, start=horn_pos_first): #acredito que essa parte precise ser modificada para funcionar corretamente com horn_pos_first!=1
		isbatch = "/".join(("sbatch {}","{}")).format(file_pathdir, ifilename) #string with the command to execute i-th sbatch file. For instance: sbatch /scratch/bingo/thiago.pena/sbatch_files/sbatch_HIDE_position_6.sh
		ijob    = "{}_{}".format(job_first_part_name,i)                      #job name to be identified on SD list
		squeue  = "squeue -u {} | grep {}".format(user_SD_name, ijob)          #command to identify i-th job by a USER on SD system
		print(i)
		print(horn_pos_first)
		if i==horn_pos_first: 
			start = str(time.ctime())                                    #just to save the first time running
		print('Running       : {}'.format(ifilename))
		print('Sbatch script : {}'.format(isbatch))
		print('Run {} job at: {}\n'.format( ijob, time.ctime()))
		#os.system(isbatch) #submit a string with a command to submit a sbatch in a row
		time.sleep(2)
		keep = True        #condition to keep waiting the job to finish
		
		if i> horn_pos_first:
			print(i)
			print(horn_pos_first)
			prev_job_name = "{}_{}".format(job_first_part_name, i-1)
			prev_complete = False
			while not prev_complete:
		#while keep:        #just will left when there is no job in a SD row.
				stream = os.popen("squeue -u {} | grep {}".format(user_SD_name, prev_job_name))
				if stream.read() == "": #when you look for a job in SD list and there is nothing you get 'stream.read()' as "". Then, you'll need rerun
		#keep=False
					prev_complete = True
				else:                   #if there is still a job running, just getting wait
					print('Still running...')
					time.sleep(timestep)  #step time to (re)check
	
		os.system(isbatch)
		



	

	print('Its over.')
	print('Started at ' + start)
	print('Finished at ' + str(time.ctime()))
	print('A total of {} submissions'.format(i+1))
	return 0

if __name__ == "__main__":
	############
	#Parameters
	############
	user_SD_name        = "thiago.pena"                              #USER name on SDumont Cluster to be identified by 'squeue' command
	file_pathdir        = "/scratch/bingo/thiago.pena/sbatch_files"  #Path directing where there are ".sh" files to be submitted
	first_part_name     = "sbatch_HIDE_position"                     #First name part of ".sh" names. For instance, a filename as "sbatch_HIDE_position_1.sh", the first name part is "sbatch_HIDE_position", "1" is the position index, and "sh" is the extension of the file
	file_extension      = "sh"                                       #Extension of the file
	horn_pos_first      = 1                                          #First position index to be submitted. Ps.: each position index doenst correspond to a horn but a specific pointing on the sky. This is because one horn could pinpoint more than one position on the sky
	horn_pos_last       = 140                                        #Last position index to be submitted.
	job_first_part_name = "ID"                                  #Identification part of the jobs to be runned. For instance, if you will submit jobs as "#SBATCH --job-name=" "HIDE_ID_1", "HIDE_ID_2", ..., "HIDE_ID_140", 'job_first_part_name' will be 'HIDE_ID'. Numbers will correspond to position index.
	timestep            = 3                                        #step time btw two checking [in sec]. Comment: Once this code aims to run to submit ".sh" in sequence and SDumont cluster only accept one job per time, we need have another code such can check if there are any code running and if not will submit a new one. Each checking will be after steptime seconds
	
	############
	# Running
	############
	main(user_SD_name, file_pathdir,
	first_part_name, file_extension, horn_pos_first, horn_pos_last, 
	job_first_part_name, timestep)
