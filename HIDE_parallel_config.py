import numpy as np
import os, sys
import configparser
import tempfile
import uuid
import shutil
import subprocess
import multiprocessing as mp
from pathlib import Path

'''
Created on January 7, 2026

authors: Nicolli Soares

'''
#=========================================================================================================
# =================================== HIDE Parallel Configure ============================================
#=========================================================================================================

def init_worker():
	"""
	Initializer for multiprocessing workers.  
	Configures HOPE with a unique prefix per process to avoid conflicts.
	"""
	import hope

	hope_dir = os.path.abspath(hope.config.prefix)
	if os.path.exists(hope_dir):
		shutil.rmtree(hope_dir, ignore_errors=True)
		
	hope.config.prefix = f".hope.{sys.version_info[0]}.{sys.version_info[1]}.{mp.current_process().pid}"

	hope.config.keeptemp = False

def run_horn(i, working_path, destination_path, dfile_short):
	"""
	Runs a horn with an isolated environment.
	
	Args:
	i: Horn number
	working_path: Working directory with input files
	destination_path: Destination directory for output
	dfile_short: Base name of the configuration file
	"""
	
	# Create temporary directories
	tag = f"horn_{i}"
	home_dir = Path(tempfile.gettempdir()) / f"hide_home_{tag}_{os.getpid()}_{uuid.uuid4().hex[:6]}"
	cache_dir = home_dir / ".cache"
	tmp_dir = Path(tempfile.gettempdir()) / f"hide_tmp_{tag}_{os.getpid()}_{uuid.uuid4().hex[:6]}"
	
	home_dir.mkdir(parents=True, exist_ok=True)
	cache_dir.mkdir(exist_ok=True)
	tmp_dir.mkdir(exist_ok=True)
	
	# Configure a modified environment
	env = os.environ.copy()
	env['HOME'] = str(home_dir)
	env['XDG_CACHE_HOME'] = str(cache_dir)
	env['TMPDIR'] = str(tmp_dir)
	
	# print(f"Running horn {i} in temporary cache directory {cache_dir}")

	subprocess.run(['cp', f'{working_path}bingo_horn_{i}.py', destination_path], check=True)
	subprocess.run(['hide', f'hide.config.{dfile_short}_{i}'], env=env, check=True)
	
	shutil.rmtree(home_dir, ignore_errors=True)
	shutil.rmtree(tmp_dir, ignore_errors=True)

#=========================================================================================================