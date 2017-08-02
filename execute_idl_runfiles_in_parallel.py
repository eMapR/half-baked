# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 10:07:19 2017

@author: braatenj
"""

import subprocess
import multiprocessing


def run_cmd(cmd):
  print(cmd)
  subprocess.call(cmd, shell=True)

#################################################################################################### 
####################################################################################################
####################################################################################################

runFiles = [
  '/vol/v1/proj/nccn/2017/scripts/02_v1_lewi_labfilt_run.pro',
  '/vol/v1/proj/nccn/2017/scripts/02_v1_olym_labfilt_run.pro',
  '/vol/v1/proj/nccn/2017/scripts/02_v1_noca_labfilt_run.pro',
  '/vol/v1/proj/nccn/2017/scripts/02_v1_mora_labfilt_run.pro'
]

#################################################################################################### 
####################################################################################################
####################################################################################################

# make a list of commands
cmds = ['idl -e @'+runFile for runFile in runFiles]
  
# run the commands in parallel
processes = len(cmds) if len(cmds)<=10 else 10
pool = multiprocessing.Pool(processes=processes)
pool.map(run_cmd, cmds)  
pool.close()