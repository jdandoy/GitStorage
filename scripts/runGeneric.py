#!/usr/bin/env python

#######################################################################################
# runSingleFit.py
#
# This is a script for running singleFit.py on all the expected mjj histograms.
# Lists to loop over include:
#   files: a set of different input files
#   fitsToUse: The fits to use
#   systematics: Systematic directories in the files
#   histNameVariations:  Potential variations to the input histograms (i.e. luminosity)
#   fitVariations:  Potential variations to fits (i.e. tolerance)
#       This variation should both change the fitting and include an addition
#       to the output name.
#
#######################################################################################


import subprocess, os, time, sys, glob
import datetime

def main():

  pids = []
  logFiles = []
  NCORES = 5
  if not os.path.exists("parallelLogs/"):
    os.makedirs("parallelLogs/")

  #######   Setup the list of commands  ######
  runCommand = sys.argv[1]
  sys.argv = sys.argv[2:]

  iArg = 0
  boolOptions, valOptions, values = [],[],[]
  while iArg < len(sys.argv):
    if not sys.argv[iArg].startswith('-'):
      print "Error, arguement ", sys.argv[iArg], " doesn't start with '-'"
      exit(1)
    if iArg == len(sys.argv)-1:
      boolOptions.append( sys.argv[iArg] )
      iArg+=1
    elif sys.argv[iArg+1].startswith('-'):
      boolOptions.append( sys.argv[iArg] )
      iArg+=1
    else:
      valOptions.append( sys.argv[iArg] )
      values.append( sys.argv[iArg+1].split(',') )
      iArg+=2

  print valOptions
  print values
  ## Add all option strings to their option values ##
  commandsList = []
  for iOption, option in enumerate(valOptions):
    commandsList.append( [] )
    for iValue, value in enumerate(values[iOption]):
      commandsList[iOption].append( option+' '+str(value) )


  import itertools
  while len(commandsList) > 1:
    list1 = commandsList.pop(0)
    list2 = commandsList.pop(0)
    commands =  map(' '.join, itertools.chain(itertools.product(list1, list2), itertools.product(list2, list1)))
    commands = [x for x in commands if not any( x.startswith(y) for y in list2 )]
    commandsList.append( commands )


  commandsList = commandsList[0]
  for iC, command in enumerate(commandsList):
    commandsList[iC] = runCommand+' '+command
    for boolOpt in boolOptions:
      commandsList[iC] += ' '+boolOpt


  print commandsList
  exit(1)

  for iC, command in enumerate(commandsList):
    while len(pids) >= NCORES:
      wait_completion(pids, logFiles)

    logFile = "parallelLogs/runLog_"+str(iC)+'.txt'
    res = submit_local_job( command, logFile )
    pids.append(res[0])
    logFiles.append(res[1])

  wait_all(pids, logFiles)


def submit_local_job(exec_sequence, logfilename):
  output_f=open(logfilename, 'w')
  pid = subprocess.Popen(exec_sequence, shell=True, stderr=output_f, stdout=output_f)
  time.sleep(0.1)  #Wait to prevent opening / closing of several files

  return pid, output_f

def wait_completion(pids, logFiles):
  print """Wait until the completion of one of the launched jobs"""
  while True:
    for pid in pids:
      if pid.poll() is not None:
        print "\nProcess", pid.pid, "has completed"
        logFiles.pop(pids.index(pid)).close()  #remove logfile from list and close it
        pids.remove(pid)

        return
    print ".",
    sys.stdout.flush()
    time.sleep(3) # wait before retrying

def wait_all(pids, logFiles):
  print """Wait until the completion of all remaining launched jobs"""
  while len(pids)>0:
    wait_completion(pids, logFiles)
  print "All jobs finished!"

if __name__ == "__main__":
  beginTime = datetime.datetime.now().time()
  main()
  print "Starting fitting at", beginTime
  print "Finished fitting at", datetime.datetime.now().time()



