#!/usr/bin/python

import os, math, sys
from time import strftime

import argparse
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--run", dest='run', default="1", help="Run number to submit")
args = parser.parse_args()

test = False # does not run the jobs
if not test:
  if not os.path.exists("gridOutput"):
    os.system("mkdir gridOutput")
  if not os.path.exists("gridOutput/gridJobs"):
    os.system("mkdir gridOutput/gridJobs")

config_name = "$ROOTCOREBIN/data/DijetResonanceAlgo/DijetResonanceAlgo.config"

production_name = "phys-exotics"

#output_tag = "Group"+args.run
#file_in = "DijetResonanceAlgo/scripts/sampleLists/GroupRuns/GroupRun_"+run

outputTags = []
filesIn = []
sysList = []

outputTags.append("MainE2_C_v2")  # CHANGE BEFORE SUBMITTING - no periods!
filesIn.append("DijetResonanceAlgo/scripts/sampleLists/Data13TeV_Main_EXOT2_gridSamples.txt")
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/Data13TeV_Main_gridSamples.txt")
sysList.append( { 'Nominal' : "0" } )

#outputTags.append("Main_C2_fGRL_Calib2")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/Data13TeV_Main_EXOT2_C2_gridSamples.txt")
#outputTags.append("Main_C3_fGRL_Calib2")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/Data13TeV_Main_EXOT2_C3_gridSamples.txt")
#outputTags.append("Main_C4_fGRL_Calib2")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/Data13TeV_Main_EXOT2_C4_gridSamples.txt")
#sysList.append( { 'Nominal' : "0"  } )

#outputTags.append("QStar")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append( "DijetResonanceAlgo/scripts/sampleLists/MC15a_EXOT2_QStar_gridSamples.txt" )
#outputTags.append("QStar_FS")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append( "DijetResonanceAlgo/scripts/sampleLists/MC15a_QStar_FullSim_gridSamples.txt")
#sysList.append({ 'All' : "0.5,1,1.5,2,2.5,3.0" } )

#outputTags.append("QBH_v3")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append( "DijetResonanceAlgo/scripts/sampleLists/QBH_gridSamples.txt" )
#sysList.append( { 'All' : "0.5,1,1.5,2,2.5,3.0" } )

#outputTags.append("BlackMax")
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/BlackMax_gridSamples.txt")
#sysList.append( { 'All' : "0.5,1,1.5,2,2.5,3.0" } )

#outputTags = "Week1"  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn = "DijetResonanceAlgo/scripts/sampleLists/Week1_EXOT2_gridSamples.txt"
#sysList.append( { 'All' : "1"  } )

outputTags.append("QCD_v2")  # CHANGE BEFORE SUBMITTING - no periods!
filesIn.append("DijetResonanceAlgo/scripts/sampleLists/MC15a_EXOT2_gridSamples.txt")
#sysList.append( { 'All' : "1"  } )
sysList.append( { 'Nominal' : "0" } )
#
##outputTags.append("CI")  # CHANGE BEFORE SUBMITTING - no periods!
##filesIn.append("DijetResonanceAlgo/scripts/sampleLists/MC15_CI_gridSamples.txt")
#sysList.append( { 'All' : "1"  } )
#
#outputTags.append("Powheg")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/MC15a_EXOT2_Powheg_gridSamples.txt")
#sysList.append( { 'All' : "1"  } )
#
#outputTags.append("PowhegHerwig")  # CHANGE BEFORE SUBMITTING - no periods!
#filesIn.append("DijetResonanceAlgo/scripts/sampleLists/MC15a_PowhegHerwig_grid.txt")
#sysList.append( { 'All' : "1"  } )



timestamp = strftime("_%Y%m%d")

for iFile, file_in in enumerate(filesIn):
  systDict = sysList[iFile]
  output_tag = outputTags[iFile]

  output_tag += timestamp
  ##### list of systematics
  #if "QStar" in file_in or "QBH" in file_in or "BlackMax" in file_in:
  #  systDict = { 'All' : "0.5,1,1.5,2,2.5,3.0" }  # run all systematics up/down 3 sigma and nominal
  #elif "Data13TeV" in file_in:
  #  systDict = { 'Nominal' : "0" }  # run nominal only
  #else:
  #  systDict = { 'All' : "1"  }  # run all systematics up/down 3 sigma and nominal



  #############
  # run over systematic variations
  #############
  for systName, value in systDict.iteritems():
    syst_tag = systName
    if not systName == "Nominal" and not "," in value:
      syst_tag += '_'+value
      syst_tag = syst_tag.replace(".","-") # do not put more periods in name
    send_output_tag = syst_tag+"."+output_tag
    print syst_tag
    print output_tag
    print send_output_tag

    submit_dir = "gridOutput/gridJobs/submitDir_"+os.path.basename(file_in).rsplit('.',1)[0]+"."+send_output_tag # define submit dir name

    command = "runDijetResonance -inFile "+file_in+" -outputTag "+send_output_tag+" -submitDir "+submit_dir+" -configName "+config_name+" -syst "+systName+" "+str(value)
    if len(production_name) > 0:
      command += ' -production '+production_name
    print command
    if not test: os.system(command)


