import os

#Get lastest GRL run
#wget http://atlasdqm.web.cern.ch/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v63-pro18-01_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml
os.system('wget http://atlasdqm.web.cern.ch/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v64-pro19_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml -O data15_13TeV.periodAllYear_DetStatus-v64-pro19_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml')
#os.system('wget http://atlasdqm.web.cern.ch/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v63-pro18-01_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml -O data15_13TeV.periodAllYear_DetStatus-v63-pro18-01_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml')
GRLFile = 'data15_13TeV.periodAllYear_HEAD_DQDefects-00-01-02_PHYS_StandardGRL_Atlas_Ready.xml'
os.system('wget http://atlasdqm.web.cern.ch/atlasdqm/grlgen/Atlas_Ready/'+GRLFile+' -O '+GRLFile)

Files = ["Data13TeV_Main_gridSamples.txt",
    "Data13TeV_Main_EXOT2_D1_gridSamples.txt",
    "Data13TeV_Main_EXOT2_D2_gridSamples.txt",
    "Data13TeV_Main_EXOT2_D3_gridSamples.txt",
    "Data13TeV_Main_EXOT2_D4_gridSamples.txt",
    "Data13TeV_Debug_gridSamples.txt"
    ]

GRLNumbers = [line for line in open(GRLFile, "r") if 'RunList' in line]
if not len(GRLNumbers) == 1:
  print "Error, more than one RunList in the GRL file??"
  exit(1)

GRLNumbers = GRLNumbers[0].replace('<Metadata Name="RunList">','').replace('</Metadata>','').rstrip().split(',')

# Check for new runs in the GRL
missingRun = False
for GRLNumber in GRLNumbers:
  if not GRLNumber in open(Files[0]).read():
    missingRun = True
    print "------Missing new file ", GRLNumber
if missingRun == False:
  print "------No new runs exist in GRL\n"

# Check if any missing Derivation runs have been finished
for file in Files:
  print "checking ", file
  with open( file , "r") as inFile:
    for line in inFile:
      if line.startswith('#?'):
        cont = line[2:]
        print "-----Checking ", cont
        os.system('rucio list-dids '+cont)
