import os
import sys

CorrNames=['noCorrTotalRegion','CorrTotalRegion']
#CorrNames=['noCorrTotalRegion']
sampleName = ['Muon']
dirName = ['Zpt_ZmassPlotsEtaBins_Gaus']
for Corr in CorrNames:
  for Dir in dirName:
    for Sample in sampleName:
      cmd_string = "root -l -b -q Zpt_ZmassCompEtaBins_Gaus.C+\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\)" %(Sample,Corr,Dir)
      os.system(cmd_string)
      
      cmd_string = "rm -f *.d *.so"
      os.system(cmd_string)
