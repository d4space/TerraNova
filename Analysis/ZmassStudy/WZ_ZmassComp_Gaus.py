import os
import sys

CorrNames=['noCorrTotalRegion','CorrTotalRegion']
#CorrNames=['noCorrTotalRegion']
sampleName = ['Muon']
dirName = ['WZ_ZmassPlots_Gaus']
for Corr in CorrNames:
  for Dir in dirName:
    for Sample in sampleName:
      cmd_string = "root -l -b -q WZ_ZmassComp_Gaus.C+\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\)" %(Sample,Corr,Dir)
      os.system(cmd_string)
      
      cmd_string = "rm -f *.d *.so"
      os.system(cmd_string)
