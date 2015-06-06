import os
import sys

CorrNames=['noCorrTotalRegion','noCorrBB','noCorrBE','noCorrEE','CorrTotalRegion','CorrBB','CorrBE','CorrEE']
#CorrNames=['noCorrBB']
#CorrNames=['CorrTotalRegion']
sampleName = ['Muon']
dirName = ['Wpt_ZmassPlots_Gaus']
for Corr in CorrNames:
  for Dir in dirName:
    for Sample in sampleName:
      cmd_string = "root -l -b -q Wpt_ZmassComp_Gaus.C+\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\)" %(Sample,Corr,Dir)
      os.system(cmd_string)
      
      cmd_string = "rm -f *.d *.so"
      os.system(cmd_string)
