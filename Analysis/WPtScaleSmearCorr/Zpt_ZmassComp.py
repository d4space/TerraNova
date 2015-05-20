import os
import sys

#CorrNames=['noCorrTotalRegion','noCorrBB','noCorrBE','noCorrEE','CorrTotalRegion','CorrBB','CorrBE','CorrEE']
CorrNames=['noCorrTotalRegion','CorrTotalRegion']
sampleName = ['Muon']
dirName = ['Zpt_ZmassPlots_Gaus_MoreBins']
for Corr in CorrNames:
  for Dir in dirName:
    for Sample in sampleName:
      cmd_string = "root -l -b -q Zpt_ZmassComp_Gaus.C+\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\)" %(Sample,Corr,Dir)
      os.system(cmd_string)
      
      cmd_string = "rm -f *.d *.so"
      os.system(cmd_string)
