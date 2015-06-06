import os
import sys

CorrNames=['noCorrTotalRegion','CorrTotalRegion']
#CorrNames=['CorrTotalRegion']
sampleName = ['Muon']
dirName = ['Zpt_ZmassPlotsPtBins_Gaus']
for Corr in CorrNames:
  for Dir in dirName:
    for Sample in sampleName:
      cmd_string = "root -l -b -q Zpt_ZmassCompPtBins_Gaus.C+\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\)" %(Sample,Corr,Dir)
      os.system(cmd_string)
      
      cmd_string = "rm -f *.d *.so"
      os.system(cmd_string)
