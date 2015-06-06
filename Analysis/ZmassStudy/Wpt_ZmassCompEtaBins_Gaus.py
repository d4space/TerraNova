import os
import sys

CorrNames=['noCorrTotalRegion','CorrTotalRegion']
sampleName = ['Muon']
#dirName = ['Wpt_ZmassPlotsEtaBins_Gaus','Wpt_ZmassPlotsEtaBins_noOverLap_Gaus','Wpt_ZmassPlotsEtaBins_LeadingLept_noOverLap_Gaus','Wpt_ZmassPlotsEtaBins_LeadingLept_Gaus','Wpt_ZmassPlotsEtaBins_TrailingLept_noOverLap_Gaus','Wpt_ZmassPlotsEtaBins_TrailingLept_Gaus']
dirName = ['Wpt_ZmassPlotsEtaBins_Gaus']
for Corr in CorrNames:
  for Dir in dirName:
    for Sample in sampleName:
      cmd_string = "root -l -b -q Wpt_ZmassCompEtaBins_Gaus.C+\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\)" %(Sample,Corr,Dir)
      os.system(cmd_string)
      
      cmd_string = "rm -f *.d *.so"
      os.system(cmd_string)
