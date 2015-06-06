import os
import sys

dirName = 'Wpt_muPtPlots'
cmd_string = "root -l -b -q Wpt_muPtComp.C+\(\\\"%s\\\"\)" %(dirName)
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)
