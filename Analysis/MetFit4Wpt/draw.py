import os
import glob
import subprocess
import sys
HFLUMI =0.018729
PIXLUMI=0.018479
#WPTLUMI=0.018977
WPTLUMI=0.018429
LUMI=WPTLUMI

#cmd_string = "root -l -q Wpt_PASformat.C+\(\\\"Wpt_plots\\\",\\\"Muon\\\",%f,0\)" %LUMI
#os.system(cmd_string)
#cmd_string = "rm *.so *.d"
#print "command %s is running" %cmd_string
#os.system(cmd_string)

#cmd_string = "root -l -q Wpt_PASformat.C+\(\\\"Wpt_plots\\\",\\\"Electron\\\",%f,0\)" %LUMI
#os.system(cmd_string)
#cmd_string = "rm *.so *.d"
#print "command %s is running" %cmd_string
#os.system(cmd_string)


cmd_string = "root -l -q Wpt_PASformat_withRatio.C+\(\\\"Wpt_plots\\\",\\\"Muon\\\",%f,0\)" %LUMI
os.system(cmd_string)
cmd_string = "rm *.so *.d"
print "command %s is running" %cmd_string
os.system(cmd_string)

cmd_string = "root -l -q Wpt_PASformat_withRatio.C+\(\\\"Wpt_plots\\\",\\\"Electron\\\",%f,0\)" %LUMI
os.system(cmd_string)
cmd_string = "rm *.so *.d"
print "command %s is running" %cmd_string
os.system(cmd_string)
