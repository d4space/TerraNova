import os
#import sys
BaseName="WpToMuNu"
cmd_string = "root -l -q resbosComp.C+\(\\\"%s\\\"\)" %BaseName
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)

BaseName="WmToMuNu"
cmd_string = "root -l -q resbosComp.C+\(\\\"%s\\\"\)" %BaseName
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)

BaseName="WInclToMuNu"
cmd_string = "root -l -q resbosComp.C+\(\\\"%s\\\"\)" %BaseName
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)

BaseName="WpToEleNu"
cmd_string = "root -l -q resbosComp.C+\(\\\"%s\\\"\)" %BaseName
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)

BaseName="WmToEleNu"
cmd_string = "root -l -q resbosComp.C+\(\\\"%s\\\"\)" %BaseName
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)

BaseName="WInclToEleNu"
cmd_string = "root -l -q resbosComp.C+\(\\\"%s\\\"\)" %BaseName
os.system(cmd_string)

cmd_string = "rm -f *.d *.so"
os.system(cmd_string)
