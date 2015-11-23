# Edits the FORTRAN source code to incorporate
# the custom model of the user

import string
import sys
import re

with open("../core/lin_mod.f90") as f:
    for num,line in enumerate(f,start=0):
        if '@configOnB' in line:
            bStart=num
        if '@configOffB' in line:
            bEnd=num
        if '@configOnF' in line:
            FStart=num
        if '@configOffF' in line:
            FEnd=num
        if '@configOnH' in line:
            HStart=num
        if '@configOffH' in line:
            HEnd=num
with open("../core/lin_mod.f90") as f:
    sourceFile=f.readlines()

# parse the config file as set by the user
with open(sys.argv[1]) as f:
    configFile=f.readlines()



nF=0
nb=0
nH=0
k=0
while k<len(configFile):
    if configFile[k][0] == 'F':
        nF+=1
        sourceFile.insert(FStart+nF,configFile[k])
    if configFile[k][0] == 'b':
        nb+=1
        sourceFile.insert(bStart+nb+nF,configFile[k])
    if configFile[k][0] == 'H':
        nH+=1
        sourceFile.insert(HStart+nb+nF+nH,configFile[k])
    k+=1
if nF != nb**2:
    print('F and b are not compatibile')
    sys.exit(0)

with open("../core/lin_mod_changed.f90",'w') as f:
    f.writelines(sourceFile)
