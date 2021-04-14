import os
from pople import calculator as calc
from pople import au2ev

#=== Remove old files
os.system("rm -f ORCA.com ORCA.out Thermochemistry.out input.*")

geom = '''
0  1
H  0.0  0.0  0.0
H  0.0  0.0  0.7
'''

exe='/home/Lib/ORCA_420/orca'

out = calc(code='orca', code_exe=exe, method='g4mp2', xyz=geom)

H=out[2]

print(' Standard enthalpy (at 298.15 K) is ', H,' hartree')

#Reference output
#Standard enthalpy (at 298.15 K) is  -1.1671036392034166  hartree


