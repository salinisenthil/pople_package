import os
from pople import calculator as calc
from pople import au2ev

#=== Remove old files
os.system("rm -f ORCA.com ORCA.out Thermochemistry.out input.*")

exe='/home/Lib/ORCA_420/orca'

#=== neutral
geom_1 = '''
0  2
Cl  0.0  0.0  0.0
'''
#=== anion
geom_2 = '''
-1  1
Cl  0.0  0.0  0.0
'''

out1 = calc(code='orca', code_exe=exe, method='g4mp2', xyz=geom_1)
out2 = calc(code='orca', code_exe=exe, method='g4mp2', xyz=geom_2)

U0_neutral=out1[0]
U0_anion=out2[0]

EA = (U0_neutral-U0_anion) * au2ev

print(' Electron affinity of Cl atom is is ', EA,' eV')

#Reference output
#Electron affinity of Cl atom is is  3.6421649462081915  eV
#Exp. value is 3.6166 eV

