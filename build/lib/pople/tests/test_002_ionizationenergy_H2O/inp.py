import os
from pople import calculator as calc
from pople import au2ev

#=== Remove old files
os.system("rm -f ORCA.com ORCA.out Thermochemistry.out input.*")

#=== Common inputs
exe='/home/Lib/ORCA_420/orca'

#=== neutral H2O
geom = '''
0  1
O      0.0         0.0      0.0
H      0.96854     0.0      0.0
H     -0.22812     0.0     -0.94129
'''

out = calc(code='orca', code_exe=exe, method='g4mp2', xyz=geom)
U0_neutral=out[0]

#=== cation H2O
geom = '''
+1  2
O      0.0         0.0      0.0
H      1.01249     0.0      0.0
H     -0.34159     0.0     -0.95313
'''

out = calc(code='orca', code_exe=exe, method='g4mp2', xyz=geom)
U0_cation=out[0]

print(' Ionization potential of H2O is ', (U0_cation-U0_neutral)*au2ev, ' eV')

#Reference output
#Ionization potential of H2O is  12.590811264819967  eV
#Exp. value is 12.619 eV

