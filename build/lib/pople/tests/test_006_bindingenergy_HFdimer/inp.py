import os
from pople import calculator as calc
from pople import au2kcm

#=== Remove old files
os.system("rm -f ORCA.com ORCA.out Thermochemistry.out input.*")

exe='/home/Lib/ORCA_420/orca'

#=== HF-dimer
dimer = '''
0  1
 F        0.00000000     0.00000000     0.00000000
 H        0.92614000     0.00000000     0.00000000
 F        1.72307000     0.00000000    -1.82971000
 H        0.79706000    -0.00010000    -1.83479000
'''

out = calc(code='orca', code_exe=exe, method='g4mp2-xp', xyz=dimer)
U0_A2=out[0]

#=== cation H2O
geom = '''
0  1
 F        0.00000000     0.00000000     0.00000000
 H        0.93389000     0.00000000     0.00000000
'''

out = calc(code='orca', code_exe=exe, method='g4mp2-xp', xyz=geom)
U0_A=out[0]

print(' Binding energy of HF-dimer is ', (U0_A2 - 2*U0_A)*au2kcm, ' kcal/mol')

#Reference output
#Binding energy of HF-dimer is  0.3457599563589871  kcal/mol
#Exp. value is -2.97 kcal/mol
#G4(MP2) value is 0.25 kcal/mol
