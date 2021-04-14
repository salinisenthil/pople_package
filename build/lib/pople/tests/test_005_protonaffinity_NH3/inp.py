import os
from pople import calculator as calc
from pople import au2ev, au2kcm, au2kjm

#=== Remove old files
os.system("rm -f ORCA.com ORCA.out Thermochemistry.out input.*")

exe='/home/Lib/ORCA_420/orca'

#=== NH3
geom = '''
           0           1
 N        0.00000000     0.00000000     0.00000000
 H       -0.39793000     0.00000000    -0.93852000
 H       -0.39793000    -0.81279000     0.46926000
 H       -0.39793000     0.81279000     0.46926000
'''

out = calc(code='orca', code_exe=exe, method='g4mp2-xp', xyz=geom)
H_mol=out[2]

#=== NH4+
geom = '''
           1           1
 N        0.00000000     0.00000000     0.00000000
 H        0.59435200     0.59435200     0.59435200
 H       -0.59435200    -0.59435200     0.59435200
 H        0.59435200    -0.59435200    -0.59435200
 H       -0.59435200     0.59435200    -0.59435200
'''

out = calc(code='orca', code_exe=exe, method='g4mp2-xp', xyz=geom)
H_cation_protonated=out[2]

PA = H_mol - H_cation_protonated

print(' Standard proton affinity (298.15 K) of NH3 is ', (PA), ' hartree')
print(' Standard proton affinity (298.15 K) of NH3 is ', (PA)*au2ev, ' eV')
print(' Standard proton affinity (298.15 K) of NH3 is ', (PA)*au2kcm, ' kcal/mol')
print(' Standard proton affinity (298.15 K) of NH3 is ', (PA)*au2kjm, ' kJ/mol')

#Reference output
#Standard proton affinity (298.15 K) of NH3 is  0.322931995145062  hartree
#Standard proton affinity (298.15 K) of NH3 is  8.787427911725857  eV
#Standard proton affinity (298.15 K) of NH3 is  202.64291773307903  kcal/mol
#Standard proton affinity (298.15 K) of NH3 is  847.8579677952026  kJ/mol
#Exp. value is 202.5 kcal/mol

