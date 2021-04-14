import os
from pople import au2kcm
from pople import getatoms, uniqatoms, getmultip
from pople import HOF_atoms, dH_atoms

#=== CH4 
geom = '''
 0  1
 C        0.00000000     0.00000000     0.00000000
 H        1.09336000     0.00000000     0.00000000
 H       -0.36445000     0.00000000    -1.03083000
 H       -0.36445000    -0.97188000     0.34361000
 H       -0.36445000     0.97188000     0.34361000
'''

#=== from test_007_atomizationenergy_CH4/Thermochemistry.out 
HTmol = -40.42387018  
U0mol = -40.42768573

#=== from test_007_atomizationenergy_CH4/out
U0molAE = 392.2148180200836 / au2kcm


#=== Atoms
atoms = getatoms(geom)  # List of atoms

uniq = uniqatoms(atoms) # Unique atom data

N_ua = uniq["N_ua"]       # No. of unique atoms
ua   = uniq["uniq_sym"]   # unique atoms
ua_N = uniq["uan"]        # unique atom frequencies

multip=getmultip(ua)    # Multiplicities of unique atoms


#=== Enthalpy of formation
HT_form = HTmol - U0mol - U0molAE

for i_ua in range(N_ua):

   HOF_0K     = HOF_atoms(ua[i_ua])
   dH_298K_0K = dH_atoms(ua[i_ua])

   HT_form = HT_form + ua_N[i_ua] * ( HOF_0K - dH_298K_0K )

print( ' Standard formation enthalpy (at 298.15 K) of CH4 is ', HT_form * au2kcm, ' kcal/mol')

#Reference output
#Standard formation enthalpy (at 298.15 K) of CH4 is  -17.61052387648575  kcal/mol
#Exp. value is -17.90

