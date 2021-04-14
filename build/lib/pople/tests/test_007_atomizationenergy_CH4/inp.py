import os
from pople import calculator as calc
from pople import au2kcm
from pople import getatoms, uniqatoms, getmultip

#=== Remove old files
os.system("rm -f ORCA.com ORCA.out Thermochemistry.out input.*")

prog='orca'
exe='/home/Lib/ORCA_420/orca'
meth='g4mp2'

#=== CH4 
geom = '''
 0  1
 C        0.00000000     0.00000000     0.00000000
 H        1.09336000     0.00000000     0.00000000
 H       -0.36445000     0.00000000    -1.03083000
 H       -0.36445000    -0.97188000     0.34361000
 H       -0.36445000     0.97188000     0.34361000
'''

out = calc(code=prog, code_exe=exe, method=meth, xyz=geom)
U0_mol = out[0]         # Internal energy (at 0 K) of molecule




#=== Atoms
atoms = getatoms(geom)  # List of atoms

uniq = uniqatoms(atoms) # Unique atom data

N_ua = uniq["N_ua"]       # No. of unique atoms
ua   = uniq["uniq_sym"]   # unique atoms
ua_N = uniq["uan"]        # unique atom frequencies

multip=getmultip(ua)    # Multiplicities of unique atoms

U0_atoms_total=0        # Sum of internal energy (at 0 K) of atoms
for i_ua in range(N_ua):

    print(' Atom ', ua[i_ua], ' with multiplicity ', multip[i_ua])

    #=== calculate for each unique atom
    geom = (f'''0  {multip[i_ua]} 
{ua[i_ua]}  0.0 0.0 0.0 ''')

    out = calc(code=prog, code_exe=exe, method=meth, xyz=geom)
    U0=out[0]

    U0_atoms_total = U0_atoms_total + U0 * ua_N[i_ua] 




#=== Atomization energy
print( ' Atomization energy of CH4 is ', (U0_atoms_total-U0_mol) * au2kcm, ' kcal/mol')

#Reference output
#Atomization energy of CH4 is  392.2148180445219  kcal/mol
#Exp. value is 397.94 kcal.mol

