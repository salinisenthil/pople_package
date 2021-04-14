kB     = 1.38064852000000*(10**-23) 
h      = 6.626070040*(10**-34)    
Ry     = 10973731.568508  
c      = 299792458    
T      = 298.15  
e      = 1.60217648700000*(10**-19)   
N_avo  = 6.02214179000000*(10**+23)  

au2j   = 2 * Ry * h * c
j2au   = 1 / au2j

au2kcm = 2 * Ry * h * c * N_avo / 4184.0
kcm2au = 1 / au2kcm

au2cmi = 2 * Ry / 100.0
cmi2au = 1 / au2cmi

au2ev  = 2 * Ry * h * c / e
au2kjm = 2 * Ry * h * c * N_avo / 1000.0

from  .utilities          import  start_time, final_time, print_header
from  .params_atomic      import  atno, sym2mass, nanb, NFC
from  .params_thermal     import  dH_atoms, HOF_atoms
from  .params_spinorbit   import  At_SO, Mol_SO
from  .energy_corrections import  thermal, HLC_g4mp2
from  .uniqatoms          import  uniqatoms, getatoms, getmultip
from  .principal_coord    import  principal_coord
from  .orca_g4mp2_ctrl    import  orca_g4mp2_ctrl
from  .orca_printbas      import  orca_printbas
from  .orca_run           import  orca_run 
from  .orca_g4mp2         import  orca_g4mp2
from  .calculator         import  calculator

#print("LOADED REQUIREMENTS FOR POPLE2")
