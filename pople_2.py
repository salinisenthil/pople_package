import os, sys, string
import linecache, math
import numpy as np
import datetime , time


current_dir_path = os.getcwd()

input_f = current_dir_path + "/initial.inp"
job_f = current_dir_path + "/job.inp"
xyz_path = current_dir_path + "/geom.xyz"

### declaring units - S
kB = 1.38064852000000*(10**-23) #kB = 1.38064852000000d-23
h = 6.626070040*(10**-34)    #6.626070040d-34
Ry = 10973731.568508  #10973731.568508d0
c = 299792458    #299792458d0
T = 298.15  #298.15d0 
e = 1.60217648700000*(10**-19)   #1.60217648700000d-19
N_avo = 6.02214179000000*(10**+23)  #6.02214179000000d+23
au2j = 2 * Ry * h * c

au2ev  = 2 * Ry * h * c          / e
au2kjm = 2 * Ry * h * c * N_avo / 1000
au2kcm = 2 * Ry * h * c * N_avo / 4184
kcm2au = 1 / au2kcm

j2au   = 1 / au2j
au2cmi = 2 * Ry / 100
cmi2au = 1 / au2cmi
### declaring units - E

###### Mol_SO - S
def Mol_SO(Nat, multip, charge, sym, SO_3rdrow_mols_val):   # number of atoms, multiplicity, charge, array of atoms in molecule, value of SO_3rdrow_mols (from orca.inp file)
    Mol_SO = 0

    # Special Case - Acetleyne - S
    if Nat == 4 and multip == 2 and charge == 1:
        countH_temp =0
        countC_temp =0
        for tmp in range(len(sym)):
            if sym[tmp] == "H":
                countH_temp= countH_temp +1
            if sym[tmp] == "C":
                countC_temp = countC_temp +1
        if countH_temp == 2 and countC_temp == 2:
            Mol_SO = -0.07  #-0.07d0

        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("DETECTED A C2H2+ SYSTEM: Using SO parameters for acetylene cation\n")
            ther_chem.write("Ref: JCP 114, 9287, 2001\n\n")
    # Special Case - Acetleyne - E

    # For diatomics with multip = 2
    if Nat == 2 and multip == 2 :
        sort_sym = sorted(sym, reverse=True)
        if SO_3rdrow_mols_val == "true":  # for 3rd_row elements

            if charge == 0:
                if sort_sym[0] == 'O' and sort_sym[1] == 'Br':  Mol_SO=-2.20 

                # COMMMENT: paper has it for cation, but it looks like it is for neutral
                if sort_sym[0] == 'Se' and sort_sym[1] == 'H':  Mol_SO=-4.21 

            if charge == +1:  ### RECHECK what the values of charge is!!!!!!!!!!!!!!!!!!!!!IMPORTANT
                if sort_sym[0] == 'K' and sort_sym[1] == 'Br':  Mol_SO=-2.99 
                if sort_sym[0] == 'H' and sort_sym[1] == 'As':  Mol_SO=-3.54 
                if sort_sym[0] == 'H' and sort_sym[1] == 'Br':  Mol_SO=-6.26 
                if sort_sym[0] == 'F' and sort_sym[1] == 'Br':  Mol_SO=-6.10 
                if sort_sym[0] == 'Na' and sort_sym[1] == 'Br':  Mol_SO=-3.93 
                if sort_sym[0] == 'Br' and sort_sym[1] == 'Br':  Mol_SO=-6.55 
        
        else:  # for non 3rd row elements, first and second rows
            if charge == 0:
                if sort_sym[0] == 'H' and sort_sym[1] == 'C':  Mol_SO=-0.07 
                if sort_sym[0] == 'O' and sort_sym[1] == 'H':  Mol_SO=-0.30 
                if sort_sym[0] == 'O' and sort_sym[1] == 'N':  Mol_SO=-0.27 
                if sort_sym[0] == 'O' and sort_sym[1] == 'Cl':  Mol_SO=-0.61 
                if sort_sym[0] == 'S' and sort_sym[1] == 'H':  Mol_SO=-1.01 
                if sort_sym[0] == 'P' and sort_sym[1] == 'O':  Mol_SO=-0.53 
                if sort_sym[0] == 'Si' and sort_sym[1] == 'H':  Mol_SO=-0.34 
            
            if charge == -1:
                if sort_sym[0] == 'N' and sort_sym[1] == 'H':  Mol_SO=-0.12 
                if sort_sym[0] == 'P' and sort_sym[1] == 'H':  Mol_SO=-0.45 
                if sort_sym[0] == 'O' and sort_sym[1] == 'O':  Mol_SO=-0.34 
                if sort_sym[0] == 'S' and sort_sym[1] == 'S':  Mol_SO=-1.12 
            
            if charge == +1:
                if sort_sym[0] == 'H' and sort_sym[1] == 'F':  Mol_SO=-0.62 
                if sort_sym[0] == 'P' and sort_sym[1] == 'H':  Mol_SO=-0.67 
                if sort_sym[0] == 'H' and sort_sym[1] == 'Cl':  Mol_SO=-1.60 
                if sort_sym[0] == 'N' and sort_sym[1] == 'N':  Mol_SO=-0.17 
                if sort_sym[0] == 'O' and sort_sym[1] == 'O':  Mol_SO=-0.43 
                if sort_sym[0] == 'P' and sort_sym[1] == 'P':  Mol_SO=-0.57 
                if sort_sym[0] == 'S' and sort_sym[1] == 'S':  Mol_SO=-1.25 
                if sort_sym[0] == 'Cl' and sort_sym[1] == 'Cl':  Mol_SO=-1.77 
                if sort_sym[0] == 'F' and sort_sym[1] == 'Cl':  Mol_SO=-1.60 

    return(Mol_SO)
###### Mol_SO - E

###### At_SO - S
def At_SO(at, charge):  # element symbol, charge
    At_SO = 0
    if charge == 0:
        S_dict_0 = {"B" : -0.05, "C" : -0.14, "O" : -0.36, "F" : -0.61, "Al" : -0.34, "Si" : -0.68, "S" : -0.89, "Cl" : -1.34, "Ga" : -2.51, "Ge" : -4.41, "Se" : -4.3,  \
            "Br" : -5.6, "Fe" : -1.84, "I" : -11.548}  
        if at in S_dict_0:  # this "if" executes only of the incoming atom "at" is present in the dict, else it returns 0 (this is equivalent to default in the f90 script)
            At_SO = S_dict_0[at]
            return(At_SO)
        else:
            return(At_SO)

    elif charge == +1:
        S_dict_1 = {"C"  :-0.2, "N"  :-0.43, "F"  :-0.67, "Ne" :-1.19, "Si" :-0.93, "P"  :-1.43, "Cl" :-1.68, "Ar" :-2.18, "Ge" :-5.37, "As" :-8.04, "Br" :-6.71, "Kr" :-8.16, "I"  :-14.028}
        if at in S_dict_1:  # this "if" executes only of the incoming atom "at" is present in the dict, else it returns 0 (this is equivalent to default in the f90 script)
            At_SO = S_dict_1[at]
            return(At_SO)
        else:
            return(At_SO)

    elif charge == -1:
        S_dict_n1 = {"B": -0.03, "O": -0.26, "Al": -0.28, "P": -0.45, "S": -0.88}
        if at in S_dict_n1:  # this "if" executes only of the incoming atom "at" is present in the dict, else it returns 0 (this is equivalent to default in the f90 script)
            At_SO = S_dict_n1[at]
            return(At_SO)
        else:
            return(At_SO)
    
    else:
        return(At_SO)
###### At_SO - E

###### NFC - S
def NFC(at):   # returns no. of forzen core electrons  # Why not use elif??? revisit
    if at in ['H','He']:
        return(0)
    if at in ['Li','Be','B','C','N','O','F','Ne','Na','Mg']:
        return(2)
    if at in ['Al','Si','P','S','Cl','Ar','K','Ca','Fe']:
        return(10)
    if at in ['Ga','Ge','As','Se','Br','Kr']:
        return(18)
    if at in ['I']:
        return(36)
###### NFC - E

###### sym2mass - S  # returns atomic masses ?
def sym2mass(at):
    mass_dict_1 ={"H" : 1.00783,"He" : 4.00260,"Li" : 7.01600,"Be" : 9.01218,"B" : 11.00931,"C" : 12.0000000,"N" : 14.00307,"O" : 15.99491,"F" : 18.99840,"Ne" : 19.99244,"Na" : 22.98977, \
        "Mg" : 23.98505,"Al" : 26.98154,"Si" : 27.97693,"P" : 30.97376,"S" : 31.97207,"Cl" : 34.96885,"Ar" : 39.96238,"K" : 38.96371,"Ca" : 39.96259,"Ga" : 68.92558,"Ge" : 73.92118, \
        "As" : 74.92160,"Se" : 79.91652,"Br" : 78.91834, "Kr" : 83.91151, "Fe" : 55.845, "I" : 126.90447}
    if at in mass_dict_1:
        sym2mass = mass_dict_1[at]
        return(sym2mass)
    else:   # this is equivalent of default in switch case, the write statements needs editing , revisit
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in sym2mass in module_geom.f90:  " + str(at)+ "\n")

###### sym2mass - E

###### nanb - S  # returns no. of alpha and beta electrons  ## CAUTION, RETURNS list, not 2 seperate numbers!!!
def nanb(at):
    na_dict_1 = {"H" : 1,"He" : 1,"Li" : 1,"Be" : 1,"B" : 2,"C" : 3,"N" : 4,"O" : 4,"F" : 4,"Ne" : 4,"Na" : 1,"Mg" : 1,"Al" : 2,"Si" : 3,"P" : 4,"S" : 4,"Cl" : 4,"Ar" : 4,"K" : 1,"Ca" : 1,\
        "Ga" : 2,"Ge" : 3,"As" : 4,"Se" : 4,"Br" : 4,"Kr" : 4,"Fe" : 10,"I" : 4}
    nb_dict_1 = {"H": 0,"He": 1,"Li": 0,"Be": 1,"B": 1,"C": 1,"N": 1,"O": 2,"F": 3,"Ne": 4,"Na": 0,"Mg": 1,"Al": 1,"Si": 1,"P": 1,"S": 2,"Cl": 3,"Ar": 4,"K": 0,"Ca": 1,"Ga": 1,"Ge": 1,"As": 1,\
        "Se": 2,"Br": 3,"Kr": 4,"Fe": 6,"I": 3,}
    if at in na_dict_1:
        na = na_dict_1[at]
        nb = nb_dict_1[at]
        return(na,nb)
    else:  # equivalent to case default, if at not in na_dict_1, this is printed, revist
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in nanb in module_geom.f90:  " + str(at) + " \n")

###### nanb - E


######## HOF_atoms - S
def HOF_atoms(at):
    h = 6.626070040*(10**-34)    #6.626070040d-34
    Ry = 10973731.568508  #10973731.568508d0
    c = 299792458    #299792458d0
    N_avo = 6.02214179000000*(10**+23)  #6.02214179000000d+23
    au2kcm = 2 * Ry * h * c * N_avo / 4184
    kcm2au = 1 / au2kcm

    hof_dict_1 = {"H": 51.63,"Li": 37.69,"Be": 76.48,"B": 136.2,"C": 169.98,"N": 112.53,"O": 58.99 ,"F": 18.47,"Na": 25.69,"Mg": 34.87,"Al": 78.23,"Si": 106.6,"P": 75.42,"S": 65.66,"Cl": 28.59, \
        "K": 21.48303059,"Ca": 42.38503824,"Fe": 98.7373,"Ga": 64.763384321,"Ge": 89.354,"As": 68.86,"Se": 57.8991,"Br": 28.1836,"I": 25.62}
#case ('K ') # todo collect reference for 3-rd row from Sambit
# case ('Br') !Phys. Chem. Chem. Phys., 2015, 17, 3584--3598
# case ('I') ! JANAF
# many other cases have numbers commented out...take a look at OG f90 file, thanks
    if at in hof_dict_1:
        HOF_atoms = hof_dict_1[at] * kcm2au
        return(HOF_atoms)
    else: # equivalent to case default
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in HOF_atoms: " + str(at)+ " \n")

######## HOF_atoms - E

######## dH_atoms - S
def dH_atoms(at):
    h = 6.626070040*(10**-34)    #6.626070040d-34
    Ry = 10973731.568508  #10973731.568508d0
    c = 299792458    #299792458d0
    N_avo = 6.02214179000000*(10**+23)  #6.02214179000000d+23
    au2kcm = 2 * Ry * h * c * N_avo / 4184
    kcm2au = 1 / au2kcm

# case ('Br')   !Phys. Chem. Chem. Phys., 2015, 17, 3584--3598
# case ('I')  !JANAF 
#case ('K ')  !1.4811185d0
# case ('Ca') !1.481118547d0
# case ('Ga') !1.5657266d0
    dH_dict_1 = {"H": 1.01,"Li": 1.10,"Be": 0.46,"B": 0.29,"C": 0.25,"N": 1.04,"O": 1.04,"F": 1.05,"Na": 1.54,"Mg": 1.19,"Al": 1.08,"Si": 0.76,"P": 1.28,"S": 1.05,"Cl": 1.10,"K": 1.6926,\
        "Ca": 1.3709,"Fe": 1.08,"Ga": 1.3291,"Ge": 1.104,"As": 1.23,"Se": 1.319,"Br": 2.930,"I": 1.58}
    if at in dH_dict_1:
        dH_atoms = dH_dict_1[at] * kcm2au
        return(dH_atoms)
    else:  # equivalent to case default
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in dH_atoms: " + str(at)+ " \n")

######## dH_atoms - E

######### atno - S
def atno(at):
    at_num_dict ={"H": 1,"He": 2,"Li": 3,"Be": 4,"B": 5,"C": 6,"N": 7,"O": 8,"F": 9,"Ne": 10,"Na": 11,"Mg": 12,"Al": 13,"Si": 14,"P": 15,"S": 16,"Cl": 17,"Ar": 18,"K": 19,"Ca": 20, \
        "Sc": 21,"Ti": 22,"V": 23,"Cr": 24,"Mn": 25,"Fe": 26,"Co": 27,"Ni": 28,"Cu": 29,"Zn": 30,"Ga": 31,"Ge": 32,"As": 33,"Se": 34,"Br": 35,"Kr": 36, "I":53}

    if at in at_num_dict:
        atno = at_num_dict[at]
        return(atno)
    else:
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in atno in module_geom.f90: " + str(at)+ " \n")

######### atno - E

######## uniqatoms - S
def uniqatoms(sym):  # to fine uniq atoms and the number of occurance of uniq atoms in the molecule
# we want 'ua' unique atom types
# we want 'uan' no. of unique atoms of each type
# Ex. C20H42
# N_ua = 2
# ua = 'C' 'H
# uan = 20  42
    uan = []
    uniq_sym = list(set(sym))  # ua
    N_ua = len(uniq_sym)
    for tmp_i in range(len(uniq_sym)):
        count_a = 0
        for tmp_j in range(len(sym)):
            if uniq_sym[tmp_i] == sym[tmp_j]:
                count_a = count_a + 1
        uan.append(count_a)
        uniqat_d = {}
        uniqat_d["N_ua"] = int(N_ua)
        uniqat_d["uniq_sym"] = uniq_sym
        uniqat_d["uan"] = uan
        print(uniqat_d)
    #return(N_ua, uniq_sym, uan)  ### revisit, how do you want variables returned?
    return(uniqat_d)
######## uniqatoms - E

######## printbas - S  ### new_file.com needs to be renamed
def printbas(fname, at, install_dir_val):  # fname = value of the  basis set name from the inp file,  install_dir_val = value from initial input file by user
    #print("===Entering printbas with: ", fname, at, " ===")
    basisSet_fpath = fname
    start_phrase = "NewGTO   "+ at
    num_lines_bas = sum(1 for line_tmp1 in open(basisSet_fpath,"r"))

    for temp_num, temp_l in enumerate(open(basisSet_fpath,"r")):
        if start_phrase in temp_l.strip():
            bas_start_lno = temp_num+1
            break

    #print(bas_start_lno)
    with open("input.com", "a") as new_f:

        for l1 in range(bas_start_lno,num_lines_bas):
            #print(l1)
            req_line_1 = linecache.getline(basisSet_fpath, l1)
            if "end" in req_line_1.strip():
                break
            else:
                new_f.write(req_line_1)

        new_f.write("   end\n")


    

##### IMPORTANT !!! where are the contents being written to????
######## printbas - E

######## principal_coord - S
def principal_coord():
    with open("geom.xyz","r") as xyz_f:
        num_l_xyz = sum(1 for l in xyz_f)
        Nat = int(linecache.getline("geom.xyz",1).strip())
        l1_chr_mul = linecache.getline("geom.xyz",2)
        sym = []   # will contain list of atom symbols in the mol # same order
        Rlist = []
        for tmp_l in range(3,num_l_xyz):   ### what does NH do ???
            l1_xyz_1 = linecache.getline("geom.xyz",tmp_l)
            l1_lsp1 = l1_xyz_1.split()
            sym.append(l1_lsp1[0])
            Rlist.append(float(l1_lsp1[1]))
            Rlist.append(float(l1_lsp1[2]))
            Rlist.append(float(l1_lsp1[3]))
        conv = np.array(Rlist)
        R_coord = np.resize(conv,[Nat,3])

        rCM_pre = []
        mass_list = []
        for tmp_j in range(Nat-1):
            print(sym[tmp_j])
            mass_atom = sym2mass(sym[tmp_j])
            mass_list.append(mass_atom)
            new1 = R_coord[tmp_j] * mass_atom
            rCM_pre.append(new1)
            #rCM_pre = R_coord[tmp_j][0] + rCM_pre[0])
        rCM_sum = np.sum(rCM_pre,axis=0)
        tot_mass = sum(mass_list)
        rCM = rCM_sum/tot_mass

        new_coord1 = []
        r_skew_symm = np.zeros((3,3))
        momin = 0.0
        for tmp_k in range(Nat-1):
            sub_cm = R_coord[tmp_k] - rCM
            new_coord1.append(sub_cm)
            r_skew_symm[0][0] = 0.0
            r_skew_symm[0][1] = -sub_cm[2]
            r_skew_symm[0][2] = sub_cm[1]

            r_skew_symm[1][0] = sub_cm[2]
            r_skew_symm[1][1] = 0.0
            r_skew_symm[1][2] = sub_cm[0]

            r_skew_symm[2][0] = sub_cm[1]
            r_skew_symm[2][1] = sub_cm[0]
            r_skew_symm[2][2] = 0.0

            momin = momin + (sym2mass(sym[tmp_k]) * np.matmul(np.transpose(r_skew_symm), r_skew_symm) )

        eig_val, eig_vec = np.linalg.eig(momin)
        Ievals = eig_val
        return(Ievals)


######## principal_coord - E

####### rung4mp2 - S
def rung4mp2(values, start_time_main):  
    with open("geom.xyz","r") as xyz_f:
        num_l_xyz = sum(1 for l in xyz_f)
        Nat = int(linecache.getline("geom.xyz",1).strip())
        l1_chr_mul = linecache.getline("geom.xyz",2)
        tria = l1_chr_mul.split()
        charge = tria[0]
        multip = tria[1]
        sym = []   # will contain list of atom symbols in the mol # same order
        Rlist = []
        sym_coord = []
        for tmp_l in range(3,num_l_xyz):   ### what does NH do ???
            l1_xyz_1 = linecache.getline("geom.xyz",tmp_l)
            l1_lsp1 = l1_xyz_1.split()
            full_l = l1_xyz_1.strip()
            sym_coord.append(full_l)  ## contains atom and coords
            sym.append(l1_lsp1[0])
            Rlist.append(l1_lsp1[1])
            Rlist.append(l1_lsp1[2])
            Rlist.append(l1_lsp1[3])
            conv = np.array(Rlist)
            R_coord = np.resize(conv,[Nat,3])


    if values["isatom"] != "true":
        if values["FROZEN_GEOM"] == "true":
### NEED TO STORE R, write principal coords here!!!!!
### NOT DONE : line 696 to 705

            with open("read_geom_freq.dat", "r") as froz_geom:
                num_l_fg = sum(1 for l in froz_geom)
                Nat = int(linecache.getline("read_geom_freq.dat",1).strip())
                l2_chr_mul = linecache.getline("read_geom_freq.dat",2)
                sym = []
                Rlist = []
                sym_coord = []
                for tmp_l in range(3,Nat+2):    ## check range
                    l1_xyz_1 = linecache.getline("read_geom_freq.dat",tmp_l)
                    l1_lsp1 = l1_xyz_1.split()
                    full_l = l1_xyz_1.strip()
                    sym_coord.append(full_l)  ## contains atom and coords
                    sym.append(l1_lsp1[0])
                    Rlist.append(l1_lsp1[1])
                    Rlist.append(l1_lsp1[2])
                    Rlist.append(l1_lsp1[3])
                    conv = np.array(Rlist)
                    R_coord = np.resize(conv,[Nat,3])
                nfreq =[]
                for tmp_j in range(Nat+3, num_l_fg):
                    tmp_line = float(linecache.getline("read_geom_freq.dat",tmp_j).strip())
                    nfreq.append(tmp_line)

## How does read_geom_freq.dat look like????
#       allocate(freq(1:nfreq), theta(1:nfreq))
#       do iat = 1, nfreq
#         !print *, iat
#         read(10,*)freq(iat)
#       enddo

            if Nat > 2:
                Ievals = principal_coord() # check
                linnonlin = "NL"
                if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                    linnonlin = "L"
            if Nat == 2:    linnonlin = "L"

            if linnonlin == "L":    nfreq = 3 * Nat - 5
            if linnonlin == "NL":   nfreq = 3 * Nat - 6


### revisit line 710 
        else:

    #call system_clock(count_start1)
            start_time1 = time.time()
            runorca(values["method_opt_freq"], values["basis_opt_freq"], "true", values["custombasis_opt_freq"], "false", values)
            end_time1 = time.time()
    #call system_clock(count_end1, count_rate1)
            if values["IPss"] != "true" or values["verticalIP"] != "true":
                with open("input.xyz", "r") as new_i_xyz:
                    num_l_in = sum(1 for l in new_i_xyz)
                    Nat = int(linecache.getline("input.xyz",1).strip())
                    title = linecache.getline("input.xyz",2)
                    sym = []
                    Rlist = []
                    sym_coord = []
                    for tmp_l in range(3,num_l_in+1):    ## check range
                        l1_xyz_1 = linecache.getline("input.xyz",tmp_l)
                        l1_lsp1 = l1_xyz_1.split()
                        full_l = l1_xyz_1.strip()
                        sym_coord.append(full_l)  ## contains atom and coords
                        sym.append(l1_lsp1[0])
                        Rlist.append(l1_lsp1[1])
                        Rlist.append(l1_lsp1[2])
                        Rlist.append(l1_lsp1[3])
                        conv = np.array(Rlist)
                        R_coord = np.resize(conv,[Nat,3])

############################################# REWRITE???
        if Nat > 2:
            Ievals = principal_coord() # check
            linnonlin = "NL"
            if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                linnonlin = "L"
        if Nat == 2:    linnonlin = "L"

############################################# REWRITE???
        
        os.system("nfreq=$( grep 'The total number of vibrations considered is ' input.out  | tail -1 | awk '{print $8}' > scr_lc0);  \
             nfreq=$( cat scr_lc0); grep -A$(( $nfreq+4 )) 'IR SPECTRUM' input.out | tail -$nfreq | awk '{print $2}' > freq.txt")

    ### need to extract nfreq from scr file
        tr = linecache.getline("scr_lc0",1).strip()
        lsp111 = tr.split()
        nfreq = int(lsp111[0])
#       open(unit=10, file='scr')
#       read(10,*) nfreq
#       close(10)
        if linnonlin == "NL":
            nfreq = (3*Nat)-6
            os.system("nfreq=$( cat scr_lc0); grep -A$(( $nfreq+4 )) 'IR SPECTRUM' input.out | tail -$(( "+str(nfreq)+" )) | awk '{print $2}' > freq.txt")
            # write(*,'(a)') trim(cmd)  - IMPORTANTT!!!!!!!!!!!!!!
## revist ...write statement might be missing
        freq =[]
        with open("freq.txt", "r") as fr:
            num_l_frq = sum(1 for l in fr)
            for tmp_p in range(1,num_l_frq+1):  # check range
                tmp_val = float(linecache.getline("freq.txt",tmp_p).strip())
                freq.append(tmp_val)


        iopt = 999
        if values["isatom"] != "true":
            os.system("grep 'OPTIMIZATION RUN DONE' input.out | wc | awk '{print $1}' > scr_lc1")
        iopt = int(linecache.getline("scr_lc1",1).strip())
 
        os.system("rm -f input* scr freq.txt")
        #os.system("rm -f scr ")


    ################################################################################################

    if linnonlin == "NL":
        if nfreq != ((3*Nat)-6):
            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write(" ************************************\n")
                ther_chem.write(" * WARNING: No. of freq. != 3N - 6 *\n")  # WARNING!!
                ther_chem.write(" ************************************\n")
   # MANUAL
    #freq = [30.0,324.0, 3433, 4, 4, 4, 5, 5, 6,7 , 8, 8,30.0,324.0, 3433, 4, 4, 4, 5, 5, 6,7 , 8, 8, 1, 2,3] 
    os.system("rm -f freq0.txt")
    with open("freq0.txt", "a") as new_freq0:
        for tmp_i in range(len(freq)):
            new_freq0.write(str(freq[tmp_i]) + "\n")  # check if this is printed in proper format


    if values["isatom"] != "true":
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write(" Optimized atomic coordinates (Angstrom)\n")
            for tmp in range(Nat):  ## CHECK range
                ther_chem.write(str(sym_coord[tmp]) + "\n")  # writes atoms + coord as from geo,xyz
            ther_chem.write("\n")

            ther_chem.write(" Unscaled harmonic wavenumbers (cm^-1)\n")
            tmp = 0
            print(nfreq)
            for tmp in range(nfreq): # check range
                #print(type(values["scalfac"]))# * values["scalfac"])
                ther_chem.write(str( freq[tmp] ))  
                ther_chem.write("\n")

            ther_chem.write("\n")
            ther_chem.write(" Scaled harmonic wavenumbers (cm^-1)\n")
            tmp = 0
            print(nfreq)
            for tmp in range(nfreq): # check range
                #print(type(values["scalfac"]))# * values["scalfac"])
                ther_chem.write(str( freq[tmp] * float(values["scalfac"]) ) )  
                ther_chem.write("\n")

            ther_chem.write("\n")
            ther_chem.write(" Scaling factor used:  " + str(values["scalfac"]) + "\n\n")
            #call system_clock(count_end, count_rate)
            end_time2 = time.time()
            ther_chem.write(" * Geometry optimization/frequencies done in " + str(end_time1 - start_time1) + " s\n")
            ther_chem.write(" ** Elapsed time =  " + str(end_time2 -start_time_main) + " s\n")
    else:
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write(" * Geometry optimization skipped for atom\n\n")
    
    if values["isatom"] != "true":
        if values["FROZEN_GEOM"] != "true":
            if iopt == 0:
                with open("Thermochemistry.out", "a") as ther_chem:
                    ther_chem.write("E R R O R: Geometry Optimization Failed")

# !=== SINGLE POINT CCSD(T)
    if values["DLPNO_CCSDT"] == "true":
        values["switch_DLPNO_CCSDT"] = "false"
    
    start_time1 = time.time()
    if values["Ntotale"] > 0:
        if Nat == 1:
            if values["basis_ccsdt"]  == "GTBAS1": values["basis_ccsdt"] = 'GTBAS1atm'
        if values["restart_cc"] == "true":
            charge = charge + 2
            tmp_method = "HF "
            runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
            charge = charge - 2
    
        if values["DLPNO_CCSDT"] == "true":
            values["switch_DLPNO_CCSDT"] == "true"
        
        runorca(values["method_ccsdt"], values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
        end_time1= time.time()
        
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc2")
        else:
            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc2")
        
        tmp = linecache.getline("scr_lc2",1).strip()
        E_ccsdt = float(tmp)
        print("== CCSD(T) Energy ==", E_ccsdt)
        
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            E_mp2 = E_ccsdt
        else:
            os.system("grep 'Initial E(tot)' input.out | awk '{print $4}' > scr_lc20")
            tmp = linecache.getline("scr_lc20",1).strip()
            E_mp2 = float(tmp)
            os.system("rm -f input* scr ")

    ### IMPORTANT --- indentation ERROR!! revisit
    else:
        E_ccsdt = 0
        E_mp2 = 0

    if values["G4MP2TM"] == "true":
        values["switch_load_rel_file"] = "true"
        values["switch_guess"] = "true"
        values["ALLELE"] = "true"

    ## POSSIBLE indentation faults from f90 ...ask to recheck , revisit
        if values["Ntotale"] > 0:
            if (values["Notale"] < 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1 ):
                tmp_method = "HF "
                if Nat == 1:
                    if values["basis_ccsdt"] == "GTBAS1": values["basis_ccsdt"] = "GTBAS1atm"
                    runorca(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                else:  # ASK
                    runorca(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
            elif (values["Ntotale"] == 2) or (values["Ntotale"] - values["Ntotalecore"] == 1 ):
                if Nat == 1:
                    if values["basis_ccsdt"] == "GTBAS1": 
                        values["basis_ccsdt"] = "GTBAS1atm"  # check
                        runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                    else:
                        runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
            else:
                if values["restart_cc"] == "true":
                    charge = charge + 1
                    tmp_method = "HF "
                    runorca(tmp_method,values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                    charge = charge - 2
                runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
            end_time1 = time.time()

            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc3")
            tmp = linecache.getline("scr_lc3",1).strip()
            E_ccsdt_rel = float(tmp)
    # indent issues IMPORTANT !!!!!!! revisit
        else:
            E_ccsdt_rel = 0

        values["switch_load_rel_file"] = "false"
        values["switch_guess"] = "false"
        values["ALLELE"] = "false"
            
        E_rel = E_ccsdt_rel - E_ccsdt

    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * CCSD(T) done in " + str(end_time1 - start_time1) + " s\n")
        ther_chem.write(" ** Elapsed time = " + str(end_time_n - start_time_main) + " s \n\n")
        
    if values["switch_read_RIMP2_small"] != "true":
        os.system("rm -f input* scr ")
    
    if values["DLPNO_CCSDT"] == "true":
        values["switch_DLPNO_CCSDT"] == "false"

## == TESTING ( comm from f90)
    if values["DLPNO_CCSDT"] == "true":
        if values["switch_read_RIMP2_small"] != "true":
            if (multip != 1) or (values["switch_RIMP2_small"] == "true") :
                start_time1 = time.time()
                if (values["Ntotale"] == "1") or (values["Ntotale"] - values["Ntotalecore"] <= 1) :
                    tmp_method = "HF "
                    runorca(values["method_mp2_s"], values["basis_mp2_s"], "false", values["custombasis_mp2_s"], "true", values)
                else:
                    runorca(values["method_mp2_s"], values["basis_mp2_s"], "false", values["custombasis_mp2_s"], "true", values)
                end_time1= time.time()
## why not club thse two???????/
                if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1) :
                    os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc4")
                else:
                    os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc4")
            tmp = linecache.getline("scr_lc4",1).strip()
            E_mp2 = float(tmp)
            os.system("rm -f scr ")

            end_time_n = time.time()
            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write(" * MP2/S done (for non-singlets) in   " + str(end_time1 - start_time1) + " s \n")
                ther_chem.write(" ** Elapsed time = " + str(end_time_n - start_time_main) + " s \n\n")
#MP2 TOTAL ENERGY
    if values["switch_read_RIMP2_small"] == "true":
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc5")
        else:
            os.system("grep 'MP2 TOTAL ENERGY' input.out  | awk '{print $4}' > scr_lc5")
        tmp = linecache.getline("scr_lc5",1).strip()
        E_mp2 = float(tmp)
    os.system("rm -f input* scr ")

# === SINGLE POINT MP2L
    start_time1 = time.time()  # can this be moved into the if statement?
    if (values["Ntotale"] > 0) or (values["ccsdt_cbs"] != "true") : ## where is ccsdt_cbs declared?
        if values["restart_mp2"] == "true":
            charge = charge + 2
            tmp_method = "HF "
            runorca(tmp_method, values["basis_mp2"], "false", values["custombasis_mp2"], "true", values)
            charge = charge - 2
        runorca(values["method_mp2"], values["basis_mp2"], "false", values["custombasis_mp2"], "true", values)
        end_time1= time.time()
    
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc6")
        else:
            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc6")
        tmp = linecache.getline("scr_lc6",1).strip()
        E_mp2L = float(tmp)

        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            E_hfL = E_mp2L
        else:
            if values["flag_RIMP2"] == "true":
                os.system("grep 'RI-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            elif values["flag_DLPNOMP2"] == "true":
                os.system("grep 'DLPNO-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            else:
                os.system("grep 'MP2 CORRELATION ENERGY ' input.out | awk '{print $5}' > scr_lc7")
        tmp = linecache.getline("scr_lc7",1).strip()
        E_hfL = float(tmp)
        E_hfL = E_mp2L - E_hfL

        os.system("rm -f input* scr ")
# VERIFY - indentation issues, which else this goes to??? revisit
    else:
        E_hfL= 0
        E_mp2L = 0

    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * MP2/L done in   " + str(end_time1 - start_time1) + " s \n")
        ther_chem.write(" ** Elapsed time = " + str(end_time_n - start_time_main) + " s \n\n")

#== SINGLE POINT HF/VTZ
    start_time1 = time.time()
    if (values["Ntotale"] > 0) and (values["ccsdt_cbs"] != "true"):
        if values["restart_hf3"] == "true":
            charge = charge + 2
            tmp_method = "HF "
            runorca(tmp_method, values["basis_hf3"], "false", values["custombasis_hf3"], "false", values)
            charge = charge - 2
        runorca(values["method_hf3"], values["basis_hf3"], "false", values["custombasis_hf3"], "false", values)
        end_time1 = time.time()

        os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc8")
        tmp = linecache.getline("scr_lc8",1).strip()
        E_hfT = float(tmp)
        
        os.system("rm -f input* scr ")
    else:
        E_hfT = 0
    
    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * HF/VTZ done in  " + str(end_time1 - start_time1) + " s\n")
        ther_chem.write(" ** Elapsed time = " + str(end_time_n - start_time_main) + " s\n\n")

# == SINGLE POINT HF/VQZ
    start_time1 = time.time()
    if (values["Ntotale"] > 0 ) and (values["ccsdt_cbs"] != "true"):
        if values["restart_hf4"] == "true":
            charge = charge + 2
            tmp_method = "HF "
            runorca(tmp_method, values["basis_hf4"], "false", values["custombasis_hf4"], "false", values)
            charge = charge - 2
        runorca(values["method_hf4"], values["basis_hf4"], "false", values["custombasis_hf4"], "false", values)
        end_time1 = time.time()

        os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc9")
        tmp = linecache.getline("scr_lc9",1).strip()
        E_hfQ = float(tmp)
        
        os.system("rm -f input* scr ")
    else:
        E_hfQ = 0
    
    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * HF/VQZ done in  " + str(end_time1 - start_time1) + " s\n")
        ther_chem.write(" ** Elapsed time = " + str(end_time_n  - start_time_main) + " s\n\n")

# Thermal
    kB = 1.38064852000000*(10**-23) #kB = 1.38064852000000d-23
    h = 6.626070040*(10**-34)    #6.626070040d-34
    Ry = 10973731.568508  #10973731.568508d0
    c = 299792458.0    #299792458d0
    T = 298.15  #298.15d0 
    e = 1.60217648700000*(10**-19)   #1.60217648700000d-19
    N_avo = 6.02214179000000*(10**+23)  #6.02214179000000d+23
    au2j = 2 * Ry * h * c

    au2ev  = 2 * Ry * h * c          / e
    au2kjm = 2 * Ry * h * c * N_avo / 1000.0
    au2kcm = 2 * Ry * h * c * N_avo / 4184.0
    kcm2au = 1 / au2kcm

    j2au   = 1 / au2j
    au2cmi = 2 * Ry / 100.0
    cmi2au = 1 / au2cmi
    Theta = []
    if values["isatom"] != "true":
        for tmp_i in range(len(freq)):
            freq[tmp_i] = freq[tmp_i] * float(values["scalfac"])
            Theta_new = c * 100.0 *  h * freq[tmp_i] / kB
            Theta.append(Theta_new)
        dE_vib = 0
        for tmp in range(nfreq):  # check range
            dE_vib = dE_vib + kB * Theta[tmp] * j2au * ( 0.5 + 1 / ( math.exp(Theta[tmp]/T) - 1) )  ## check
        ZPVE = 0.5 * sum(freq) * cmi2au
        if linnonlin == "L": 
            dE_rot = kB * T * j2au
        elif linnonlin == "NL":
            dE_rot = kB * T * j2au * (3.0/2.0)
        else:
            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write("ERROR: unknown entry for linear/nonlinear")
    else:
        ZPVE    = 0
        dE_vib  = 0
        dE_rot  = 0

    dE_tra = kB * T * j2au * (3.0/2.0)
    Ethermal = dE_vib + dE_rot + dE_tra

#  !=== SO
    HLC_SO = 0
    if values["isatom"] != "true":
        HLC_SO = Mol_SO(Nat, multip, charge, sym, values["SO_3rdrow_mols"])
        print(HLC_SO)
    else:
        HLC_SO = At_SO(sym[0], charge)
    
    if values["HLCeqZERO"] == "true":
        HLC = 0
    else:
        if values["isatom"] != "true":
            if values["Na"] == values["Nb"]:
                HLC = - values["AA"] * values["Nb"]
            else:
                HLC = - ( values["ApAp"] * values["Nb"] ) - (values["BB"] * (values["Na"]- values["Nb"]) )  # check
        else:
             HLC = - (values["CC"] * values["Nb"] ) - (values["DD"] * (values["Na"] - values["Nb"]) )  # + At_SO(sym(1), charge)
        
        if values["isatom"] == "true":
            if sym[0] == "Be" and charge == 0 : HLC == - values["EE"]
            if sym[0] == "Mg" and charge == 0 : HLC == - values["EE"]
            if sym[0] == "Ca" and charge == 0 : HLC == - values["EE"]
            if sym[0] == "Li" and charge == -1 : HLC == - values["EE"]
            if sym[0] == "Na" and charge == -1 : HLC == - values["EE"]
            if sym[0] == "K" and charge == -1 : HLC == - values["EE"]
        
        if Nat == 2:
            if charge == 0:
                if sym[0] == "Li" and sym[1] == "Li": HLC = - values["EE"]
                if sym[0] == "Na" and sym[1] == "Na": HLC = - values["EE"]
                if sym[0] == "K" and sym[1] == "K": HLC   = - values["EE"]
                if sym[0] == "Li" and sym[1] == "Na": HLC = - values["EE"]
                if sym[0] == "Na" and sym[1] == "Li": HLC = - values["EE"]
                # if ( (trim(sym(1)) .eq. 'Be') .and. (trim(sym(2)) .eq. 'H' ) ) HLC = -EE
                # if ( (trim(sym(1)) .eq. 'H' ) .and. (trim(sym(2)) .eq. 'Be') ) HLC = -EE
    print(values["AA"],values["ApAp"],values["BB"],values["CC"],values["DD"],values["EE"],values["Na"],values["Nb"])
    HLC  = HLC / 1000.0
    HLC_SO = HLC_SO / 1000.0
    HLC0 = HLC

# CBS extrapolation hard-coded for HF (AVQZ, AV5Z)
    if values["HF_CBS_default"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-1.63)) / (1 - math.exp(-1.63))
    if values["HF_CBS_orca_23_def2"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-10.39*(math.sqrt(3.0)- math.sqrt(2.0)))) / (1 - math.exp(-10.39*(math.sqrt(3.0)-math.sqrt(2.0))))
    if values["HF_CBS_orca_34_def2"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-7.88*( math.sqrt(4.0)- math.sqrt(3.0)))) /  (1 - math.exp(-7.88*(math.sqrt(4.0)-math.sqrt(3.0))))
    if values["HF_CBS_orca_23_cc"] == "true":   E_HF_limit = (E_hfQ - E_hfT * math.exp(-4.42*( math.sqrt(3.0)- math.sqrt(2.0)))) / (1 - math.exp(-4.42* (math.sqrt(3.0)-math.sqrt(2.0))))
    if values["HF_CBS_orca_34_cc"] == "true":   E_HF_limit = (E_hfQ - E_hfT * math.exp(-5.46*( math.sqrt(4.0)- math.sqrt(3.0)))) / (1 - math.exp(-5.46* (math.sqrt(4.0)-math.sqrt(3.0))))
    
    E_HF_CBS   = E_HF_limit - E_hfL
    print(E_mp2L,E_mp2)
    E_MP2_CBS  = E_mp2L - E_mp2
    if values["noHLC"] == "true": HLC = 0
    
    if values["ccsdt_cbs"] == "true":
        E_HF_CBS = 0
        E_MP2_CBS = 0
        eps0 = E_ccsdt +  HLC + HLC_SO
    else:
        eps0 = E_ccsdt + E_HF_CBS + E_MP2_CBS + HLC + HLC_SO
    
    if values["G4MP2TM"] == "true":
        eps0 = eps0 + E_rel

    U0 = eps0 + ZPVE
    UT = eps0 + Ethermal
    HT = eps0 + Ethermal + kB * T * j2au

    if values["IPss"] == "true":   ## post check reduce all indices by 1
        IPbreakdown = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC,HLC_SO,ZPVE  ]
        EAbreakdown = [ E_ccsdt, E_HF_CBS, E_MP2_CBS,HLC,HLC_SO,ZPVE  ]
        PAbreakdown = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE, (dE_vib - ZPVE) , dE_rot, dE_tra , (kB * T * j2au) ]
    else:
        IPbreakdown = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE ]
        EAbreakdown = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE ]
        PAbreakdown = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE, (dE_vib - ZPVE), dE_rot, dE_tra, (kB * T * j2au)  ]

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write("Temperature=  " + str(T) +  "   Pressure=  1.0  \n")
        ther_chem.write("E(ZPE)=  " +str(ZPVE) + "  E(Thermal)=  " + str(Ethermal) + "\n")  # revisit
        if values["G4MP2TM"] == "true":
            ther_chem.write("E(CCSD(T),rel)=  "  + str(E_ccsdt_rel) + "\n")
        ther_chem.write("HLC=  " + str(HLC) + "   SO=  " + str(HLC_SO) + "  \n")
        ther_chem.write("E(CCSD(T))=  " + str(E_ccsdt)  + "\n")
        ther_chem.write("DE(MP2)=   "  + str(E_MP2_CBS) + "  DE(HF)=  " + str(E_HF_CBS)  + "\n")
        ther_chem.write("G4MP2(0 K)=  "  + str(U0)  + " G4MP2 Energy= " + str(UT) + "\n")
        ther_chem.write("G4MP2 Enthalpy=  " + str(HT) + "  G4MP2 Free Energy= 0.0 \n\n")  # verify
        
        end_time_n = time.time()
        ther_chem.write(" * G4(MP2) done\n")
        ther_chem.write(" ** Elapsed time = "  + str(end_time_n - start_time_main)  + "  s \n\n")
        
    os.system("rm -f scr scrd")

####### rung4mp2 - E


####### runorca - S
def runorca(method, basis,optfreq,custombasis, correlated, values):
    #print("=== Entering runorca ===")
    with open("input.com", "w") as com_f:
        if optfreq == "true":
            if values["verticalIP"] != "true" or values["IPss"] != "true":   # IPss not defined
                if values["MGGA"] == "true":
                    Freqstr="NumFreq"
                else:
                    Freqstr="Freq"
        
                if custombasis == "true":
                    com_f.write("! " +str(method) + "  " + values["String_Opt"] + "  " + Freqstr + "  \n")
                else:
                    com_f.write("! " +str(method) + "  " + str(basis) +"   "+values["String_Opt"] + "  " + Freqstr + "  \n")
            else:
                if custombasis == "true":
                    com_f.write("! " +str(method) + "  " + Freqstr + "  \n")
                else:
                    com_f.write("! " +str(method) + "  " + str(basis) + "  " + Freqstr + "  \n")
        else:
            if custombasis == "true":
                com_f.write("! " +str(method) + "  \n")
            else:
                com_f.write("! " +str(method) + "  " + str(basis) + "  \n")

        with open("geom.xyz","r") as xyz_f:
            num_l_xyz = sum(1 for l in xyz_f)
            l1_chr_mul = linecache.getline("geom.xyz",2)
            xyz_l1 = "* xyz  " + l1_chr_mul.strip() + " \n"
            com_f.write(xyz_l1)
            sym = []   # will contain list of atom symbols in the mol # same order
            for tmp_l in range(3,num_l_xyz+1):   ### what does NH do ???
                l1_xyz_1 = linecache.getline("geom.xyz",tmp_l)
                l1_lsp1 = l1_xyz_1.split()
                sym.append(l1_lsp1[0])
                com_f.write(l1_xyz_1)
        com_f.write("*\n")
        com_f.write("%MaxCore  " + values["maxcore_mb"] + "\n")
        com_f.write("%scf\n   MaxIter 500 \n")
        com_f.write("   Convergence  " + values["conv_scf"] + "\n")
        com_f.write("end\n")

        if values["switch_guess"] == "true":   ### this is not part of the inp file!!!
            if values["guess_TM"] == "true" and values["G4MP2TM"]:
                com_f.write("   Guess = " + values["option_guess"] + "\n")
                com_f.write("end\n")
        
        if values["switch_load_rel_file"] == "true":
            ### line 1499  IMPORTANT !!!!!!!!!!!!!
            f1 = open("rel_file.txt", "r")
            com_f.write(f1.read())
            f1.close()
            print("check if rel_file.txt exists!!")
            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write("check if rel_file.txt exists!!")
# ! if ( G4MP2TM .and. load_rel_file) then
#! %basis
#! NewGTO I "old-DKH-SVP" end
#! end
        if values["SCFDIIS"] == "true":
            com_f.write("%scf\n  DIISMaxEq 15\n")
            com_f.write("  directresetfreq 1\n")
            com_f.write("end\n")

        if values["LSHIFT"] == "true":
            com_f.write("%scf\n")
            com_f.write("  Shift Shift 0.1 ErrOff 0.1 end\n")
            com_f.write("end\n")

        if values["SOSCF"] == "true":
            com_f.write("%scf\n")
            com_f.write("  soscfmaxit 12\n")
            com_f.write("  directresetfreq 1\n")
            com_f.write("end\n")
        if values["switch_DLPNO_CCSDT"] == "true":
          com_f.write("%mdci\n")
          com_f.write("  UseFullLMP2Guess true\n")
          com_f.write("  TcutDOPre =  " + str(values["TcutDOPre"]) +"\n")
          com_f.write("end\n")
        ### revisit for datatype check!!!!!! IMPORTANT
        if ( float(values["Ntotale"]) <= float(values["nproc"]) ) or ( (float(values["Ntotale"])-float(values["Ntotalecore"]))  < float(values["nproc"]) ):
            com_f.write("%pal nprocs  1 \n")
        else:
            com_f.write("%pal nprocs  "+values["nproc"]+" \n")
        
        com_f.write("end\n")
        com_f.write("%method\n")  ## CHECK
        com_f.write("  IntAcc 7.0\n")
        
        if values["optdiis"] == "true":
            com_f.write("  Z_solver DIIS\n")
            com_f.write("  Z_MaxIter 300\n")
        
        if correlated == "true":
            if values["ALLELE"] == "true":  ### CHECK!!!!
                uniq_atom_res = uniqatoms(sym)
                for iat in range(int(uniq_atom_res["N_ua"])):
                    pre1 = uniq_atom_res["uniq_sym"]
                    at_pr1 = pre1[iat]
                    com_f.write("  NewNCore " + at_pr1 + "  " + " 0  end\n")
                for iat in range(int(uniq_atom_res["N_ua"])):
                    pre1 = uniq_atom_res["uniq_sym"]
                    at_pr1 = pre1[iat]
                    NFC_res = NFC(at_pr1)
                    com_f.write("  NewNCore " + at_pr1 + "  " + NFC_res +"  end\n")
        
        com_f.write("end\n")
        
        if optfreq == "true":
            com_f.write("%geom\n")
            if values["MGGA"] == "true":
                com_f.write("  Calc_Hess true; NumHess true\n")
            else:
                com_f.write("  Calc_Hess true\n")
            com_f.write("  Recalc_Hess " + str(values["iterhess"]) +" \n")  ## revisit !!!! CHECK!!! IMPORTANT
            com_f.write("end\n")
            com_f.write("%freq  Temp  273.15, 298.15\n")
            com_f.write("end\n")
            #com_f.write("end")
        if custombasis == "true":
            com_f.write("%basis \n")
        
    if custombasis == "true":
        with open("input.com", "a") as com_f:
            uniq_atom_res = uniqatoms(sym)
            fname = values["install_dir"] + basis
            #print("fname is ", fname)
            ##write(fname,'(a,a)')trim(adjustl(install_dir)), trim(adjustl(basis))  ## IMPORTANT
            #print(uniq_atom_res)
            if Nat == 1:  # check datatype
                printbas(fname, sym[0], values["install_dir"])    ## what is fname, it should be the basis set keyword. but which one??
            else:
                for iat1 in range(int(uniq_atom_res["N_ua"])):
                    pre2 = uniq_atom_res["uniq_sym"]
                    at_pr2 = pre2[iat1]
                    #print(iat1, at_pr2)
                    printbas(fname, at_pr2, values["install_dir"])
            com_f.write("end\n")
#   os.system("cat basis.com >> input.com")

    print(values["orca_dir"] + "orca input.com > input.out")
    os.system(values["orca_dir"] + "orca input.com > input.out")
    os.system("cat input.com >> ORCA_G4MP2.com")
    os.system("cat input.out >> ORCA_G4MP2.out")
####### runorca - E


### Reading input settings from input_f - end
start_time_main = time.time()
### Writing intro Thermochemistry.out - start
with open("Thermochemistry.out", "w") as ther_chem:
    day = datetime.datetime.today().strftime("%A")
    date = datetime.date.today()
    cur_date = date.strftime("%B %d, %Y")
    now = datetime.datetime.now()
    cur_time = now.strftime(" %H:%M:%S")


    ther_chem.write("==============================================\n")
    ther_chem.write("|           G4MP2 THERMOCHEMISTRY            |\n")
    ther_chem.write("==============================================\n")
    ther_chem.write("Starting time:  " + day + "  "+ cur_date + "  " + cur_time + "  " + time.tzname[0] +"  \n\n")
    ther_chem.write("* Constants and conversion factors *\n\n")
    ther_chem.write("h      =   6.6260700399999999E-34 J s\n")
    ther_chem.write("kB     =   1.3806485199999999E-23 J K^-1\n")
    ther_chem.write("Ry     =   1.0973731568507999E+07 m^-1\n")
    ther_chem.write("e      =   1.6021764870000000E-19 C\n")
    ther_chem.write("c      =   2.9979245800000000E+08 m s^-1\n")
    ther_chem.write("N_avo  =   6.0221417899999999E+23 mol^-1\n\n")

    ther_chem.write("au2ev  =   2.7211388291762535E+01\n")
    ther_chem.write("au2kcm =   6.2750957099203276E+02\n")
    ther_chem.write("au2kjm =   2.6255000450306652E+03\n")
    ther_chem.write("au2cmi =   2.1947463137015997E+05\n\n")
    
    ther_chem.write(" * Parameterized for the following atoms  *\n")
    ther_chem.write("    ____                                    \n")
    ther_chem.write("    |H |                                    \n")
    ther_chem.write("    --------    ---------------------       \n")
    ther_chem.write("    |Li| Be|    |B  |C  |N  |O  |F  |       \n")
    ther_chem.write("    |Na| Mg|    |Al |Si |P  |S  |Cl |       \n")
    ther_chem.write("    |Na| Mg|    |Al |Si |P  |S  |Cl |       \n")
    ther_chem.write("    |K | Ca|    |Ga |Ge |As |Se |Br |       \n")
    ther_chem.write("    --------    ---------------------       \n")
    ther_chem.write("\n")
### Writing intro Thermochemistry.out - end

val = {}
### List of predefined variables ?? does the order and manner of occurance matter or can all of them be declared at the begining?
val["ALLELE"] = "false"
val["switch_load_rel_file"] = "false"
val["switch_guess"] = "false"
val["FROZEN_GEOM"] = "false"


val["SO_3rdrow_mols"] = "false"
val["iterhess"] = 5

val["AA"] = 9.472
val["BB"] = 3.102
val["CC"] = 9.741
val["DD"] = 2.115
val["ApAp"] = 9.769
val["EE"] = 2.379

val["TcutDOPre"] = 3e-2

val["HF_CBS_default"]      = "true"
val["HF_CBS_orca_23_def2"] = "false" 
val["HF_CBS_orca_34_def2"] = "false" 
val["HF_CBS_orca_23_cc"]   = "false" 
val["HF_CBS_orca_34_cc"]   = "false" 

val["switch_DLPNO_CCSDT"]  = "false"

val["G4MP2TM"] = "false"
val["ccsdt_cbs"] = "false"
val["noHLC"] = "false"


#  open(unit=999, file='orca_g4mp2.inp')
#  read(999,nml=ORCA_G4MP2_INPUT) 
# not written - needed? revisit

### Reading input settings from input_f - start
inp_par = ["openmpi_dir","orca_dir","install_dir","maxcore_mb","nproc","method_opt_freq", "basis_opt_freq", "custombasis_opt_freq", "String_Opt", "MGGA", "FROZEN_GEOM", "method_ccsdt", \
    "basis_ccsdt", "custombasis_ccsdt", "DLPNO_CCSDT", "switch_RIMP2_small", "switch_read_RIMP2_small", "method_mp2_s", "basis_mp2_s", "custombasis_mp2_s", "method_mp2", "basis_mp2", \
    "custombasis_mp2", "flag_RIMP2", "flag_DLPNOMP2", "method_hf3", "basis_hf3", "custombasis_hf3", "method_hf4", "basis_hf4", "custombasis_hf4", "scalfac", "calc_IP", "verticalIP", "calc_EA", \
    "verticalEA", "calc_PA", "calc_AE", "calc_BE", "calc_HF", "conv_scf", "HLCeqZERO", "SOSCF", "SCFDIIS", "SO_3rdrow_mols", "LSHIFT", "optdiis", "HF_CBS_default", "HF_CBS_orca_23_def2", \
    "HF_CBS_orca_34_def2", "HF_CBS_orca_23_cc", "HF_CBS_orca_34_cc", "restart_cc", "restart_mp2", "restart_hf3", "restart_hf4"]

with open(input_f, "r") as ini_inp:
    u = 0
    while u < len((inp_par)):
        for line in ini_inp:
            if inp_par[u] in line:
                val_line = line.split("=")
                if "#" in val_line[1]: ### CHECK!!!!!!
                    val_split = val_line[1].split("#")   ## this condition is to make sure if some comments are added to the input line, the comment (that follows # ) is ignored
                    val[inp_par[u]] = val_split[0].strip()
                else:
                    val[inp_par[u]] = val_line[1].strip()
                break
        u = u +1


with open("geom.xyz","r") as xyz_f:
    num_l_xyz = sum(1 for l in xyz_f)
    Nat = int(linecache.getline("geom.xyz",1).strip())
    l1_chr_mul = linecache.getline("geom.xyz",2)
    tria = l1_chr_mul.split()
    charge = int(tria[0])
    multip = int(tria[1])
    sym = []   # will contain list of atom symbols in the mol # same order
    Rlist = []
    sym_coord = []
    for tmp_l in range(3,num_l_xyz+1):   ### what does NH do ???
        l1_xyz_1 = linecache.getline("geom.xyz",tmp_l)
        l1_lsp1 = l1_xyz_1.split()
        full_l = l1_xyz_1.strip()
        sym_coord.append(full_l)  ## contains atom and coords
        sym.append(l1_lsp1[0])
        Rlist.append(l1_lsp1[1])
        Rlist.append(l1_lsp1[2])
        Rlist.append(l1_lsp1[3])
        conv = np.array(Rlist)
        R_coord = np.resize(conv,[Nat,3])

val["Ntotal"] = 0
val["Ntotale"] = 0
val["Ntotalecore"] = 0
for tmp_j in range(Nat):
    na_nb_list = nanb(sym[tmp_j])
    a = na_nb_list[0]
    b = na_nb_list[1]
    val["Ntotal"] = val["Ntotal"] + a + b
    print(tmp_j,a,b)
    val["Ntotale"] = val["Ntotale"] + atno(sym[tmp_j])
    val["Ntotalecore"] = val["Ntotalecore"] + NFC(sym[tmp_j])


val["Ntotal"] = val["Ntotal"] - charge
val["Ntotale"] = val["Ntotale"] - charge
val["Na"] = (multip + val["Ntotal"] - 1 ) / 2.0
val["Nb"] = val["Na"] + 1 - multip

print(val["Na"],val["Nb"],val["Ntotal"])

val["isatom"] = "false"
if Nat == 1: val["isatom"] == "true"

uniqat_d = uniqatoms(sym)

val["IPss"] = "false"
rung4mp2(val, start_time_main)
os.system("cp freq0.txt freq.txt")

# revisit
U0mol = U0
UTmol = UT
HTmol = HT
HLCmol = HLC0

 #=== JOB-TYPE-1 IS IP, from here
if val["calc_IP"] == "true":
    val["IPss"] == "true"
    with open("geom_cation_IP.xyz", "r") as g_cat:  # need to look at this file
        with open("input.xyz", "w") as inp_x:
            Nat= int(linecache.getline("geom_cation_IP.xyz",1).strip())
            tr = linecache.getline("geom_cation_IP.xyz",2).strip()
            lsp11 = tr.split()
            charge = lsp11[0]
            multip = lsp11[1]
            inp_x.write(str(Nat) + " \n")
            inp_x.write("0  " + str(multip) + " \n")
            #inp_x.write(sym_coord)  # check  what coords are written here? from geom.xyz file or from geom_cation_IP files??? -- this line writes from geom.xyz
            sym = []   # will contain list of atom symbols in the mol # same order
            Rlist = []
            sym_coord = []
            for tmp_l in range(3,num_l_xyz+1):   
                l1_xyz_2 = linecache.getline("geom_cation_IP.xyz",tmp_l)
                l1_lsp2 = l1_xyz_2.split()
                full_l = l1_xyz_2.strip()
                sym_coord.append(full_l)  ## contains atom and coords
                sym.append(l1_lsp2[0])
                Rlist.append(l1_lsp2[1])
                Rlist.append(l1_lsp2[2])
                Rlist.append(l1_lsp2[3])
                conv = np.array(Rlist)
                R_coord = np.resize(conv,[Nat,3])
            inp_x.write(sym_coord)  

            val["Ntotal"] = 0
            val["Ntotale"] = 0
            val["Ntotalecore"] = 0
            # write(10,'(a2,5x,3f15.8)') sym(iat), R(iat,1:3)
            ###read(90,*)sym(iat), R(iat,1:3)  -- whats happening?
            for tmp_k in range(Nat):
                na_nb_l = nanb(sym[tmp_k])
                na = na_nb_l[0]
                nb = na_nb_l[1]
                val["Ntotal"] = val["Ntotal"] + a + b
                val["Ntotale"] = val["Ntotale"] + atno(sym[tmp_j])
                val["Ntotalecore"] = val["Ntotalecore"] + NFC(sym[tmp_j])

    val["Ntotal"] = val["Ntotal"] - charge
    val["Ntotale"] = val["Ntotale"] - charge
    val["Na"] = (multip + val["Ntotal"] - 1 ) / 2.0
    val["Nb"] = val["Na"] + 1 - multip

    val["isatom"] == "false"
    if Nat == 1: val["isatom"] == "true"

    rung4mp2(val, start_time_main)
    os.system("cp freq0.txt freq_cation_IP.txt")

    U0mol2 = U0mol - U0
    UTmol2 = UTmol - UT
    HTmol2 = HTmol - HT
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *---------*\n")
        ther_chem.write(" * SUMMARY *\n")
        ther_chem.write(" *---------*\n\n")
        ther_chem.write(" * Ionization Energy (mol+ - mol) = " + str(-U0mol2) + " Hartree\n")
        ther_chem.write(" * Ionization Energy (mol+ - mol) = " + str(-U0mol2*au2ev) + " eV\n")
        ther_chem.write(" * Ionization Energy (mol+ - mol) = " + str(-U0mol2*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * Ionization Energy (mol+ - mol) = " + str(-U0mol2*au2kjm) + " kj/mol\n\n")
        
        ther_chem.write("* Breakdown *\n")
        ther_chem.write(" * CCSD(T) (mol+ - mol)           = " + str((IPbreakdown[2][1]-IPbreakdown[1][1])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * DE (MP2)(mol+ - mol)           = " + str((IPbreakdown[2][3]-IPbreakdown[1][3])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * DE (HF) (mol+ - mol)           = " + str((IPbreakdown[2][2]-IPbreakdown[1][2])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * HLC     (mol+ - mol)           = " + str((IPbreakdown[2][4]-IPbreakdown[1][4])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * SO      (mol+ - mol)           = " + str((IPbreakdown[2][5]-IPbreakdown[1][5])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * ZPVE    (mol+ - mol)           = " + str((IPbreakdown[2][6]-IPbreakdown[1][6])*au2kcm) + " kcal/mol\n")
        ther_chem.write("                                    -------------------------\n")
        #ther_chem.write("TOTAL                          = " + str(sum(IPbreakdown(2,:)-IPbreakdown(1,:))*au2kcm) +  " kcal/mol\n")
        ther_chem.write("                                    -------------------------\n\n")
        ther_chem.write("** Elapsed time = " + str(end_time_n - start_time_main) + "  s\n")

if val["calc_EA"] == "true":
    if val["vertical_EA"] == "true": val["verticalIP"] = val["verticalEA"]
    val["IPss"] == "true"

    with open("geom_anion_EA.xyz", "r") as g_cat:  # need to look at this file
        with open("input.xyz", "w") as inp_x:
            Nat= int(linecache.getline("geom_anion_EA.xyz",1).strip())
            tr = linecache.getline("geom_anion_EA.xyz",2).strip()
            lsp11 = tr.split()
            charge = lsp11[0]
            multip = lsp11[1]
            inp_x.write(str(Nat) + " \n")
            inp_x.write("0  " + str(multip) + " \n")
#            inp_x.write(sym_coord)  # check
            sym = []   # will contain list of atom symbols in the mol # same order
            Rlist = []
            sym_coord = []
            for tmp_l in range(3,num_l_xyz+1):   
                l1_xyz_2 = linecache.getline("geom_anion_EA.xyz",tmp_l)
                l1_lsp2 = l1_xyz_2.split()
                full_l = l1_xyz_2.strip()
                sym_coord.append(full_l)  ## contains atom and coords
                sym.append(l1_lsp2[0])
                Rlist.append(l1_lsp2[1])
                Rlist.append(l1_lsp2[2])
                Rlist.append(l1_lsp2[3])
                conv = np.array(Rlist)
                R_coord = np.resize(conv,[Nat,3])
            inp_x.write(sym_coord)  


            val["Ntotal"] = 0
            val["Ntotale"] = 0
            val["Ntotalecore"] = 0
            ###read(90,*)sym(iat), R(iat,1:3)  -- whats happening?
            for tmp_k in range(Nat):
                na_nb_l = nanb(sym[tmp_k])
                na = na_nb_l[0]
                nb = na_nb_l[1]
                val["Ntotal"] = val["Ntotal"] + a + b
                val["Ntotale"] = val["Ntotale"] + atno(sym[tmp_j])
                val["Ntotalecore"] = val["Ntotalecore"] + NFC(sym[tmp_j])

    val["Ntotal"] = val["Ntotal"] - charge
    val["Ntotale"] = val["Ntotale"] - charge
    val["Na"] = (multip + val["Ntotal"] - 1 ) / 2
    val["Nb"] = val["Na"] + 1 - multip

    val["isatom"] == "false"
    if Nat == 1: val["isatom"] == "true"

    rung4mp2(val, start_time_main)
    os.system("cp freq0.txt freq_anion_EA.txt")

    U0mol2 = U0mol - U0
    UTmol2 = UTmol - UT
    HTmol2 = HTmol - HT
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *---------*\n")
        ther_chem.write(" * SUMMARY *\n")
        ther_chem.write(" *---------*\n\n")
        ther_chem.write(" * Electron affinity (mol+ - mol -) = " + str(U0mol2) + " Hartree\n")
        ther_chem.write(" * Electron affinity (mol+ - mol -) = " + str(U0mol2*au2ev) + " eV\n")
        ther_chem.write(" * Electron affinity (mol+ - mol -) = " + str(U0mol2*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * Electron affinity (mol+ - mol -) = " + str(U0mol2*au2kjm) + " kj/mol\n\n")
        
        ther_chem.write("* Breakdown *\n")
        ther_chem.write(" * CCSD(T) (mol+ - mol)           = " + str((EAbreakdown[2][1]-EAbreakdown[1][1])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * DE (MP2)(mol+ - mol)           = " + str((EAbreakdown[2][3]-EAbreakdown[1][3])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * DE (HF) (mol+ - mol)           = " + str((EAbreakdown[2][2]-EAbreakdown[1][2])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * HLC     (mol+ - mol)           = " + str((EAbreakdown[2][4]-EAbreakdown[1][4])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * SO      (mol+ - mol)           = " + str((EAbreakdown[2][5]-EAbreakdown[1][5])*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * ZPVE    (mol+ - mol)           = " + str((EAbreakdown[2][6]-EAbreakdown[1][6])*au2kcm) + " kcal/mol\n")
        ther_chem.write("                                    -------------------------\n")
        #ther_chem.write("TOTAL                          = " + str(sum(IPbreakdown(1,:)-IPbreakdown(2,:))*au2kcm) +  " kcal/mol\n")
        ther_chem.write("                                    -------------------------\n\n")
        ther_chem.write("** Elapsed time = " + str(end_time_n - start_time_main) + "  s\n")

if val["calc_PA"] == "true":
    val["verticalIP"] = "false"
    val["IPss"] = "true"

    with open("geom_cation_PA.xyz", "r") as g_cat:  # need to look at this file
        Nat= int(linecache.getline("geom_cation_PA.xyz",1).strip())
        tr = linecache.getline("geom_cation_PA.xyz",2).strip()
        lsp11 = tr.split()
        charge = lsp11[0]
        multip = lsp11[1]
    
#    deallocate(R, sym)
#    allocate(R(1:Nat,1:3), sym(1:Nat))  # is this section needed ! IMPORTANT revisit
        sym = []   
        Rlist = []
        sym_coord = []
        for tmp_l in range(3,num_l_xyz+1):   
            l1_xyz_2 = linecache.getline("geom_cation_PA.xyz",tmp_l)
            l1_lsp2 = l1_xyz_2.split()
            full_l = l1_xyz_2.strip()
            sym_coord.append(full_l)  ## contains atom and coords
            sym.append(l1_lsp2[0])
            Rlist.append(l1_lsp2[1])
            Rlist.append(l1_lsp2[2])
            Rlist.append(l1_lsp2[3])
            conv = np.array(Rlist)
            R_coord = np.resize(conv,[Nat,3])
        #inp_x.write(sym_coord)  

    val["Ntotal"] = 0
    val["Ntotale"] = 0
    val["Ntotalecore"] = 0
    ###read(90,*)sym(iat), R(iat,1:3)  -- whats happening?
    for tmp_k in range(Nat):
        na_nb_l = nanb(sym[tmp_k])
        na = na_nb_l[0]
        nb = na_nb_l[1]
        val["Ntotal"] = val["Ntotal"] + a + b
        val["Ntotale"] = val["Ntotale"] + atno(sym[tmp_j])
        val["Ntotalecore"] = val["Ntotalecore"] + NFC(sym[tmp_j])

    val["Ntotal"] = val["Ntotal"] - charge
    val["Ntotale"] = val["Ntotale"] - charge
    val["Na"] = (multip + val["Ntotal"] - 1 ) / 2.0
    val["Nb"] = val["Na"] + 1 - multip

    val["isatom"] == "false"
    if Nat == 1: val["isatom"] == "true"

    uniq_atom_res = uniqatoms(sym)

    rung4mp2(val, start_time_main)
    os.system("cp freq0.txt freq_cation_PA.txt")

    U0mol2 = U0mol - U0
    UTmol2 = UTmol - UT
    HTmol2 = HTmol - HT

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *---------* \n")
        ther_chem.write(" * SUMMARY * \n")
        ther_chem.write(" *---------* \n\n")
        
        ther_chem.write(" * Proton affinity (mol - molH+) = " + str(HTmol2)        + " Hartree\n")
        ther_chem.write(" * Proton affinity (mol - molH+) = " + str(HTmol2*au2e)   + " eV\n")
        ther_chem.write(" * Proton affinity (mol - molH+) = " + str(HTmol2*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * Proton affinity (mol - molH+) = " + str(HTmol2*au2kjm) + " kj/mol\n\n")

        ther_chem.write(" * Breakdown *\n")
        ther_chem.write(" * CCSD(T)         (mol - molH+) = " + str((PAbreakdown(1,1)-PAbreakdown(2,1))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * DE (MP2)        (mol - molH+) = " + str((PAbreakdown(1,3)-PAbreakdown(2,3))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * DE (HF)         (mol - molH+) = " + str((PAbreakdown(1,2)-PAbreakdown(2,2))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * HLC             (mol - molH+) = " + str((PAbreakdown(1,4)-PAbreakdown(2,4))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * SO              (mol - molH+) = " + str((PAbreakdown(1,5)-PAbreakdown(2,5))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * ZPVE            (mol - molH+) = " + str((PAbreakdown(1,6)-PAbreakdown(2,6))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * Vib-thermal     (mol - molH+) = " + str((PAbreakdown(1,7)-PAbreakdown(2,7))*au2kcm)   + " kcal/mol\n")  # rewrite wont work
        ther_chem.write(" * Rot-thermal     (mol - molH+) = " + str((PAbreakdown(1,8)-PAbreakdown(2,8))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * Trans-thermal   (mol - molH+) = " + str((PAbreakdown(1,9)-PAbreakdown(2,9))*au2kcm)   + " kcal/mol\n")
        ther_chem.write(" * kBT             (mol - molH+) = " + str((PAbreakdown(1,10)-PAbreakdown(2,10))*au2kcm) + " kcal/mol\n")
        ther_chem.write("                                   -------------------------\n")
        #ther_chem.write("* TOTAL                         = " + str(sum(PAbreakdown(1,:)-PAbreakdown(2,:))*au2kcm) + " kcal/mol\n"  # wont work
        ther_chem.write("                                   -------------------------\n\n")
        ther_chem.write(" ** Elapsed time = " + str(end_time - start_time_main) + " s\n\n")

if val["calc_BE"] == "true":
    val["verticalIP"] = "false"
    val["IPss"] = "true"

    with open("geom_monomer_A.xyz", "r") as g_cat:  # need to look at this file
        Nat= int(linecache.getline("geom_monomer_A.xyz",1).strip())
        tr = linecache.getline("geom_monomer_A.xyz",2).strip()
        lsp11 = tr.split()
        charge = lsp11[0]
        multip = lsp11[1]
#    deallocate(R, sym)
#    allocate(R(1:Nat,1:3), sym(1:Nat))  # is this section needed ! IMPORTANT revisit
        sym = []   
        Rlist = []
        sym_coord = []
        for tmp_l in range(3,num_l_xyz+1):   
            l1_xyz_2 = linecache.getline("geom_monomer_A.xyz",tmp_l)
            l1_lsp2 = l1_xyz_2.split()
            full_l = l1_xyz_2.strip()
            sym_coord.append(full_l)  ## contains atom and coords
            sym.append(l1_lsp2[0])
            Rlist.append(l1_lsp2[1])
            Rlist.append(l1_lsp2[2])
            Rlist.append(l1_lsp2[3])
            conv = np.array(Rlist)
            R_coord = np.resize(conv,[Nat,3])
        #inp_x.write(sym_coord)  

    val["Ntotal"] = 0
    val["Ntotale"] = 0
    val["Ntotalecore"] = 0
    ###read(90,*)sym(iat), R(iat,1:3)  -- whats happening?
    for tmp_k in range(Nat):
        na_nb_l = nanb(sym[tmp_k])
        na = na_nb_l[0]
        nb = na_nb_l[1]
        val["Ntotal"] = val["Ntotal"] + a + b
        val["Ntotale"] = val["Ntotale"] + atno(sym[tmp_j])
        val["Ntotalecore"] = val["Ntotalecore"] + NFC(sym[tmp_j])

    val["Ntotal"] = val["Ntotal"] - charge
    val["Ntotale"] = val["Ntotale"] - charge
    val["Na"] = (multip + val["Ntotal"] - 1 ) / 2.0
    val["Nb"] = val["Na"] + 1 - multip

    val["isatom"] == "false"
    if Nat == 1: val["isatom"] == "true"

    uniq_atom_res = uniqatoms(sym)

    rung4mp2(val, start_time_main)
    os.system("cp freq0.txt freq_monomer_A.txt")

    U0mol2 = U0mol - (2*U0)
    UTmol2 = UTmol - (2*UT)
    HTmol2 = HTmol - (2*HT)

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *---------* \n")
        ther_chem.write(" * SUMMARY * \n")
        ther_chem.write(" *---------* \n\n")
        
        ther_chem.write(" * Binding energy (A2 - 2 * A) = " +  str(U0mol2)        + " Hartree\n")
        ther_chem.write(" * Binding energy (A2 - 2 * A) = " +  str(U0mol2*au2ev)  + " eV\n")
        ther_chem.write(" * Binding energy (A2 - 2 * A) = " +  str(U0mol2*au2kcm) + " kcal/mol\n")
        ther_chem.write(" * Binding energy (A2 - 2 * A) = " +  str(U0mol2*au2kjm) + " kj/mol\n\n")
        
        ther_chem.write(" ** Elapsed time = " + str(end_time_n - start_time_main) + " s\n\n")  # which end time!! ??
        
# !=== Atomic calculations

if val["calc_HF"] == "true":
    val["calc_AE"] == "true"

if val["calc_AE"] == "true":
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * G4MP2 for reference atoms*\n\n")

    U0molAE = U0mol
    UTmolAE = UTmol
    HTmolAE = HTmol
    #uniq_atom_res = uniqatoms(sym)
    for tmp_o in range(uniq_atom_res["N_ua"]):
        Nat = 1
        # deallocate(sym, R)
        sym =[]
        ua = uniq_atom_res["uniq_sym"]
        if Nat == 1: val["isatom"] = "true"
        sym.append(ua[tmp_o])
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write(" * ATOM:" + str(sym) + "\n")
        
        na_nb_l = nanb(sym[0]) # check
        na = na_nb_l[0]
        na = na_nb_l[1]
        R = 0 ## ask
        multip = na - nb + 1
        charge = 0
        with open("input.xyz", "w") as new_in:
            new_in.write(str(Nat) + "\n")
            new_in.write("0  " + str(multip) + "\n")
            new_in.write(str(sym[0]) + "0.0 0.0 0.0 \n")  # ask

        val["Ntotale"] = atno(ua[tmp_o])
        val["Ntotalecore"] = NFC(ua[tmp_o])
        rung4mp2(val,start_time_main)
        uan_l = uniq_atom_res["uan"]
        U0molAE = U0molAE - uan_l[tmp_o] * U0
        UTmolAE = UTmolAE - uan_l[tmp_o] * UT
        HTmolAE = HTmolAE - uan_l[tmp_o] * HT
        HLCmol =  HLCmol -  uan_l[tmp_o] * HLC0

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *---------* \n")
        ther_chem.write(" * SUMMARY * \n")
        ther_chem.write(" *---------* \n\n")

        ther_chem.write(" * Atomization energy (Atoms - mol) = " + str(-U0molAE)  + " Hartree\n")
        ther_chem.write(" * Atomization energy (Atoms - mol) = " + str(-U0molAE*au2ev) + " eV\n")
        ther_chem.write(" * Atomization energy (Atoms - mol) = " + str(-U0molAE*au2kcm) +  " kcal/mol\n")
        ther_chem.write(" * Atomization energy (Atoms - mol) = " + str(-U0molAE*au2kjm) + " kj/mol\n")
        ther_chem.write(" * Atomization energy (Atoms - mol), NO HLC = " + str(-(U0molAE-HLCmol)*au2kcm) +  " kcal/mol\n\n")
        
        ther_chem.write(" ** Elapsed time = " + str(end_time_n - start_time_main) + " s\n\n")  # which end time!! ??

if val["calc_HF"] == "true":
    # allocate(HOF_0K(1:N_ua), dH_298K_0K(1:N_ua))
    HOF = HTmol - U0mol + U0molAE
    #uniq_atom_res = uniqatoms(sym)
    HOF_0K = []
    dH_298K_0K = []
    ua = uniq_atom_res["uniq_sym"]
    for tmp_u in range(uniq_atom_res["N_ua"]):
        HOF_0K.append(HOF_atoms(ua[tmp_u]))
        dH_298K_0K.append(dH_atoms(ua[tmp_u]))
    
#    if val["read_HOF_params"] == "true":
#      readloop: do
#        read(999,*,iostat=stat) tmp_ch, tmp_dp1, tmp_dp2
#        if ( stat .ne. 0) exit readloop
#        readloop2: do iat = 1, N_ua
#          if ( trim(tmp_ch) .eq. trim(ua(iat)) ) then
#            HOF_0K(iat) = tmp_dp1 * kcm2au
#            dH_298K_0K(iat) = tmp_dp2 * kcm2au
#            cycle readloop2
#          endif
#        enddo readloop2
#      enddo readloop
#    endif
#   close(999)
#### IMPORTANT revisit
    uan_l = uniq_atom_res["uan"]
    HOF = HTmol - U0mol + U0molAE
    for tmp_u in range(uniq_atom_res["N_ua"]):
        HOF = HOF + uan_l[tmp_u] * ( HOF_0K[tmp_u] - dH_298K_0K[tmp_u] )

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *---------* \n")
        ther_chem.write(" * SUMMARY * \n")
        ther_chem.write(" *---------* \n\n")

        ther_chem.write("* Heat of formation = " + str(HOF) + " Hartree\n")
        ther_chem.write("* Heat of formation = " + str(HOF*au2ev) + " eV\n")
        ther_chem.write("* Heat of formation = " + str(HOF*au2kcm) + " kcal/mol\n")
        ther_chem.write("* Heat of formation = " + str(HOF*au2kjm) + " kj/mol\n")
        ther_chem.write("* Heat of formation, NO HLC = " +  str((HOF-HLCmol)*au2kcm) + " kcal/mol\n\n")
        
        ther_chem.write("** Elapsed time = " +str(end_time_n - start_time_main) + " s\n\n")  # which end time!! ??

#deallocate(R, sym, ua, uan)
#if val["isatom"] != "true":
    #if ( allocated(freq) ) deallocate(freq, theta)


day = datetime.datetime.today().strftime("%A")
date = datetime.date.today()
cur_date = date.strftime("%B %d, %Y")
now = datetime.datetime.now()
cur_time = now.strftime(" %H:%M:%S")

with open("Thermochemistry.out", "a") as ther_chem:
    ther_chem.write("==============================================\n")
    ther_chem.write("Final time:  " + day + "  "+ cur_date + "  " + cur_time + "  " + time.tzname[0] +"  \n\n")
    ther_chem.write("==============================================\n")

    
    
