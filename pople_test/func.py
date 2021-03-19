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

def atno(at):
    at_num_dict ={"H": 1,"He": 2,"Li": 3,"Be": 4,"B": 5,"C": 6,"N": 7,"O": 8,"F": 9,"Ne": 10,"Na": 11,"Mg": 12,"Al": 13,"Si": 14,"P": 15,"S": 16,"Cl": 17,"Ar": 18,"K": 19,"Ca": 20, \
        "Sc": 21,"Ti": 22,"V": 23,"Cr": 24,"Mn": 25,"Fe": 26,"Co": 27,"Ni": 28,"Cu": 29,"Zn": 30,"Ga": 31,"Ge": 32,"As": 33,"Se": 34,"Br": 35,"Kr": 36, "I":53}

    if at in at_num_dict:
        atno = at_num_dict[at]
        return(atno)
    else:
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in atno in module_geom.f90: " + str(at)+ " \n")

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
def nanb(at):
    na_dict_1 = {"H" : 1,"He" : 1,"Li" : 1,"Be" : 1,"B" : 2,"C" : 3,"N" : 4,"O" : 4,"F" : 4,"Ne" : 4,"Na" : 1,"Mg" : 1,"Al" : 2,"Si" : 3,"P" : 4,"S" : 4,"Cl" : 4,"Ar" : 4,"K" : 1,"Ca" : 1,\
        "Ga" : 2,"Ge" : 3,"As" : 4,"Se" : 4,"Br" : 4,"Kr" : 4,"Fe" : 10,"I" : 4}
    nb_dict_1 = {"H": 0,"He": 1, "Li" : 0,"Be": 1,"B": 1,"C": 1,"N": 1,"O": 2,"F": 3,"Ne": 4,"Na": 0,"Mg": 1,"Al": 1,"Si": 1,"P": 1,"S": 2,"Cl": 3,"Ar": 4,"K": 0,"Ca": 1,"Ga": 1,"Ge": 1,"As": 1,\
        "Se": 2,"Br": 3,"Kr": 4,"Fe": 6,"I": 3,}
    if at in na_dict_1:
        na = na_dict_1[at]
        nb = nb_dict_1[at]
        return(na,nb)
    else:  # equivalent to case default, if at not in na_dict_1, this is printed, revist
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in nanb in module_geom.f90:  " + str(at) + " \n")

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

import os, sys, string
import linecache, math
import numpy as np


def printbas(fname, at, install_dir_val):  # fname = value of the  basis set name from the inp file,  install_dir_val = value from initial input file by user
    #print("===Entering printbas with: ", fname, at, " ===")
    basisSet_fpath = fname
    start_phrase = "NewGTO   "+ at
    print(start_phrase)
    num_lines_bas = sum(1 for line_tmp1 in open(basisSet_fpath,"r"))

    for temp_num, temp_l in enumerate(open(basisSet_fpath,"r")):
        if start_phrase in temp_l.strip():
            if start_phrase == temp_l.strip():
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
def principal_coord():
    with open("input.xyz","r") as xyz_f:
        num_l_xyz = sum(1 for l in xyz_f)
        Nat = int(linecache.getline("input.xyz",1).strip())
        l1_chr_mul = linecache.getline("input.xyz",2)
        sym = []   # will contain list of atom symbols in the mol # same order
        Rlist = []
        for tmp_l in range(3,num_l_xyz+1):   ### what does NH do ???
            l1_xyz_1 = linecache.getline("input.xyz",tmp_l)
            l1_lsp1 = l1_xyz_1.split()
            sym.append(l1_lsp1[0])
            Rlist.append(float(l1_lsp1[1]))
            Rlist.append(float(l1_lsp1[2]))
            Rlist.append(float(l1_lsp1[3]))
        conv = np.array(Rlist)
        R_coord = np.resize(conv,[Nat,3])

        rCM_pre = []
        mass_list = []
        for tmp_j in range(Nat):
            #print(sym[tmp_j])
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
        for tmp_k in range(Nat):
            sub_cm = R_coord[tmp_k] - rCM
            new_coord1.append(sub_cm)
            r_skew_symm[0][0] = 0.0
            r_skew_symm[0][1] = -sub_cm[2]
            r_skew_symm[0][2] = sub_cm[1]

            r_skew_symm[1][0] = sub_cm[2]
            r_skew_symm[1][1] = 0.0
            r_skew_symm[1][2] =-sub_cm[0]

            r_skew_symm[2][0] =-sub_cm[1]
            r_skew_symm[2][1] = sub_cm[0]
            r_skew_symm[2][2] = 0.0

            momin = momin + (sym2mass(sym[tmp_k]) * np.matmul(np.transpose(r_skew_symm), r_skew_symm) )

        eig_val, eig_vec = np.linalg.eig(momin)

        idx = eig_val.argsort()[::+1]
        eig_val = eig_val[idx]
        eig_vec = eig_vec[:,idx]

        Ievals = eig_val


        print("MOM== ", Ievals)
        return(Ievals)
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
        #print(uniqat_d)
    #return(N_ua, uniq_sym, uan)  ### revisit, how do you want variables returned?
    return(uniqat_d)
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
import os, sys, string
import linecache, math
import numpy as np
import datetime , time


from pople_test import NFC
from pople_test import uniqatoms
from pople_test import printbas


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

        with open("inp.xyz","r") as xyz_f:
            num_l_xyz = sum(1 for l in xyz_f)
            linecache.clearcache()
            Nat = int(linecache.getline("inp.xyz",1).strip())
            l1_chr_mul = linecache.getline("inp.xyz",2)
            #xyz_l1 = "* xyz  " + l1_chr_mul.strip() + " \n"
            if values["restart_check"] == "true":
                 test_l = "* xyz  " + str(values["restart_charge"]) + "  " + str(multip) + " \n"
            else:
                 test_l = "* xyz  " + l1_chr_mul.strip() + " \n"
            #test_l = "* xyz  " + l1_chr_mul.strip() + " \n"
            com_f.write(test_l)
            #com_f.write(xyz_l1)
            sym = []   # will contain list of atom symbols in the mol # same order
            
            for tmp_l in range(3,num_l_xyz+1):   ### what does NH do ???
                l1_xyz_1 = linecache.getline("inp.xyz",tmp_l)
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
            uniq_atom_res = uniqatoms.uniqatoms(sym)
            if values["ALLELE"] == "true":  ### CHECK!!!!
                for iat in range(int(uniq_atom_res["N_ua"])):
                    pre1 = uniq_atom_res["uniq_sym"]
                    at_pr1 = pre1[iat]
                    com_f.write("  NewNCore " + at_pr1 + "  " + " 0  end\n")
            else:
                for iat in range(int(uniq_atom_res["N_ua"])):
                    pre1 = uniq_atom_res["uniq_sym"]
                    at_pr1 = pre1[iat]
                    NFC_res = NFC.NFC(at_pr1)
                    com_f.write("  NewNCore " + at_pr1 + "  " + str(NFC_res) +"  end\n")
        
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
            uniq_atom_res = uniqatoms.uniqatoms(sym)
            fname = values["install_dir"] + basis
            #print("fname is ", fname)
            ##write(fname,'(a,a)')trim(adjustl(install_dir)), trim(adjustl(basis))  ## IMPORTANT
            #print(uniq_atom_res)
            if Nat == 1:  # check datatype
                printbas.printbas(fname, sym[0], values["install_dir"])    ## what is fname, it should be the basis set keyword. but which one??
            else:
                for iat1 in range(int(uniq_atom_res["N_ua"])):
                    pre2 = uniq_atom_res["uniq_sym"]
                    at_pr2 = pre2[iat1]
                    #print(iat1, at_pr2)
                    printbas.printbas(fname, at_pr2, values["install_dir"])
            com_f.write("end\n")
#   os.system("cat basis.com >> input.com")

    print(values["orca_dir"] + "orca input.com > input.out")
    os.system(values["orca_dir"] + "orca input.com > input.out")
    os.system("cat input.com >> ORCA_G4MP2.com")
    os.system("cat input.out >> ORCA_G4MP2.out")
####### runorca - E
import os, sys, string
import linecache, math
import numpy as np
import datetime , time


from pople_test import Mol_SO
from pople_test import At_SO
from pople_test import principal_coord
from pople_test import runorca


####### rung4mp2 - S
def rung4mp2(values, start_time_main):  
    if values["isatom"] == "true": Nat=1
    os.system("cat inp.xyz")


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
                    #full_l = l1_xyz_1.strip()
                    full_l1 = l1_xyz_1.strip().split()
                    full_l = full_l1[0] + "     "+str(round(float(full_l1[1]),8)) + "     "+str(round(float(full_l1[2]),8)) + "     "+str(round(float(full_l1[3]),8))
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

                print("freq check")
                print(nfreq)

## How does read_geom_freq.dat look like????
#       allocate(freq(1:nfreq), theta(1:nfreq))
#       do iat = 1, nfreq
#         !print *, iat
#         read(10,*)freq(iat)
#       enddo

            if Nat > 2:
                Ievals = principal_coord.principal_coord() # check
                linnonlin = "NL"
                print("FROZEN_GEOM", FROZEN_GEOM)
                print("Eigenvalues of Moment of inertia:", Ievals)
                if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                    linnonlin = "L"
            if Nat == 2:    linnonlin = "L"

            if linnonlin == "L":    nfreq = 3 * Nat - 5
            if linnonlin == "NL":   nfreq = 3 * Nat - 6


### revisit line 710 
        else:

    #call system_clock(count_start1)
            start_time1 = time.time()
            print("==check isatom==",values["isatom"])
            runorca.runorca(values["method_opt_freq"], values["basis_opt_freq"], "true", values["custombasis_opt_freq"], "false", values)
            end_time1 = time.time()
    #call system_clock(count_end1, count_rate1)
            if values["IPss"] != "true" or values["verticalIP"] != "true":
                with open("inp.xyz", "r") as new_i_xyz:
                    num_l_in = sum(1 for l in new_i_xyz)
                    Nat = int(linecache.getline("inp.xyz",1).strip())
                    title = linecache.getline("inp.xyz",2).strip()
                    tl1 = title.split()
                    charge = int(tl1[0])
                    multip = int(tl1[1])
                with open("input.xyz", "r") as new_i_xyz:
                    num_l_in = sum(1 for l in new_i_xyz)
                    Nat = int(linecache.getline("input.xyz",1).strip())
                    title = linecache.getline("input.xyz",2).strip()
                    tl1 = title.split()
                   #charge = int(tl1[0])
                   #multip = int(tl1[1])
                    sym = []
                    Rlist = []
                    sym_coord = []
                    for tmp_l in range(3,num_l_in+1):    ## check range
                        l1_xyz_1 = linecache.getline("input.xyz",tmp_l)
                        l1_lsp1 = l1_xyz_1.split()
                        #full_l = l1_xyz_1.strip()
                        full_l1 = l1_xyz_1.strip().split()
                        full_l = full_l1[0] + "     "+str(round(float(full_l1[1]),8)) + "     "+str(round(float(full_l1[2]),8)) + "     "+str(round(float(full_l1[3]),8))
                        sym_coord.append(full_l)  ## contains atom and coords
                        sym.append(l1_lsp1[0])
                        Rlist.append(l1_lsp1[1])
                        Rlist.append(l1_lsp1[2])
                        Rlist.append(l1_lsp1[3])
                        conv = np.array(Rlist)
                        R_coord = np.resize(conv,[Nat,3])
                with open("temp.xyz", "w") as t_xyz:
                   e1= linecache.getline("inp.xyz",1).strip()
                   t_xyz.write(e1 + "\n")
                   e2= linecache.getline("inp.xyz",2).strip()
                   t_xyz.write(e2 + "\n")
                   for e3 in range(len(sym_coord)):
                       t_xyz.write(sym_coord[e3] + "\n")
            os.system("mv temp.xyz  inp.xyz")

#head -2 inp.xyz > inp.xyz
#tail -Nat input.xyz >> inp.xyz

############################################# REWRITE???
        if Nat > 2:
            Ievals = principal_coord.principal_coord() # check
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
            if nfreq > (3*Nat)-6:
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
        print("freq check 2", freq)
   

        iopt = 999
        if values["isatom"] != "true":
            os.system("grep 'OPTIMIZATION RUN DONE' input.out | wc | awk '{print $1}' > scr_lc1")
        iopt = int(linecache.getline("scr_lc1",1).strip())
 
        os.system("rm -f input* scr freq.txt")
        #os.system("rm -f scr freq.txt")
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
            #print(nfreq)
            for tmp in range(nfreq): # check range
                #print(type(values["scalfac"]))# * values["scalfac"])
                ther_chem.write(str( freq[tmp] ))  
                ther_chem.write("\n")

            ther_chem.write("\n")
            ther_chem.write(" Scaled harmonic wavenumbers (cm^-1)\n")
            tmp = 0
            #print(nfreq)
            for tmp in range(nfreq): # check range
                #print(type(values["scalfac"]))# * values["scalfac"])
                ther_chem.write(str( freq[tmp] * float(values["scalfac"]) ) )  
                ther_chem.write("\n")

            ther_chem.write("\n")
            ther_chem.write(" Scaling factor used:  " + str(values["scalfac"]) + "\n\n")
            #call system_clock(count_end, count_rate)
            end_time2 = time.time()
            ther_chem.write(" * Geometry optimization/frequencies done in              " + str(end_time1 - start_time1) + " s\n")
            ther_chem.write(" ** Elapsed time =              " + str(round(end_time2 -start_time_main,2)) + " s\n")
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
# major addition - not in f90 , check if other numbers change
#  if HOF, if not atoms
    with open("inp.xyz", "r") as new_i_xyz:
        linecache.clearcache()
        num_l_in = sum(1 for l in new_i_xyz)
        Nat = int(linecache.getline("inp.xyz",1).strip())
        title = linecache.getline("inp.xyz",2).strip()
        tl1 = title.split()
        charge = int(tl1[0])
        multip = int(tl1[1])
  # with open("inp.xyz", "r") as new_i_xyz:
  #     num_l_in = sum(1 for l in new_i_xyz)
  #     Nat = int(linecache.getline("inp.xyz",1).strip())
  #     title = linecache.getline("inp.xyz",2).strip()
  #     tl1 = title.split()
       #charge = int(tl1[0])
       #multip = int(tl1[1])
        sym = []
        Rlist = []
        sym_coord = []
        for tmp_l in range(3,num_l_in+1):    ## check range
            l1_xyz_1 = linecache.getline("inp.xyz",tmp_l)
            l1_lsp1 = l1_xyz_1.split()
            #full_l = l1_xyz_1.strip()
            full_l1 = l1_xyz_1.strip().split()
            print(full_l1)
            full_l = full_l1[0] + "     "+str(round(float(full_l1[1]),8)) + "     "+str(round(float(full_l1[2]),8)) + "     "+str(round(float(full_l1[3]),8))
            sym_coord.append(full_l)  ## contains atom and coords
            sym.append(l1_lsp1[0])
            Rlist.append(l1_lsp1[1])
            Rlist.append(l1_lsp1[2])
            Rlist.append(l1_lsp1[3])
            conv = np.array(Rlist)
            R_coord = np.resize(conv,[Nat,3])
#    os.system("rm input*")


    start_time1 = time.time()
    if values["Ntotale"] > 0:
        if Nat == 1:
            if values["basis_ccsdt"]  == "GTBAS1": values["basis_ccsdt"] = 'GTBAS1atm'
        if values["restart_cc"] == "true":
            #charge = charge + 2
            tmp_method = "HF "
            values["restart_check"] = "true"
            values["restart_charge"] = charge + 2
            runorca.runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
            #charge = charge - 2
            values["restart_charge"] = values["restart_charge"] - 2
    
        if values["DLPNO_CCSDT"] == "true":
            values["switch_DLPNO_CCSDT"] = "true"
        
        runorca.runorca(values["method_ccsdt"], values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
        end_time1= time.time()

        
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc2")
#            print(values["Ntotale"], values["Ntotalecore"], "==HF==")
        else:
            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc2")
        #    print(values["Ntotale"], values["Ntotalecore"], "==CC==")
        
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
                    runorca.runorca(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                else:  # ASK
                    runorca.runorca(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
            elif (values["Ntotale"] == 2) or (values["Ntotale"] - values["Ntotalecore"] == 1 ):
                if Nat == 1:
                    if values["basis_ccsdt"] == "GTBAS1": 
                        values["basis_ccsdt"] = "GTBAS1atm"  # check
                        runorca.runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                    else:
                        runorca.runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
            else:
                if values["restart_cc"] == "true":
                    #charge = charge + 2
                    tmp_method = "HF "
                    values["restart_check"] = "true"
                    values["restart_charge"] = charge + 2
                    runorca.runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
                    #charge = charge - 2
                    values["restart_charge"] = values["restart_charge"] - 2

                    #charge = charge + 2
                    #tmp_method = "HF "
                    #runorca(tmp_method,values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                    #charge = charge - 2
                runorca.runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
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
    end_time1 = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write("\n")
        ther_chem.write(" * CCSD(T) done in              " + str(round(end_time1 - start_time1,2)) + " s\n")
        ther_chem.write(" ** Elapsed time =              " + str(round(end_time_n - start_time_main,2)) + " s \n\n")
        
    if values["switch_read_RIMP2_small"] != "true":
        os.system("rm -f input* scr ")
    
    if values["DLPNO_CCSDT"] == "true":
        values["switch_DLPNO_CCSDT"] = "false"

    print("block dlpno ccsd check", values["switch_read_RIMP2_small"], values["DLPNO_CCSDT"], "outside")
## == TESTING ( comm from f90)
    if values["DLPNO_CCSDT"] == "true":
        print("block dlpno ccsd check", values["switch_read_RIMP2_small"], values["DLPNO_CCSDT"], "inside")
        if values["switch_read_RIMP2_small"] != "true":
            print( multip , values["Ntotale"], values["Ntotalecore"], "ntotale ntotalecore")
            if (multip != 1) or (values["switch_RIMP2_small"] == "true") :
                start_time1 = time.time()
                print(values["Ntotale"], values["Ntotalecore"], "ntotale ntotalecore")
                if (values["Ntotale"] == "1") or (values["Ntotale"] - values["Ntotalecore"] <= 1) :
                    tmp_method = "HF "
                    print("block dlpno ccsd check  1")
                    runorca.runorca(values["method_mp2_s"], values["basis_mp2_s"], "false", values["custombasis_mp2_s"], "true", values)
                else:
                    print("block dlpno ccsd check  2")
                    runorca.runorca(values["method_mp2_s"], values["basis_mp2_s"], "false", values["custombasis_mp2_s"], "true", values)
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
    if (values["Ntotale"] > 0) and (values["ccsdt_cbs"] != "true") : ## where is ccsdt_cbs declared?
        if values["restart_mp2"] == "true":
            #charge = charge + 2
            tmp_method = "HF "
            values["restart_check"] = "true"
            values["restart_charge"] = charge + 2
            runorca.runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
            #charge = charge - 2
            values["restart_charge"] = values["restart_charge"] - 2

            #charge = charge + 2
            #tmp_method = "HF "
            #runorca(tmp_method, values["basis_mp2"], "false", values["custombasis_mp2"], "true", values)
            #charge = charge - 2
        runorca.runorca(values["method_mp2"], values["basis_mp2"], "false", values["custombasis_mp2"], "true", values)
        end_time1= time.time()
    
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc6")
        else:
            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc6")
        tmp = linecache.getline("scr_lc6",1).strip()
        E_mp2L = float(tmp)

        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            E_hfL = E_mp2L
            print("=== E_hfL === 00", E_hfL, E_mp2L, values["Ntotale"], values["Ntotalecore"])
        else:
            if values["flag_RIMP2"] == "true":
                os.system("grep 'RI-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            elif values["flag_DLPNOMP2"] == "true":
                os.system("grep 'DLPNO-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            else:
                os.system("grep 'MP2 CORRELATION ENERGY ' input.out | awk '{print $5}' > scr_lc7")
            tmp = linecache.getline("scr_lc7",1).strip()
            E_hfL = float(tmp)

            print("=== E_hfL ===", E_hfL, E_mp2L, values["Ntotale"], values["Ntotalecore"])
            E_hfL = E_mp2L - E_hfL

        os.system("rm -f input* scr ")
# VERIFY - indentation issues, which else this goes to??? revisit
    else:
        E_hfL= 0
        E_mp2L = 0

    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * MP2/L done in                " + str(round(end_time1 - start_time1,2)) + " s \n")
        ther_chem.write(" ** Elapsed time =              " + str(round(end_time_n - start_time_main,2)) + " s \n\n")

#== SINGLE POINT HF/VTZ
    start_time1 = time.time()
    if (values["Ntotale"] > 0) and (values["ccsdt_cbs"] != "true"):
        if values["restart_hf3"] == "true":
            #charge = charge + 2
            tmp_method = "HF "
            values["restart_check"] = "true"
            values["restart_charge"] = charge + 2
            runorca.runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
            #charge = charge - 2
            values["restart_charge"] = values["restart_charge"] - 2

            #charge = charge + 2
            #tmp_method = "HF "
            #runorca(tmp_method, values["basis_hf3"], "false", values["custombasis_hf3"], "false", values)
            #charge = charge - 2
        runorca.runorca(values["method_hf3"], values["basis_hf3"], "false", values["custombasis_hf3"], "false", values)
        end_time1 = time.time()

        os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc8")
        tmp = linecache.getline("scr_lc8",1).strip()
        E_hfT = float(tmp)
        
        os.system("rm -f input* scr ")
    else:
        E_hfT = 0
    
    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * HF/VTZ done in              " + str(round(end_time1 - start_time1,2)) + " s\n")
        ther_chem.write(" ** Elapsed time =              " + str(round(end_time_n - start_time_main,2)) + " s\n\n")

# == SINGLE POINT HF/VQZ
    start_time1 = time.time()
    if (values["Ntotale"] > 0 ) and (values["ccsdt_cbs"] != "true"):
        if values["restart_hf4"] == "true":
            #charge = charge + 2
            tmp_method = "HF "
            values["restart_check"] = "true"
            values["restart_charge"] = charge + 2
            runorca.runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
            #charge = charge - 2
            values["restart_charge"] = values["restart_charge"] - 2

            #charge = charge + 2
            #tmp_method = "HF "
            #runorca(tmp_method, values["basis_hf4"], "false", values["custombasis_hf4"], "false", values)
            #charge = charge - 2
        runorca.runorca(values["method_hf4"], values["basis_hf4"], "false", values["custombasis_hf4"], "false", values)
        end_time1 = time.time()

        os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc9")
        tmp = linecache.getline("scr_lc9",1).strip()
        E_hfQ = float(tmp)
        
        os.system("rm -f input* scr ")
    else:
        E_hfQ = 0
    
    end_time_n = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" * HF/VQZ done in              " + str(round(end_time1 - start_time1,2)) + " s\n")
        ther_chem.write(" ** Elapsed time =              " + str(round(end_time_n  - start_time_main,2)) + " s\n\n")

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
    print("Vib Rot Tra Tot = ",dE_vib,dE_rot,dE_tra,Ethermal)

#  !=== SO
    HLC_SO = 0
    if values["isatom"] != "true":
        HLC_SO = Mol_SO.Mol_SO(Nat, multip, charge, sym, values["SO_3rdrow_mols"])
        #print(HLC_SO)
    else:
        HLC_SO = At_SO.At_SO(sym[0], charge)
    
    print("==check near HLC ",values["isatom"],Nat)
    if values["HLCeqZERO"] == "true":
        HLC = 0
    else:
        if values["isatom"] != "true":
            if values["Na"] == values["Nb"]:
                print("HLC atoms-1:",values["Na"],values["Nb"])
                HLC = - values["AA"] * values["Nb"]
            else:
                print("HLC atoms-2:",values["Na"],values["Nb"])
                HLC = - values["ApAp"] * values["Nb"] - values["BB"] * (values["Na"]- values["Nb"])   # check
        else:
             print("HLC atoms-3:",values["Na"],values["Nb"])
             HLC = - (values["CC"] * values["Nb"] ) - (values["DD"] * (values["Na"] - values["Nb"]) )  # + At_SO(sym(1), charge)
        
        if values["isatom"] == "true":
            if sym[0] == "Be" and charge == 0 : HLC = - values["EE"]
            if sym[0] == "Mg" and charge == 0 : HLC = - values["EE"]
            if sym[0] == "Ca" and charge == 0 : HLC = - values["EE"]
            if sym[0] == "Li" and charge == -1 : HLC = - values["EE"]
            if sym[0] == "Na" and charge == -1 : HLC = - values["EE"]
            if sym[0] == "K" and charge == -1 : HLC = - values["EE"]
        elif Nat == 2:
            if charge == 0:
                if sym[0] == "Li" and sym[1] == "Li": HLC = - values["EE"]
                if sym[0] == "Na" and sym[1] == "Na": HLC = - values["EE"]
                if sym[0] == "K" and sym[1] == "K": HLC   = - values["EE"]
                if sym[0] == "Li" and sym[1] == "Na": HLC = - values["EE"]
                if sym[0] == "Na" and sym[1] == "Li": HLC = - values["EE"]
                # if ( (trim(sym(1)) .eq. 'Be') .and. (trim(sym(2)) .eq. 'H' ) ) HLC = -EE
                # if ( (trim(sym(1)) .eq. 'H' ) .and. (trim(sym(2)) .eq. 'Be') ) HLC = -EE
    #print(values["AA"],values["ApAp"],values["BB"],values["CC"],values["DD"],values["EE"],values["Na"],values["Nb"])
    HLC  = HLC / 1000.0
    HLC_SO = HLC_SO / 1000.0
    values["HLC0"] = HLC

# CBS extrapolation hard-coded for HF (AVQZ, AV5Z)
    if values["HF_CBS_default"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-1.63)) / (1 - math.exp(-1.63))
    if values["HF_CBS_orca_23_def2"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-10.39*(math.sqrt(3.0)- math.sqrt(2.0)))) / (1 - math.exp(-10.39*(math.sqrt(3.0)-math.sqrt(2.0))))
    if values["HF_CBS_orca_34_def2"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-7.88*( math.sqrt(4.0)- math.sqrt(3.0)))) /  (1 - math.exp(-7.88*(math.sqrt(4.0)-math.sqrt(3.0))))
    if values["HF_CBS_orca_23_cc"] == "true":   E_HF_limit = (E_hfQ - E_hfT * math.exp(-4.42*( math.sqrt(3.0)- math.sqrt(2.0)))) / (1 - math.exp(-4.42* (math.sqrt(3.0)-math.sqrt(2.0))))
    if values["HF_CBS_orca_34_cc"] == "true":   E_HF_limit = (E_hfQ - E_hfT * math.exp(-5.46*( math.sqrt(4.0)- math.sqrt(3.0)))) / (1 - math.exp(-5.46* (math.sqrt(4.0)-math.sqrt(3.0))))
    
    E_HF_CBS   = E_HF_limit - E_hfL
    #print(E_mp2L,E_mp2)
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

    values["U0"] = eps0 + ZPVE
    values["UT"] = eps0 + Ethermal
    values["HT"] = eps0 + Ethermal + kB * T * j2au

    if values["IPss"] == "true":   ## post check reduce all indices by 1
        print("-----here---")
        values["IPbreakdown2"] = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC,HLC_SO,ZPVE  ]
        values["EAbreakdown2"] = [ E_ccsdt, E_HF_CBS, E_MP2_CBS,HLC,HLC_SO,ZPVE  ]
        values["PAbreakdown2"] = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE, (dE_vib - ZPVE) , dE_rot, dE_tra , (kB * T * j2au) ]
    else:
        values["IPbreakdown1"] = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE ]
        values["EAbreakdown1"] = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE ]
        values["PAbreakdown1"] = [ E_ccsdt, E_HF_CBS, E_MP2_CBS, HLC, HLC_SO, ZPVE, (dE_vib - ZPVE), dE_rot, dE_tra, (kB * T * j2au)  ]

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write("           Temperature=     " + str(round(T,6)) +  "                   Pressure=       1.000000 \n")
        ther_chem.write("                E(ZPE)=       " +str(round(ZPVE,6)) + "             E(Thermal)=       " + str(round(Ethermal,6)) + "\n")  # revisit
        if values["G4MP2TM"] == "true":
            ther_chem.write("            E(CCSD(T),rel)=       "  + str(round(E_ccsdt_rel,6)) + "\n")
        ther_chem.write("                   HLC=      " + str(round(HLC,6)) + "                     SO=       " + str(round(HLC_SO,6)) + "  \n")
        ther_chem.write("            E(CCSD(T))=     " + str(round(E_ccsdt,6))  + "\n")
        ther_chem.write("               DE(MP2)=      "  + str(round(E_MP2_CBS,6)) + "                 DE(HF)=      " + str(round(E_HF_CBS,6))  + "\n")
        ther_chem.write("            G4MP2(0 K)=     "  + str(round(values["U0"],6))  + "           G4MP2 Energy=     " + str(round(values["UT"],6)) + "\n")
        ther_chem.write("        G4MP2 Enthalpy=     " + str(round(values["HT"],6)) + "      G4MP2 Free Energy=       0.000000 \n\n")  # verify
        
        end_time_n = time.time()
        ther_chem.write(" * G4(MP2) done\n")
        ther_chem.write(" ** Elapsed time =              "  + str(round(end_time_n - start_time_main,2))  + "  s \n\n")
        
    os.system("rm -f scr scrd")

####### rung4mp2 - E
