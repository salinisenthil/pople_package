import os, sys, string
import linecache, math
import numpy as np
import datetime , time


from pople import NFC
from pople import uniqatoms
from pople import printbas


####### runorca - S
def runorca(method, basis,optfreq,custombasis, correlated, values):
    """
    Runs orca

            Parameters:
                    method (char) : Name of functional to be used
                    basis (char)  :  Basis set name
                    optfreq (char) : true/false value of the optfreq keyword in control.inp
                    custombasis (char) : true/false value of the custombasis keyword in control.inp
                    correlated (char) : true/false value of the coreleated keyword in control.inp
                    values (dict): Values of the control variables mentioned by user in control.inp 

    """

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
            if values["restart_check"] == "true":
                 test_l = "* xyz  " + str(values["restart_charge"]) + "  " + str(multip) + " \n"
            else:
                 test_l = "* xyz  " + l1_chr_mul.strip() + " \n"
            com_f.write(test_l)
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
#            print("check if rel_file.txt exists!!")
            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write("check if rel_file.txt exists!!")

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
            uniq_atom_res = uniqatoms(sym)
            if values["ALLELE"] == "true":  ### CHECK!!!!
                for iat in range(int(uniq_atom_res["N_ua"])):
                    pre1 = uniq_atom_res["uniq_sym"]
                    at_pr1 = pre1[iat]
                    com_f.write("  NewNCore " + at_pr1 + "  " + " 0  end\n")
            else:
                for iat in range(int(uniq_atom_res["N_ua"])):
                    pre1 = uniq_atom_res["uniq_sym"]
                    at_pr1 = pre1[iat]
                    NFC_res = NFC(at_pr1)
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
            uniq_atom_res = uniqatoms(sym)
            fname =  basis
            if Nat == 1:  # check datatype
                printbas(fname, sym[0])    ## what is fname, it should be the basis set keyword. but which one??
            else:
                for iat1 in range(int(uniq_atom_res["N_ua"])):
                    pre2 = uniq_atom_res["uniq_sym"]
                    at_pr2 = pre2[iat1]
                    printbas(fname, at_pr2)
            com_f.write("end\n")
#   os.system("cat basis.com >> input.com")

    print(values["orca_dir"] + "orca input.com > input.out")
    os.system(values["orca_dir"] + "orca input.com > input.out")
    os.system("cat input.com >> ORCA_G4MP2.com")
    os.system("cat input.out >> ORCA_G4MP2.out")
####### runorca - E
