import os, sys, string
import linecache, math
import numpy as np
import datetime , time


from pople import NFC
from pople import uniqatoms
from pople import orca_printbas


####### orca_run - S
def orca_run(method, basis,optfreq,custombasis, correlated, values, charge, multip, sym, R_coord):
    """
    Runs orca

            Parameters:
                    method (char) : Name of functional to be used
                    basis (char)  :  Basis set name
                    optfreq (char) : true/false value of the optfreq keyword 
                    custombasis (char) : true/false value of the custombasis keyword
                    correlated (char) : true/false value of the correlated keyword 
                    values (dict): Values of the control variables 

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

        Nat=len(sym)
        com_f.write("*xyz  "+str(charge)+" "+str(multip) + "\n")
        for tmp in range(Nat):
            R_x=float(R_coord[tmp][0])
            R_y=float(R_coord[tmp][1])
            R_z=float(R_coord[tmp][2])
            com_f.write(' {:2s}{:15.8f}{:15.8f}{:15.8f}\n'.format(sym[tmp],R_x,R_y,R_z))  
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
            f1 = open("rel_file.txt", "r")
            com_f.write(f1.read())
            f1.close()
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
            com_f.write("  TcutDOPre =  " + str(values["TcutDOPre"]) +"\n")   #TODO Is this really needed?
            com_f.write("end\n")

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
        if custombasis == "true":
            com_f.write("%basis \n")
        
    if custombasis == "true":
        uniq_atom_res = uniqatoms(sym)
        fname =  basis
        if Nat == 1:  
            orca_printbas(fname, sym[0]) 
        else:
            for iat1 in range(int(uniq_atom_res["N_ua"])):
                orca_printbas(fname,  uniq_atom_res["uniq_sym"][iat1])  # GTBAS1 C 
        with open("input.com", "a") as com_f:
            com_f.write("end\n")

    os.system(values["orca_exe"] + " input.com > input.out")
    os.system("cat input.com >> ORCA.inp")
    os.system("cat input.out >> ORCA.out")
    #os.system("rm -f input*")
####### orca_run - E
