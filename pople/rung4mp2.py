import os, sys, string
import linecache, math
import numpy as np
import datetime , time


from pople import Mol_SO
from pople import At_SO
from pople import principal_coord
from pople import runorca


####### rung4mp2 - S
def rung4mp2(values, start_time_main): 
    """
    Runs g4mp2

            Parameters:
                    values (dict): Values of the control variables mentioned by user in control.inp
                    start_time_main (float) :  Initial start time of the program pople.run()

    """
 
    if values["isatom"] == "true": Nat=1
    os.system("cat inp.xyz")


    if values["isatom"] != "true":
        if values["FROZEN_GEOM"] == "true":
            start_time1 = time.time()
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
                    full_l1 = l1_xyz_1.strip().split()
                    full_l = full_l1[0] + "     "+str(round(float(full_l1[1]),8)) + "     "+str(round(float(full_l1[2]),8)) + "     "+str(round(float(full_l1[3]),8))
                    sym_coord.append(full_l)  ## contains atom and coords
                    sym.append(l1_lsp1[0])
                    Rlist.append(l1_lsp1[1])
                    Rlist.append(l1_lsp1[2])
                    Rlist.append(l1_lsp1[3])
                    conv = np.array(Rlist)
                    R_coord = np.resize(conv,[Nat,3])
                freq =[]
                for tmp_j in range(Nat+3, num_l_fg+1):
                    tmp_line = float(linecache.getline("read_geom_freq.dat",tmp_j).strip())
                    freq.append(tmp_line)

                print("freq check")
                print(freq)

            if Nat > 2:
                Ievals = principal_coord() # check
                linnonlin = "NL"
                print("FROZEN_GEOM", FROZEN_GEOM)
                print("Eigenvalues of Moment of inertia:", Ievals)
                if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                    linnonlin = "L"
            if Nat == 2:    linnonlin = "L"

            if linnonlin == "L":    nfreq = 3 * Nat - 5
            if linnonlin == "NL":   nfreq = 3 * Nat - 6


        else:
            start_time1 = time.time()
            print("==check isatom==",values["isatom"])
            runorca.runorca(values["method_opt_freq"], values["basis_opt_freq"], "true", values["custombasis_opt_freq"], "false", values)
            end_time1 = time.time()

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
                    sym = []
                    Rlist = []
                    sym_coord = []
                    for tmp_l in range(3,num_l_in+1):    ## check range
                        l1_xyz_1 = linecache.getline("input.xyz",tmp_l)
                        l1_lsp1 = l1_xyz_1.split()
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

############################################
            if Nat > 2:
                Ievals = principal_coord() # check
                linnonlin = "NL"
                if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                    linnonlin = "L"
            if Nat == 2:    linnonlin = "L"

#############################################
        
            os.system("nfreq=$( grep 'The total number of vibrations considered is ' input.out  | tail -1 | awk '{print $8}' > scr_lc0);  \
                 nfreq=$( cat scr_lc0); grep -A$(( $nfreq+4 )) 'IR SPECTRUM' input.out | tail -$nfreq | awk '{print $2}' > freq.txt")

    ### need to extract nfreq from scr file
            tr = linecache.getline("scr_lc0",1).strip()
            lsp111 = tr.split()
            nfreq = int(lsp111[0])

            if linnonlin == "NL":
                if nfreq > (3*Nat)-6:
                    nfreq = (3*Nat)-6
                    os.system("nfreq=$( cat scr_lc0); grep -A$(( $nfreq+4 )) 'IR SPECTRUM' input.out | tail -$(( "+str(nfreq)+" )) | awk '{print $2}' > freq.txt")

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
            for tmp in range(Nat-1):  ## CHECK range
                ther_chem.write(str(sym_coord[tmp]) + "\n")  # writes atoms + coord as from geo,xyz
            ther_chem.write("\n")

            ther_chem.write(" Unscaled harmonic wavenumbers (cm^-1)\n")
            tmp = 0
            for tmp in range(nfreq): # check range
                ther_chem.write(str( freq[tmp] ))  
                ther_chem.write("\n")

            ther_chem.write("\n")
            ther_chem.write(" Scaled harmonic wavenumbers (cm^-1)\n")
            tmp = 0
            for tmp in range(nfreq): # check range
                ther_chem.write(str( freq[tmp] * float(values["scalfac"]) ) )  
                ther_chem.write("\n")

            ther_chem.write("\n")
            ther_chem.write(" Scaling factor used:  " + str(values["scalfac"]) + "\n\n")
            end_time2 = time.time()
            ther_chem.write(" * Geometry optimization/frequencies done in              " + str(end_time2 - start_time1) + " s\n")
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
        sym = []
        Rlist = []
        sym_coord = []
        for tmp_l in range(3,num_l_in+1):    ## check range
            l1_xyz_1 = linecache.getline("inp.xyz",tmp_l)
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

    else:
        E_ccsdt = 0
        E_mp2 = 0

    if values["G4MP2TM"] == "true":
        values["switch_load_rel_file"] = "true"
        values["switch_guess"] = "true"
        values["ALLELE"] = "true"

        if values["Ntotale"] > 0:
            if (values["Notale"] < 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1 ):
                tmp_method = "HF "
                if Nat == 1:
                    if values["basis_ccsdt"] == "GTBAS1": values["basis_ccsdt"] = "GTBAS1atm"
                    runorca.runorca(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                    #print(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "--1--")
                else:  # ASK
                    runorca.runorca(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                    #print(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "--2--")
            elif (values["Ntotale"] == 2) or (values["Ntotale"] - values["Ntotalecore"] == 1 ):
                if Nat == 1:
                    if values["basis_ccsdt"] == "GTBAS1": 
                        values["basis_ccsdt"] = "GTBAS1atm"  # check
                        runorca.runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                        #print(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "--3--")
                    else:
                        runorca.runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                        #print(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "--4--")
            else:
                if values["restart_cc"] == "true":
                    #charge = charge + 2
                    tmp_method = "HF "
                    values["restart_check"] = "true"
                    values["restart_charge"] = charge + 2
                    runorca.runorca(tmp_method, values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values)
                    #charge = charge - 2
                    values["restart_charge"] = values["restart_charge"] - 2

                runorca.runorca(values["method_ccsdt_rel"], values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "true", values)
                #print(tmp_method, values["basis_ccsdt_rel"], "false", values["custombasis_ccsdt"], "--5--")
            end_time1 = time.time()

            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc3")
            tmp = linecache.getline("scr_lc3",1).strip()
            E_ccsdt_rel = float(tmp)
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
        #print("block dlpno ccsd check", values["switch_read_RIMP2_small"], values["DLPNO_CCSDT"], "inside")
        if values["switch_read_RIMP2_small"] != "true":
            #print( multip , values["Ntotale"], values["Ntotalecore"], "ntotale ntotalecore")
            if (multip != 1) or (values["switch_RIMP2_small"] == "true") :
                start_time1 = time.time()
                #print(values["Ntotale"], values["Ntotalecore"], "ntotale ntotalecore")
                if (values["Ntotale"] == "1") or (values["Ntotale"] - values["Ntotalecore"] <= 1) :
                    tmp_method = "HF "
                    #print("block dlpno ccsd check  1")
                    runorca.runorca(values["method_mp2_s"], values["basis_mp2_s"], "false", values["custombasis_mp2_s"], "true", values)
                else:
                    #print("block dlpno ccsd check  2")
                    runorca.runorca(values["method_mp2_s"], values["basis_mp2_s"], "false", values["custombasis_mp2_s"], "true", values)
                end_time1= time.time()

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
            #print("=== E_hfL === 00", E_hfL, E_mp2L, values["Ntotale"], values["Ntotalecore"])
        else:
            if values["flag_RIMP2"] == "true":
                os.system("grep 'RI-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            elif values["flag_DLPNOMP2"] == "true":
                os.system("grep 'DLPNO-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            else:
                os.system("grep 'MP2 CORRELATION ENERGY ' input.out | awk '{print $5}' > scr_lc7")
            tmp = linecache.getline("scr_lc7",1).strip()
            E_hfL = float(tmp)

            #print("=== E_hfL ===", E_hfL, E_mp2L, values["Ntotale"], values["Ntotalecore"])
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
        HLC_SO = Mol_SO(Nat, multip, charge, sym, values["SO_3rdrow_mols"])
        #print(HLC_SO)
    else:
        HLC_SO = At_SO(sym[0], charge)
    
    #print("==check near HLC ",values["isatom"],Nat)
    if values["HLCeqZERO"] == "true":
        HLC = 0
    else:
        if values["isatom"] != "true":
            if values["Na"] == values["Nb"]:
                #print("HLC atoms-1:",values["Na"],values["Nb"])
                HLC = - values["AA"] * values["Nb"]
            else:
                #print("HLC atoms-2:",values["Na"],values["Nb"])
                HLC = - values["ApAp"] * values["Nb"] - values["BB"] * (values["Na"]- values["Nb"])   # check
        else:
             #print("HLC atoms-3:",values["Na"],values["Nb"])
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
