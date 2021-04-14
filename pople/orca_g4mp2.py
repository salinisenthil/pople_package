import os, sys, string
import linecache, math
import numpy as np
import datetime , time

from pople import au2ev, au2kcm, au2kjm, kcm2au, j2au, au2cmi, cmi2au, kB, h, Ry, c, T, e, N_avo, au2j

from pople import At_SO, Mol_SO
from pople import principal_coord
from pople import orca_run
from pople import thermal, HLC_g4mp2
from pople import nanb
from pople import print_header
from pople import start_time
from pople import final_time

####### orca_g4mp2 - S
def orca_g4mp2(values, start_time_main): 
    """
    Runs g4mp2

            Parameters:
                    values (dict): Values of the control variables 
                    start_time_main (float) :  Initial start time of the program pople.run()

    """
    print_header()
    start_time()
    print("\n Printing detailed output in Thermochemistry.out\n")
    print(" Opt/Freq Started: Step 1/5")

    if values["isatom"] != "true":
        # CONTAINS INPUT GEOMETRY
        new_i_xyz = open("inp.xyz", "r")
        sym =[]
        Rlist = [] 
        iline=0
        for l in new_i_xyz:
            lsplit=l.split()
            if iline == 0: Nat = int(l.strip())
            if iline == 1: 
                charge=eval(lsplit[0])
                multip=eval(lsplit[1])
            if iline > 1:
                sym.append(lsplit[0])
                Rlist.append(lsplit[1])
                Rlist.append(lsplit[2])
                Rlist.append(lsplit[3])
            iline=iline+1
        new_i_xyz.close()
 
        R_coord = np.resize(np.array(Rlist),[Nat,3])

        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write(" *  Input atomic coordinates (Angstroem)  *\n\n")
            for tmp in range(Nat): 
                R_x=float(R_coord[tmp][0])
                R_y=float(R_coord[tmp][1])
                R_z=float(R_coord[tmp][2])
                ther_chem.write(' {:2s}{:15.8f}{:15.8f}{:15.8f}\n'.format(sym[tmp],R_x,R_y,R_z))  
            ther_chem.write("\n")

        if values["FROZEN_GEOM"] == "true":

            os.system("mv inp.xyz  input.xyz")

            #==== Principal coordinates transformation
            if Nat > 2:
                Ievals = principal_coord() # check
                linnonlin = "NL"
                if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                    linnonlin = "L"
            if Nat == 2:    linnonlin = "L"

            freq_file = open("freq.txt", "r")
            freq =[]
            for l in freq_file:
                lsplit=l.split()
                freq.append(lsplit[0])
            freq_file.close()
            nfreq=len(freq)

            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write(" *  Input unscaled harmonic wavenumbers (cm^-1)  *\n\n")
                tmp = 0
                for tmp in range(nfreq): 
                    ther_chem.write('{:6d}{:15.4f}'.format(tmp+1,float(freq[tmp])))
                    ther_chem.write("\n")
    
                ther_chem.write("\n")
                ther_chem.write(" *  Performing single point calculation using input geometry and frequencies  *\n")

        else:
            start_time_new = time.time()

            orca_run(values["method_opt_freq"], values["basis_opt_freq"], "true", values["custombasis_opt_freq"], "false", values, charge, multip, sym, R_coord)

            #==== Principal coordinates transformation
            if Nat > 2:
                Ievals = principal_coord() # check
                linnonlin = "NL"
                if abs(Ievals[0]) < 1*(10**-10):  ## CHECK esp Ievals[0 or 1]
                    linnonlin = "L"
            if Nat == 2:    linnonlin = "L"

            # CONTAINS OPTIMIZED GEOMETRY
            new_i_xyz = open("input.xyz", "r")
            Rlist = [] 
            iline=0
            for l in new_i_xyz:
                if iline > 1:
                    lsplit=l.split()
                    Rlist.append(lsplit[1])
                    Rlist.append(lsplit[2])
                    Rlist.append(lsplit[3])
                iline=iline+1
            new_i_xyz.close()
 
            R_coord = np.resize(np.array(Rlist),[Nat,3])

            # Nat, charge-multiplicity are from inp.xyz
            # sym and coordinates from input.xyz
            with open("temp.xyz", "w") as t_xyz:
                t_xyz.write(str(Nat) + "\n")
                t_xyz.write(str(charge)+" "+str(multip) + "\n")
                for tmp in range(Nat): 
                    R_x=float(R_coord[tmp][0])
                    R_y=float(R_coord[tmp][1])
                    R_z=float(R_coord[tmp][2])
                    t_xyz.write(' {:2s}{:15.8f}{:15.8f}{:15.8f}\n'.format(sym[tmp],R_x,R_y,R_z))  # writes atoms + coord as from geo,xyz

            os.system("mv temp.xyz  inp.xyz") 
        
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
            
            iopt = 999
            if values["isatom"] != "true":
                os.system("grep 'OPTIMIZATION RUN DONE' input.out | wc | awk '{print $1}' > scr_lc1")
            iopt = int(linecache.getline("scr_lc1",1).strip())
            
            os.system("rm -f input* ")

            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write(" *  Optimized atomic coordinates (Angstroem)  *\n\n")
                for tmp in range(Nat): 
                    R_x=float(R_coord[tmp][0])
                    R_y=float(R_coord[tmp][1])
                    R_z=float(R_coord[tmp][2])
                    ther_chem.write(' {:2s}{:15.8f}{:15.8f}{:15.8f}\n'.format(sym[tmp],R_x,R_y,R_z))  
                ther_chem.write("\n")
    
                ther_chem.write(" *  Unscaled harmonic wavenumbers (cm^-1)  *\n\n")
                tmp = 0
                for tmp in range(nfreq): 
                    ther_chem.write('{:6d}{:15.4f}'.format(tmp+1,freq[tmp]))
                    ther_chem.write("\n")
    
                ther_chem.write("\n")
                ther_chem.write(" *  Scaled harmonic wavenumbers (cm^-1), scaling factor used is " + str(values["scalfac"]) + "  *\n\n")
                tmp = 0
                for tmp in range(nfreq): 
                    new_f = freq[tmp] * float(values["scalfac"])
                    ther_chem.write('{:6d}{:15.4f}'.format(tmp+1,new_f))
                    ther_chem.write("\n")
    
                ther_chem.write("\n")
                end_time       = time.time()
                gf1 = (f""" > Opt+Freq done in{end_time - start_time_new:20.2f} s
 > Elapsed time =  {end_time - start_time_main:20.2f} s
""")
                ther_chem.write(gf1)
    
            if iopt == 0:
                with open("Thermochemistry.out", "a") as ther_chem:
                    ther_chem.write(" *****************************************\n")
                    ther_chem.write(" *  ERROR: Geometry Optimization Failed  *\n")
                    ther_chem.write(" *****************************************")

        #=== common for both frozen_geom and relaxed
        if linnonlin == "NL":
            if nfreq != ((3*Nat)-6):
                with open("Thermochemistry.out", "a") as ther_chem:
                    nl1 = """

 **************************************************************
 * WARNING: No. of freq. != 3N - 6 for a non-linear molecule! *
 **************************************************************

"""
        elif linnonlin == "L":
            if nfreq != ((3*Nat)-5):
                with open("Thermochemistry.out", "a") as ther_chem:
                    nl1 = """

 ***********************************************************
 * WARNING: No. of freq. != 3N - 5  for a linear molecule! *
 ***********************************************************

"""
                    ther_chem.write(nl1)
    else: #=== atoms
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write(" [ NOTE: Geometry optimization and frequency calculations are skipped for atom ]\n")

    if values["isatom"] == "true":
        Nat=1
        with open("inp.xyz", "r") as new_i_xyz:
           linecache.clearcache()
           title = linecache.getline("inp.xyz",2).strip()
           tl1 = title.split()
           charge = int(tl1[0])
           multip = int(tl1[1])
           sym = []
           Rlist=[]
           iline=0
           for l in new_i_xyz:
               if iline > 1:
                   lsplit=l.split()
                   sym.append(lsplit[0])
                   Rlist.append(lsplit[1])
                   Rlist.append(lsplit[2])
                   Rlist.append(lsplit[3])
               iline=iline+1
           R_coord = np.resize(np.array(Rlist),[Nat,3])
 
    print(" Opt/Freq Completed")
# === SINGLE POINT CCSD(T)/S
    start_time_new = time.time()
    print(" CCSD(T)  Started: Step 2/5")
    if values["Ntotale"] > 0:
        if Nat == 1:
            if values["basis_ccsdt"]  == "GTBAS1": values["basis_ccsdt"] = 'GTBAS1atm'

        if values["restart_cc"] == "true":
            values["restart_check"] = "true"
            orca_run("HF ", values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "false", values, charge+2, multip, sym, R_coord)
    
        if values["DLPNO_CCSDT"] == "true":
            values["switch_DLPNO_CCSDT"] = "true"
        
        orca_run(values["method_ccsdt"], values["basis_ccsdt"], "false", values["custombasis_ccsdt"], "true", values, charge, multip, sym, R_coord)
        end_time= time.time()
        
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc2")
            tmp = linecache.getline("scr_lc2",1).strip()
            E_ccsdt = float(tmp)
            E_mp2 = E_ccsdt
        else:
            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc2")
            tmp = linecache.getline("scr_lc2",1).strip()
            E_ccsdt = float(tmp)
            if values["DLPNO_CCSDT"] == "true":
                os.system("grep 'MP2 TOTAL ENERGY' input.out  | awk '{print $4}' > scr_lc20") # not same as 'Initial E(tot)' for open-shell systems, ex. 'C' atom
            else:
                os.system("grep 'Initial E(tot)' input.out | awk '{print $4}' > scr_lc20")
            tmp = linecache.getline("scr_lc20",1).strip()
            E_mp2 = float(tmp)
    else:
        E_ccsdt = 0.0
        E_mp2 = 0.0

    os.system("rm -f input* ")

    end_time = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        print(" CCSD(T)  Completed")
        cc1 = (f"""
 > CCSD(T) done in {end_time - start_time_new:20.2f} s
 > Elapsed time =  {end_time - start_time_main:20.2f} s
""")
        ther_chem.write(cc1)
    
    if values["DLPNO_CCSDT"] == "true":
        values["switch_DLPNO_CCSDT"] = "false"

# === SINGLE POINT MP2L
    start_time_new = time.time()  # can this be moved into the if statement?
    print(" MP2/L    Started: Step 3/5")
    if (values["Ntotale"] > 0) and (values["ccsdt_cbs"] != "true") : ## where is ccsdt_cbs declared?
        if values["restart_mp2"] == "true":
            values["restart_check"] = "true"
            orca_run("HF ", values["basis_mp2"], "false", values["custombasis_mp2"], "false", values, charge+2, multip, sym, R_coord)

        orca_run(values["method_mp2"], values["basis_mp2"], "false", values["custombasis_mp2"], "true", values, charge, multip, sym, R_coord)
    
        if (values["Ntotale"] == 1) or (values["Ntotale"] - values["Ntotalecore"] <= 1):
            os.system("grep 'Total Energy ' input.out  | awk '{print $4}' > scr_lc6")
            tmp = linecache.getline("scr_lc6",1).strip()
            E_mp2L = float(tmp)
            E_hfL = E_mp2L
        else:
            os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc6")
            tmp = linecache.getline("scr_lc6",1).strip()
            E_mp2L = float(tmp)
            if values["flag_RIMP2"] == "true":
                os.system("grep 'RI-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            elif values["flag_DLPNOMP2"] == "true":
                os.system("grep 'DLPNO-MP2 CORRELATION ENERGY:' input.out | awk '{print $4}' > scr_lc7")
            else:
                os.system("grep 'MP2 CORRELATION ENERGY ' input.out | awk '{print $5}' > scr_lc7")
            tmp = linecache.getline("scr_lc7",1).strip()
            E_hfL = E_mp2L - float(tmp)

        os.system("rm -f input* ")
    else:
        E_hfL= 0
        E_mp2L = 0

    end_time = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        print(" MP2/L    Completed")
        mp12 = (f"""
 > MP2/L done in   {end_time - start_time_new:20.2f} s
 > Elapsed time =  {end_time - start_time_main:20.2f} s
""")
        ther_chem.write(mp12)

#== SINGLE POINT HF/VTZ
    start_time_new = time.time()
    print(" HF/VTZ   Started: Step 4/5")
    if (values["Ntotale"] > 0) and (values["ccsdt_cbs"] != "true"):
        if values["restart_hf3"] == "true":
            values["restart_check"] = "true"
            orca_run("HF ", values["basis_hf3"], "false", values["custombasis_hf3"], "false", values, charge+2, multip, sym, R_coord)

        orca_run(values["method_hf3"], values["basis_hf3"], "false", values["custombasis_hf3"], "false", values, charge, multip, sym, R_coord)

        os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc8")
        tmp = linecache.getline("scr_lc8",1).strip()
        E_hfT = float(tmp)
        
        os.system("rm -f input* ")
    else:
        E_hfT = 0
    
    end_time = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        print(" HF/VTZ   Completed")
        hf11 = (f"""
 > HF/VTZ done in  {end_time - start_time_new:20.2f} s
 > Elapsed time =  {end_time - start_time_main:20.2f} s
""")
        ther_chem.write(hf11)

# == SINGLE POINT HF/VQZ
    start_time_new = time.time()
    print(" HF/VQZ   Started: Step 5/5")
    if (values["Ntotale"] > 0 ) and (values["ccsdt_cbs"] != "true"):
        if values["restart_hf4"] == "true":
            values["restart_check"] = "true"
            orca_run("HF ", values["basis_hf4"], "false", values["custombasis_hf4"], "false", values, charge+2, multip, sym, R_coord)

        orca_run(values["method_hf4"], values["basis_hf4"], "false", values["custombasis_hf4"], "false", values, charge, multip, sym, R_coord)

        os.system("grep 'FINAL SINGLE POINT ENERGY ' input.out  | awk '{print $5}' > scr_lc9")
        tmp = linecache.getline("scr_lc9",1).strip()
        E_hfQ = float(tmp)
        
        os.system("rm -f input* ")
    else:
        E_hfQ = 0
    
    end_time = time.time()
    with open("Thermochemistry.out", "a") as ther_chem:
        print(" HF/VQZ   Completed")
        hf12 = ( f"""
 > HF/VQZ done in  {end_time - start_time_new:20.2f} s
 > Elapsed time =  {end_time - start_time_main:20.2f} s
""")
        ther_chem.write(hf12)

    #=== Thermal
    if values["isatom"] == "true":
        freq=[]
        linnonlin=""
    dE_ZPE, dE_vib, dE_rot, dE_tra, dE_thermal = thermal(values["isatom"], freq, values["scalfac"],linnonlin, 298.15)
    
    #=== SO
    dE_SO = 0.0
    if values["isatom"] == "true":
        dE_SO = At_SO(sym[0], charge)
    else:
        dE_SO = Mol_SO(Nat, multip, charge, sym, values["SO_3rdrow_mols"])

    dE_SO = dE_SO / 1000.0
    
    #=== HLC
    HLC_params = [  values["AA"],  values["BB"], values["CC"],values["DD"], values["ApAp"], values["EE"]  ]

    with open("Thermochemistry.out", "a") as ther_chem:
        cs1 = ( f""" 
 *  HLC parameters in milli hartree  *
 
 A  = {values["AA"]:.5f} 
 A' = {values["ApAp"]:.5f}
 B  = {values["BB"]:.5f} 
 C  = {values["CC"]:.5f} 
 D  = {values["DD"]:.5f} 
 E  = {values["EE"]:.5f} 
""")
        ther_chem.write(cs1)

    if values["HLCeqZERO"] == "true":
        dE_HLC = 0.0
    else:
        dE_HLC = HLC_g4mp2(charge, multip, sym, values["isatom"], HLC_params)

    values["HLC0"] = dE_HLC

# CBS extrapolation hard-coded for HF (AVQZ, AV5Z)
    if values["HF_CBS_default"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-1.63)) / (1 - math.exp(-1.63))
    if values["HF_CBS_orca_23_def2"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-10.39*(math.sqrt(3.0)- math.sqrt(2.0)))) / (1 - math.exp(-10.39*(math.sqrt(3.0)-math.sqrt(2.0))))
    if values["HF_CBS_orca_34_def2"] == "true": E_HF_limit = (E_hfQ - E_hfT * math.exp(-7.88*( math.sqrt(4.0)- math.sqrt(3.0)))) /  (1 - math.exp(-7.88*(math.sqrt(4.0)-math.sqrt(3.0))))
    if values["HF_CBS_orca_23_cc"] == "true":   E_HF_limit = (E_hfQ - E_hfT * math.exp(-4.42*( math.sqrt(3.0)- math.sqrt(2.0)))) / (1 - math.exp(-4.42* (math.sqrt(3.0)-math.sqrt(2.0))))
    if values["HF_CBS_orca_34_cc"] == "true":   E_HF_limit = (E_hfQ - E_hfT * math.exp(-5.46*( math.sqrt(4.0)- math.sqrt(3.0)))) / (1 - math.exp(-5.46* (math.sqrt(4.0)-math.sqrt(3.0))))
    
    dE_HF_CBS   = E_HF_limit - E_hfL
    dE_MP2_CBS  = E_mp2L - E_mp2
    
    if values["ccsdt_cbs"] == "true":
        dE_HF_CBS = 0
        dE_MP2_CBS = 0
        eps0 = E_ccsdt +  dE_HLC + dE_SO
    else:
        eps0 = E_ccsdt + dE_HF_CBS + dE_MP2_CBS + dE_HLC + dE_SO
    
    if values["G4MP2TM"] == "true":
        eps0 = eps0 + E_rel

    values["U0"] = eps0 + dE_ZPE
    values["UT"] = eps0 + dE_ZPE + dE_thermal
    values["HT"] = eps0 + dE_ZPE + dE_thermal + kB * T * j2au

    values["breakdown"] = [ E_ccsdt, dE_MP2_CBS, dE_HF_CBS, dE_HLC, dE_SO, dE_ZPE, (dE_vib - dE_ZPE), dE_rot, dE_tra, dE_thermal, (kB * T * j2au)  ] # 11 columns

    with open("Thermochemistry.out", "a") as ther_chem:
        end_time = time.time()

        nl2 = (f"""
           Temperature ={T:20.8f} kelvin                     
              Pressure =          1.00000000 atm

            E(CCSD(T)) ={E_ccsdt:20.8f} hartree            
               dE(MP2) ={dE_MP2_CBS:20.8f} hartree     <-- post-HF basis set correction                      
                dE(HF) ={dE_HF_CBS:20.8f} hartree     <-- HF-level basis set correction
                dE(SO) ={dE_SO:20.8f} hartree     <-- spin-orbit correction
               dE(HLC) ={dE_HLC:20.8f} hartree     <-- higher-level correction
               dE(ZPE) ={dE_ZPE:20.8f} hartree     <-- zero-point correction
           dE(Thermal) ={dE_thermal:20.8f} hartree     <-- thermal correction

     Electronic Energy ={E_ccsdt+dE_MP2_CBS+dE_HF_CBS+dE_SO+dE_HLC:20.8f} hartree     <-- E = E(CCSD(T)) + dE(MP2) + dE(HF) + dE(SO) + dE(HLC) 
  Internal Energy(0 K) ={values["U0"]:20.8f} hartree     <-- U0 = E + E(ZPE)
  Internal Energy(T K) ={values["UT"]:20.8f} hartree     <-- UT = U0 + E(Thermal)
         Enthalpy(T K) ={values["HT"]:20.8f} hartree     <-- HT = UT + kB*T       

 > Enthalpy calculation done
 > Elapsed time =  {end_time - start_time_main:20.2f} s

""")

        ther_chem.write(nl2) 
    os.system("rm -f scr*")

    print(" --- Terminating ---\n\n")

####### orca_g4mp2 - E
