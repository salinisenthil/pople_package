def orca_g4mp2_ctrl(user_method_type):
    """
    Returns the default control variables/values of the chosen method
            Parameters:
                    user_method_type (char): Method chosen by user in pople.inp
                    
            Returns:
                    val (dict): Dictionary of all control variables required for pople run
    """
    val = {}
#### IMPORTANT: use only "true"  or "false" here -- this section is case sensitive  ### THANKS!
    if user_method_type == "g4mp2-xp" or user_method_type == "G4MP2-XP" or user_method_type == "G4MP2-xp" or user_method_type == "g4mp2-XP":

        val = { "method_opt_freq" : "B3LYP/G Grid7 RIJCOSX def2/J GridX9 Printbasis ", 
        "basis_opt_freq" : "GTBAS3",
        "custombasis_opt_freq" : "true",

        "AA" : 9.51842 ,
        "BB" : 3.19566  ,
        "CC" : 9.80448 ,
        "DD" : 2.05800 ,
        "ApAp" : 9.79309 ,
        "EE" : 2.38393 ,

        "maxcore_mb": "16000",
        "nproc" : "1",
        "String_Opt": "TightOpt",

        "MGGA" : "false",

        "FROZEN_GEOM" : "false",

        "method_ccsdt" : "DLPNO-CCSD(T) RI-MP2 def2-SVP/C RIJK def2/JK TightPNO Printbasis ",
        "basis_ccsdt" : "GTBAS1",
        "custombasis_ccsdt" : "true",
        "DLPNO_CCSDT" : "true",

        "method_mp2_s" : "RI-MP2 def2-SVP/C RIJK def2/JK Printbasis ",
        "basis_mp2_s" : "GTBAS1",
        "custombasis_mp2_s" : "true",

        "method_mp2" :  "RI-MP2 def2-TZVP/C RIJK def2/JK Printbasis ",
        "basis_mp2" : "GTMP2largeXP",
        "custombasis_mp2" : "true",
        "flag_RIMP2" : "true",
        "flag_DLPNOMP2" : "false",

        "method_hf3" : "HF RIJK def2/JK Printbasis ",
        "basis_hf3" : "GFHFB3",
        "custombasis_hf3" : "true",

        "method_hf4" : "HF RIJK def2/JK Printbasis ",
        "basis_hf4" : "GFHFB4",
        "custombasis_hf4" : "true",

        "scalfac" : "0.9854",

        "conv_scf" : "VeryTight",

        "HLCeqZERO" : "false",

        "SOSCF" : "false",
        "SCFDIIS" : "false",

        "SO_3rdrow_mols" : "true",

        "LSHIFT" : "false",

        "optdiis" : "false",

        "HF_CBS_default"     : "true",
        "HF_CBS_orca_23_def2" : "false",
        "HF_CBS_orca_34_def2" : "false",
        "HF_CBS_orca_23_cc"   : "false",
        "HF_CBS_orca_34_cc"   : "false",

        "ALLELE": "false",
        "switch_load_rel_file": "false",
        "switch_guess": "false",
        "FROZEN_GEOM": "false",
        "SO_3rdrow_mols": "true",
        "iterhess": 5,
        "TcutDOPre": 3e-2,
        "HF_CBS_default":  "true",  ## ask!!!
        "HF_CBS_orca_23_def2": "false",
        "HF_CBS_orca_34_def2": "false",
        "HF_CBS_orca_23_cc": "false",
        "HF_CBS_orca_34_cc": "false",
        "switch_DLPNO_CCSDT": "false",
        "G4MP2TM": "false",
        "ccsdt_cbs": "false",
        "restart_check": "false",
        "switch_RIMP2_small": "false",
        "verticalIP": "false",
        "verticalEA": "false",
        "IPss": "false",
        "hof_C_hydrocarbon": "false",

        "restart_cc"  : "false",
        "restart_mp2" : "false",
        "restart_hf3" : "false",
        "restart_hf4" : "false"
        }
        
    if user_method_type == "g4mp2" or user_method_type == "G4MP2":
        val = { 
        "method_opt_freq"      : "B3LYP/G Grid5 " ,
        "basis_opt_freq"       : "GTBAS3" ,
        "custombasis_opt_freq" : "true" ,

        "AA" : 9.472 ,
        "BB" : 3.102 , 
        "CC" : 9.741 , 
        "DD" : 2.115 , 
        "ApAp" : 9.769 , 
        "EE" : 2.379 ,

        "maxcore_mb": "16000",
        "nproc" : "1",

        "String_Opt" : "TightOpt" ,
    
        "MGGA" : "false" ,
    
        "FROZEN_GEOM" : "false" ,
    
        "method_ccsdt" : "CCSD(T) " ,
        "basis_ccsdt" : "GTBAS1" ,
        "custombasis_ccsdt" : "true" ,
        "DLPNO_CCSDT" : "false" ,
        "switch_RIMP2_small" : "false" ,
    
#       "switch_read_RIMP2_small" : "false" ,
    
        "method_mp2_s" : "MP2 " ,
        "basis_mp2_s" : "GTBAS1" ,
        "custombasis_mp2_s" : "true" ,
    
        "method_mp2" : "MP2 " ,
        "basis_mp2" : "GTMP2largeXP" ,
        "custombasis_mp2" : "true" ,
        "flag_RIMP2" : "false" ,
        "flag_DLPNOMP2" : "false" ,
    
        "method_hf3" : "HF " ,
        "basis_hf3" : "GFHFB3" ,
        "custombasis_hf3" : "true" ,
    
        "method_hf4" : "HF " ,
        "basis_hf4" : "GFHFB4" ,
        "custombasis_hf4" :"true" ,
    
        "scalfac" : "0.9854" ,
    
        "conv_scf" : "VeryTight" ,
    
        "HLCeqZERO" : "false" ,
    
        "SOSCF" : "false" ,
        "SCFDIIS" : "false" ,
    
        "SO_3rdrow_mols" : "true" ,
    
        "LSHIFT" : "false" ,
    
        "optdiis" : "false" ,
    
        "HF_CBS_default"      : "true" ,
        "HF_CBS_orca_23_def2" : "false" ,
        "HF_CBS_orca_34_def2" : "false" ,
        "HF_CBS_orca_23_cc"   : "false" ,
        "HF_CBS_orca_34_cc"   : "false" ,

        "ALLELE": "false",
        "switch_load_rel_file": "false",
        "switch_guess": "false",
        "FROZEN_GEOM": "false",
        "SO_3rdrow_mols": "true",
        "iterhess": 5,
        "TcutDOPre": 3e-2,
        "HF_CBS_default": "true",  ## ask!!!
        "HF_CBS_orca_23_def2": "false",
        "HF_CBS_orca_34_def2": "false",
        "HF_CBS_orca_23_cc": "false",
        "HF_CBS_orca_34_cc": "false",
        "switch_DLPNO_CCSDT": "false",
        "G4MP2TM": "false",
        "ccsdt_cbs": "false",
        "restart_check": "false",
        "switch_RIMP2_small": "false",
        "verticalIP": "false",
        "verticalEA": "false",
        "IPss": "false",
        "hof_C_hydrocarbon": "false",
    
        "restart_cc"  : "false" ,
        "restart_mp2" : "false" ,
        "restart_hf3" : "false" ,
        "restart_hf4" : "false" 
        }

    return(val)
