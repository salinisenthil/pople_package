def At_SO(at, charge):
    """
    Returns atomic energy correction due to spin orbit coupling

            Input:
                    at (str): Atomic symbol
                    charge (int): Atomic charge

            Returns:
                    At_SO (float): Energy correction due to spin orbit coupling
    """

    S_dict_neutral = {"B" : -0.05, "C" : -0.14, "O" : -0.36, "F"  : -0.61, "Al" : -0.34, "Si" : -0.68,  "S" : -0.89, "Cl" : -1.34, "Ga" : -2.51, "Ge" : -4.41, "Se" : -4.3,  "Br" : -5.6, "Fe" : -1.84, "I" : -11.548}
    S_dict_cation  = {"C" : -0.2,  "N" : -0.43, "F" : -0.67, "Ne" :-1.19,  "Si" : -0.93, "P"  : -1.43, "Cl" : -1.68, "Ar" : -2.18, "Ge" : -5.37, "As" : -8.04, "Br" : -6.71, "Kr" : -8.16, "I" :-14.028}
    S_dict_anion   = {"B" : -0.03, "O" : -0.26, "Al": -0.28, "P"  : -0.45,  "S" : -0.88}

    At_SO = 0.0 # default

    if charge == 0:
        if at in S_dict_neutral:
            At_SO = S_dict_neutral[at]
    elif charge == +1:
        if at in S_dict_cation:
            At_SO = S_dict_cation[at]
    elif charge == -1:
        if at in S_dict_anion:
            At_SO = S_dict_anion[at]

    return(At_SO)

def Mol_SO(Nat, multip, charge, sym, SO_3rdrow_mols_val):
    """
    Returns molecular spin orbit correction

            Input:
                    Nat (int): Number of atoms in the molecule
                    multip (int): Multiplcity of the molecule
                    charge (int): Charge on the molecule
                    sym (list, str) : List of atoms in the molecule
                    SO_3rdrow_mols_val (char) : true/false value of the SO_3rdrow_mols_val keyword

            Returns:
                    Mol_SO (float): Molecular spin orbit correction
    """

    Mol_SO = 0
    
    # Special Case - Acetylene cation 2Pi state
    if Nat == 4 and multip == 2 and charge == 1:
        countH_temp =0
        countC_temp =0
        for tmp in range(len(sym)):
            if sym[tmp] == "H":
                countH_temp= countH_temp +1
            if sym[tmp] == "C":
                countC_temp = countC_temp +1
        if countH_temp == 2 and countC_temp == 2:
            Mol_SO = -0.07 
    
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("[ NOTE: Detected C2H2+: Using SO parameters for acetylene cation, Ref: J. Chem. Phys. 114 (2001) 9287 ]\n")
    # Special Case - Acetleyne - E
    
    # For diatomics with multip = 2
    if Nat == 2 and multip == 2 :

        sort_sym = sorted(sym, reverse=True)

        if SO_3rdrow_mols_val == "true":  # for 3rd_row elements
    
            if charge == 0:
                if sort_sym[0] == 'O' and sort_sym[1] == 'Br':  Mol_SO = -2.20
    
                # COMMMENT: paper has it for cation, but it looks like it is for neutral
                if sort_sym[0] == 'Se' and sort_sym[1] == 'H':  Mol_SO = -4.21
    
            if charge == +1:  
                if sort_sym[0] == 'K' and sort_sym[1] == 'Br':  Mol_SO = -2.99
                if sort_sym[0] == 'H' and sort_sym[1] == 'As':  Mol_SO = -3.54
                if sort_sym[0] == 'H' and sort_sym[1] == 'Br':  Mol_SO = -6.26
                if sort_sym[0] == 'F' and sort_sym[1] == 'Br':  Mol_SO = -6.10
                if sort_sym[0] == 'Na' and sort_sym[1] == 'Br':  Mol_SO = -3.93
                if sort_sym[0] == 'Br' and sort_sym[1] == 'Br':  Mol_SO = -6.55
    
        else:  # for non 3rd row elements, first and second rows
            if charge == 0:
                if sort_sym[0] == 'H' and sort_sym[1] == 'C':  Mol_SO = -0.07
                if sort_sym[0] == 'O' and sort_sym[1] == 'H':  Mol_SO = -0.30
                if sort_sym[0] == 'O' and sort_sym[1] == 'N':  Mol_SO = -0.27
                if sort_sym[0] == 'O' and sort_sym[1] == 'Cl':  Mol_SO = -0.61
                if sort_sym[0] == 'S' and sort_sym[1] == 'H':  Mol_SO = -1.01
                if sort_sym[0] == 'P' and sort_sym[1] == 'O':  Mol_SO = -0.53
                if sort_sym[0] == 'Si' and sort_sym[1] == 'H':  Mol_SO = -0.34
    
            if charge == -1:
                if sort_sym[0] == 'N' and sort_sym[1] == 'H':  Mol_SO = -0.12
                if sort_sym[0] == 'P' and sort_sym[1] == 'H':  Mol_SO = -0.45
                if sort_sym[0] == 'O' and sort_sym[1] == 'O':  Mol_SO = -0.34
                if sort_sym[0] == 'S' and sort_sym[1] == 'S':  Mol_SO = -1.12
    
            if charge == +1:
                if sort_sym[0] == 'H' and sort_sym[1] == 'F':  Mol_SO = -0.62
                if sort_sym[0] == 'P' and sort_sym[1] == 'H':  Mol_SO = -0.67
                if sort_sym[0] == 'H' and sort_sym[1] == 'Cl':  Mol_SO = -1.60
                if sort_sym[0] == 'N' and sort_sym[1] == 'N':  Mol_SO = -0.17
                if sort_sym[0] == 'O' and sort_sym[1] == 'O':  Mol_SO = -0.43
                if sort_sym[0] == 'P' and sort_sym[1] == 'P':  Mol_SO = -0.57
                if sort_sym[0] == 'S' and sort_sym[1] == 'S':  Mol_SO = -1.25
                if sort_sym[0] == 'Cl' and sort_sym[1] == 'Cl':  Mol_SO = -1.77
                if sort_sym[0] == 'F' and sort_sym[1] == 'Cl':  Mol_SO = -1.60
    
    return(Mol_SO)
