def Mol_SO(Nat, multip, charge, sym, SO_3rdrow_mols_val):   # number of atoms, multiplicity, charge, array of atoms in molecule, value of SO_3rdrow_mols (from orca.inp file)
    """
    Returns the molecular spin orbit correction.

            Parameters:
                    Nat (int): Number of atoms in the molecule
                    multip (int): Multiplcity of the molecule
                    charge (int): Charge on the molecule
                    sym (list, char) : List of atoms in the molecule
                    SO_3rdrow_mols_val (char) : true/false value of the SO_3rdrow_mols_val keyword

            Returns:
                    Mol_SO (float): Molecular spin orbit correction
    """

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
    
            if charge == +1:  
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
