def dH_atoms(at):
    """
    Returns the enthalpy corrections of the element.

            Parameters:
                    at (char): Symbol of the element

            Returns:
                    dH_atoms (float): Enthalpy corrections of the element
    """

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

