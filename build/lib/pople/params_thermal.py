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

#def HOF_atoms(at, values):
def HOF_atoms(at):
    """
    Returns the heat of formation of the element.

            Parameters:
                    at (char): Symbol of the element
                    values (dict): Values of the control variables 

            Returns:
                    HOF_atoms (float): Heat of formation of the element
    """

    h = 6.626070040*(10**-34)    #6.626070040d-34
    Ry = 10973731.568508  #10973731.568508d0
    c = 299792458    #299792458d0
    N_avo = 6.02214179000000*(10**+23)  #6.02214179000000d+23
    au2kcm = 2 * Ry * h * c * N_avo / 4184
    kcm2au = 1 / au2kcm

    hof_dict_1 = {"H": 51.63,   "Li": 37.69,  "Be": 76.48,   "B": 136.2,  "C": 169.98, \
                  "N": 112.53,  "O": 58.99 ,  "F": 18.47,    "Na": 25.69, "Mg": 34.87, \
                  "Al": 78.23,  "Si": 106.6,  "P": 75.42,    "S": 65.66,  "Cl": 28.59, \
                  "K": 21.48303059,  "Ca": 42.38503824,   "Fe": 98.7373,  "Ga": 64.763384321, \
                  "Ge": 89.354,   "As": 68.86,   "Se": 57.8991,   "Br": 28.1836,   "I": 25.62}

#   if values["hof_C_hydrocarbon"]=="true":
#       new_val = {"C": 170.11}
#       hof_dict_1.update(new_val)

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

