def atno(at):
    """
    Returns atom number

            Input:
                    at (str): Atomic symbol

            Returns:
                    atno (int): Atom number 
    """
    at_num_dict ={"H":  1,  "He":  2,  "Li":  3,  "Be":  4,   "B":  5,\
                  "C":  6,   "N":  7,   "O":  8,   "F":  9,  "Ne": 10,\
                 "Na": 11,  "Mg": 12,  "Al": 13,  "Si": 14,   "P": 15,\
                  "S": 16,  "Cl": 17,  "Ar": 18,   "K": 19,  "Ca": 20,\
                 "Sc": 21,  "Ti": 22,   "V": 23,  "Cr": 24,  "Mn": 25,\
                 "Fe": 26,  "Co": 27,  "Ni": 28,  "Cu": 29,  "Zn": 30,\
                 "Ga": 31,  "Ge": 32,  "As": 33,  "Se": 34,  "Br": 35,\
                 "Kr": 36,   "I": 53}

    if at in at_num_dict:
        atno = at_num_dict[at]
        return(atno)
    else:
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: Unsupported element: " + str(at)+ " \n")
def sym2mass(at):
    """
    Returns atomic mass

            Input:
                    at (str):  Atomic symbol

            Returns:
                    sym2mass (float): Atomic mass of the most abundant isotope
             
            Ref: https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
    """
    mass_dict = {"H" :  1.00782503223,  "He" :  4.00260325413, "Li" :  7.0160034366 ,  "Be" :  9.012183065  ,  "B" : 11.00930536   ,\
                 "C" : 12.0000000    ,   "N" : 14.00307400443,  "O" : 15.99491461957,  "F"  : 18.99840316273, "Ne" : 19.9924401762 ,\
                "Na" : 22.9897692820 ,  "Mg" : 24.985041697  , "Al" : 26.98153853   , "Si"  : 27.97692653465,  "P" : 30.97376199842,\
                 "S" : 31.9720711744 ,  "Cl" : 34.968852682  , "Ar" : 39.9623831237 ,  "K"  : 38.9637064864 , "Ca" : 39.962590863  ,\
                "Ga" : 68.9255735    ,  "Ge" : 73.921177761  , "As" : 74.92159457   , "Se"  : 75.9165218    , "Br" : 78.9183376    ,\
                "Kr" : 83.9114977282 ,   "I" :126.9044719}

    if at in mass_dict:
        sym2mass = mass_dict[at]
        return(sym2mass)
    else:
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: Unsupported element: " + str(at)+ " \n")
def NFC(at):
    """
    Returns the number of frozen core electrons

            Input:
                    at (str): Atomic symbol

            Returns:
                    NFC (int): Number of frozen core electrons
    """

    if at in ['H','He']:
        return(0)
    elif at in ['Li','Be','B','C','N','O','F','Ne','Na','Mg']:
        return(2)
    elif at in ['Al','Si','P','S','Cl','Ar','K','Ca','Fe']:
        return(10)
    elif at in ['Ga','Ge','As','Se','Br','Kr']:
        return(18)
    elif at in ['I']:
        return(36)
    else:
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: Unsupported element: " + str(at)+ " \n")
def nanb(at):
    """
    Returns the number of alpha and beta electron in the element.

            Parameters:
                    at (char): Symbol of the element

            Returns:
                    nanb (list, int): Number of alpha electrons, number of beta electrons in the atom 
    """

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

