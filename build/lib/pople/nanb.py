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

