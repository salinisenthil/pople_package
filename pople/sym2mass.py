def sym2mass(at):
    """
    Returns the atomic mass of the given atom

            Parameters:
                    at (char): Symbol of the element

            Returns:
                    sym2mass (int): Atomic mass of the atom
    """

    mass_dict_1 ={"H" : 1.00783,"He" : 4.00260,"Li" : 7.01600,"Be" : 9.01218,"B" : 11.00931,"C" : 12.0000000,"N" : 14.00307,"O" : 15.99491,"F" : 18.99840,"Ne" : 19.99244,"Na" : 22.98977, \
        "Mg" : 23.98505,"Al" : 26.98154,"Si" : 27.97693,"P" : 30.97376,"S" : 31.97207,"Cl" : 34.96885,"Ar" : 39.96238,"K" : 38.96371,"Ca" : 39.96259,"Ga" : 68.92558,"Ge" : 73.92118, \
        "As" : 74.92160,"Se" : 79.91652,"Br" : 78.91834, "Kr" : 83.91151, "Fe" : 55.845, "I" : 126.90447}
    if at in mass_dict_1:
        sym2mass = mass_dict_1[at]
        return(sym2mass)
    else:   # this is equivalent of default in switch case, the write statements needs editing , revisit
        with open("Thermochemistry.out", "a") as ther_chem:
            ther_chem.write("Error: unknown element type encoutered in sym2mass in module_geom.f90:  " + str(at)+ "\n")
