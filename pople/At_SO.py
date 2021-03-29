def At_SO(at, charge):  # element symbol, charge
    """
    Returns energy correction for spin orbit coupling of the given element.

            Parameters:
                    at (char): Symbol of the element
                    charge (int): Charge on the element

            Returns:
                    At_SO (float): Energy correction for spin orbit coupling of the element
    """

    At_SO = 0
    if charge == 0:
        S_dict_0 = {"B" : -0.05, "C" : -0.14, "O" : -0.36, "F" : -0.61, "Al" : -0.34, "Si" : -0.68, "S" : -0.89, "Cl" : -1.34, "Ga" : -2.51, "Ge" : -4.41, "Se" : -4.3,  \
            "Br" : -5.6, "Fe" : -1.84, "I" : -11.548}
        if at in S_dict_0:  # this "if" executes only of the incoming atom "at" is present in the dict, else it returns 0 (this is equivalent to default in the f90 script)
            At_SO = S_dict_0[at]
            return(At_SO)
        else:
            return(At_SO)

    elif charge == +1:
        S_dict_1 = {"C"  :-0.2, "N"  :-0.43, "F"  :-0.67, "Ne" :-1.19, "Si" :-0.93, "P"  :-1.43, "Cl" :-1.68, "Ar" :-2.18, "Ge" :-5.37, "As" :-8.04, "Br" :-6.71, "Kr" :-8.16, "I"  :-14.028}
        if at in S_dict_1:  # this "if" executes only of the incoming atom "at" is present in the dict, else it returns 0 (this is equivalent to default in the f90 script)
            At_SO = S_dict_1[at]
            return(At_SO)
        else:
            return(At_SO)

    elif charge == -1:
        S_dict_n1 = {"B": -0.03, "O": -0.26, "Al": -0.28, "P": -0.45, "S": -0.88}
        if at in S_dict_n1:  # this "if" executes only of the incoming atom "at" is present in the dict, else it returns 0 (this is equivalent to default in the f90 script)
            At_SO = S_dict_n1[at]
            return(At_SO)
        else:
            return(At_SO)

    else:
        return(At_SO)

