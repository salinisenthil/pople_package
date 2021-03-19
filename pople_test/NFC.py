def NFC(at):   # returns no. of forzen core electrons  # Why not use elif??? revisit
    """
    Returns the number of frozen core electrons in the atom

            Parameters:
                    at (char): Symbol of the element

            Returns:
                    NHC (int): Number of frozen core electrons in the atom
    """

    if at in ['H','He']:
        return(0)
    if at in ['Li','Be','B','C','N','O','F','Ne','Na','Mg']:
        return(2)
    if at in ['Al','Si','P','S','Cl','Ar','K','Ca','Fe']:
        return(10)
    if at in ['Ga','Ge','As','Se','Br','Kr']:
        return(18)
    if at in ['I']:
        return(36)
