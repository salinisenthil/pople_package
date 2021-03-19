def uniqatoms(sym):  # to fine uniq atoms and the number of occurance of uniq atoms in the molecule
    """
    Returns the number of unique atoms in the molecule, along with their symbols and number of occurances

            Parameters:
                    sym (list, char): List of atoms in the molecule

            Returns:
                    uniqat_d (dict): Details of number of unique atoms in the molecule
                    uniqat_d["N_ua"] (int) : Number of unique atom types in the molecule
                    uniqat_d["uniq_sym"] (list, char) : List of symbols of the unique atoms present in the molecule
                    uniqat_d["uan"] (list, int) : List of number of occurances of the unique atoms present
    """

# we want 'ua' unique atom types
# we want 'uan' no. of unique atoms of each type
# Ex. C20H42
# N_ua = 2
# ua = 'C' 'H
# uan = 20  42
    uan = []
    uniq_sym = list(set(sym))  # ua
    N_ua = len(uniq_sym)
    for tmp_i in range(len(uniq_sym)):
        count_a = 0
        for tmp_j in range(len(sym)):
            if uniq_sym[tmp_i] == sym[tmp_j]:
                count_a = count_a + 1
        uan.append(count_a)
        uniqat_d = {}
        uniqat_d["N_ua"] = int(N_ua)
        uniqat_d["uniq_sym"] = uniq_sym
        uniqat_d["uan"] = uan
        #print(uniqat_d)
    #return(N_ua, uniq_sym, uan)  ### revisit, how do you want variables returned?
    return(uniqat_d)
