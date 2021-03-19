def uniqatoms(sym):  # to fine uniq atoms and the number of occurance of uniq atoms in the molecule
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
