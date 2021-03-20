def principal_coord():
    """
    Returns eigen value
    Requires input.xyz to be present in pwd. input.xyz contains the cartesian coordinates of the molecule/radical/ion in xyz format
    Used to characterize if the system is linear/non-linear

            Returns:
                Ievals  (list, float): Eigen values
    """

    with open("input.xyz","r") as xyz_f:
        num_l_xyz = sum(1 for l in xyz_f)
        Nat = int(linecache.getline("input.xyz",1).strip())
        l1_chr_mul = linecache.getline("input.xyz",2)
        sym = []   # will contain list of atom symbols in the mol # same order
        Rlist = []
        for tmp_l in range(3,num_l_xyz+1):   ### what does NH do ???
            l1_xyz_1 = linecache.getline("input.xyz",tmp_l)
            l1_lsp1 = l1_xyz_1.split()
            sym.append(l1_lsp1[0])
            Rlist.append(float(l1_lsp1[1]))
            Rlist.append(float(l1_lsp1[2]))
            Rlist.append(float(l1_lsp1[3]))
        conv = np.array(Rlist)
        R_coord = np.resize(conv,[Nat,3])

        rCM_pre = []
        mass_list = []
        for tmp_j in range(Nat):
            #print(sym[tmp_j])
            mass_atom = sym2mass(sym[tmp_j])
            mass_list.append(mass_atom)
            new1 = R_coord[tmp_j] * mass_atom
            rCM_pre.append(new1)
            #rCM_pre = R_coord[tmp_j][0] + rCM_pre[0])
        rCM_sum = np.sum(rCM_pre,axis=0)
        tot_mass = sum(mass_list)
        rCM = rCM_sum/tot_mass

        new_coord1 = []
        r_skew_symm = np.zeros((3,3))
        momin = 0.0
        for tmp_k in range(Nat):
            sub_cm = R_coord[tmp_k] - rCM
            new_coord1.append(sub_cm)
            r_skew_symm[0][0] = 0.0
            r_skew_symm[0][1] = -sub_cm[2]
            r_skew_symm[0][2] = sub_cm[1]

            r_skew_symm[1][0] = sub_cm[2]
            r_skew_symm[1][1] = 0.0
            r_skew_symm[1][2] =-sub_cm[0]

            r_skew_symm[2][0] =-sub_cm[1]
            r_skew_symm[2][1] = sub_cm[0]
            r_skew_symm[2][2] = 0.0

            momin = momin + (sym2mass(sym[tmp_k]) * np.matmul(np.transpose(r_skew_symm), r_skew_symm) )

        eig_val, eig_vec = np.linalg.eig(momin)

        idx = eig_val.argsort()[::+1]
        eig_val = eig_val[idx]
        eig_vec = eig_vec[:,idx]

        Ievals = eig_val


        print("MOM== ", Ievals)
        return(Ievals)
