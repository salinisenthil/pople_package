import numpy as np
from pople import cmi2au, j2au, kB, c, h
from pople import nanb

def thermal(isatom, freq, scalfac,linnonlin,T):
    """
    Returns atomic energy correction due to spin orbit coupling

            Input:
                    isatom (str)    : "true" for atom, "false" for molecule
                    freq (list)     : unscaled harmonic frequencies
                    scalfac (float) : scaling factor
                    linnonlin (str) : "L" for linear, "NL" for nonlinear
                    T(float)        : Temperature

            Returns:
                    At_SO (float): Energy correction due to spin orbit coupling
    """
    if isatom != "true":
        nfreq = len(freq)

        vib_temp = []
        for ifreq in range(nfreq):
            freq[ifreq] = float(freq[ifreq]) * float(scalfac)
            vib_temp_new = c * 100.0 *  h * float(freq[ifreq]) / kB
            vib_temp.append(vib_temp_new)

        dE_vib = 0
        for ifreq in range(nfreq):
            dE_vib = dE_vib + kB * vib_temp[ifreq] * j2au * ( 0.5 + 1 / ( np.exp(vib_temp[ifreq]/T) - 1) )

        dE_ZPE = 0.5 * sum(freq) * cmi2au

        if linnonlin == "L":
            dE_rot = kB * T * j2au
        elif linnonlin == "NL":
            dE_rot = kB * T * j2au * (3.0/2.0)
        else:
            with open("Thermochemistry.out", "a") as ther_chem:
                ther_chem.write("ERROR: unknown entry for linear/nonlinear")
    else:
        dE_ZPE  = 0
        dE_vib  = 0
        dE_rot  = 0

    dE_tra = kB * T * j2au * (3.0/2.0)
    dE_thermal = (dE_vib - dE_ZPE) + dE_rot + dE_tra

    return(dE_ZPE, dE_vib, dE_rot, dE_tra, dE_thermal)

def HLC_g4mp2(charge, multip, sym, isatom, HLC_params):

    Nat = len(sym)
    Ntotal = 0

    for tmp_k in range(Nat):
        na_nb_l = nanb(sym[tmp_k])
        na = na_nb_l[0]
        nb = na_nb_l[1]
        Ntotal = Ntotal + na + nb

    Ntotal = Ntotal - charge
    Na = (multip + Ntotal - 1 ) / 2.0
    Nb = Na + 1 - multip

    AA = HLC_params[0]
    BB = HLC_params[1]
    CC = HLC_params[2]
    DD = HLC_params[3]
    ApAp = HLC_params[4]
    EE = HLC_params[5]

    if isatom != "true":
        if Na == Nb:
            HLC = - AA * Nb
        else:
            HLC = - ApAp * Nb - BB * (Na- Nb)  
    else:
         HLC = - (CC * Nb ) - (DD * (Na - Nb) )

    if isatom == "true":
        if sym[0] == "Be" and charge ==  0 : HLC = -EE
        if sym[0] == "Mg" and charge ==  0 : HLC = -EE
        if sym[0] == "Ca" and charge ==  0 : HLC = -EE
        if sym[0] == "Li" and charge == -1 : HLC = -EE
        if sym[0] == "Na" and charge == -1 : HLC = -EE
        if sym[0] == "K"  and charge == -1 : HLC = -EE
    elif Nat == 2:
        if charge == 0:
            if sym[0] == "Li" and sym[1] == "Li": HLC = -EE
            if sym[0] == "Na" and sym[1] == "Na": HLC = -EE
            if sym[0] == "K"  and sym[1] ==  "K": HLC = -EE
            if sym[0] == "Li" and sym[1] == "Na": HLC = -EE
            if sym[0] == "Na" and sym[1] == "Li": HLC = -EE
            # if ( (trim(sym(1)) .eq. 'Be') .and. (trim(sym(2)) .eq. 'H' ) ) HLC = -EE
            # if ( (trim(sym(1)) .eq. 'H' ) .and. (trim(sym(2)) .eq. 'Be') ) HLC = -EE

    HLC  = HLC / 1000.0

    return(HLC)

