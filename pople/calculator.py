import os, sys, string
import linecache, math
import numpy as np
import datetime , time

from pople import orca_g4mp2
from pople import atno
from pople import NFC 
from pople import nanb 
from pople import orca_g4mp2_ctrl

def calculator(**kwargs):

    frozengeom = 'false'
    nproc = '1' 
    mem_mb = '1000'

    if 'code' in kwargs:
        code=kwargs['code']
    if 'code_exe' in kwargs:
        code_exe=kwargs['code_exe']
    if 'method' in kwargs:
        method=kwargs['method']
    if 'xyz' in kwargs:
        geom=kwargs['xyz']
    if 'nproc' in kwargs:
        nproc=kwargs['nproc']
    if 'mem_mb' in kwargs:
        mem_mb=kwargs['mem_mb']
    if 'frozengeom' in kwargs:
        frozengeom=kwargs['frozengeom']

    if frozengeom == 'true':
        if 'freqcmi' in kwargs:
            freq=kwargs['freqcmi']
    
    #=== start timer
    start_time_main = time.time()

    #=== update val dictionary
    val={}

    val = orca_g4mp2_ctrl(method)
    #val.update(new_val1)

    val["maxcore_mb"] = str(mem_mb)
    val["nproc"] = str(nproc)

    if code == 'orca':
        val["orca_exe"]=code_exe

    #=== process geometry block
    geom_line=geom.strip().split()

    charge=int(geom_line[0])
    geom_line.pop(0)

    multip=int(geom_line[0])
    geom_line.pop(0)

    Nat = int((len(geom_line))/4)
 
    sym=[]
    for iat in range(0,len(geom_line),4):
        sym.append(geom_line[iat])

    with open("inp.xyz", "w") as inp_x:
        inp_x.write(str(Nat) + " \n")
        inp_x.write(str(charge) +" "+ str(multip) + " \n")
        for iat in range(0,len(geom_line),4):
            inp_x.write(geom_line[iat]+' '+geom_line[iat+1]+' '+geom_line[iat+2]+' '+geom_line[iat+3])
            inp_x.write("\n")

    #=== process frequencies block
    if frozengeom == 'true':
        val["FROZEN_GEOM"] = "true"
        freq_line=freq.strip().split()
        with open("freq.txt", "w") as inp_x:
            for ifreq in range(len(freq_line)):
                inp_x.write(str(freq_line[ifreq]) + " \n")

    val["Ntotal"] = 0
    val["Ntotale"] = 0
    val["Ntotalecore"] = 0

    for iat in range(Nat):
        na_nb_l = nanb(sym[iat])
        na = na_nb_l[0]
        nb = na_nb_l[1]
        val["Ntotal"] = val["Ntotal"] + na + nb
        val["Ntotale"] = val["Ntotale"] + atno(sym[iat])
        val["Ntotalecore"] = val["Ntotalecore"] + NFC(sym[iat])

    val["Ntotal"] = val["Ntotal"] - charge
    val["Ntotale"] = val["Ntotale"] - charge
    val["Ntotalecore"] = val["Ntotalecore"]

    if Nat == 1: 
        val["isatom"] = "true"
    else:
        val["isatom"] = "false"

    if code == 'orca':
        if method == 'g4mp2' or method == 'g4mp2-xp':
            orca_g4mp2(val, start_time_main)

    U0 = val["U0"]
    UT = val["UT"]
    HT = val["HT"]
    breakdown = val["breakdown"]

    os.system('rm -f inp.xyz freq.txt')
    return(U0, UT, HT, breakdown)
