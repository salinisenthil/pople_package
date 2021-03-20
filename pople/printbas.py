import os, sys, string
import linecache, math
import numpy as np


def printbas(fname, at, install_dir_val):  
    """
    Prints the basis set parameters into the input file input.com

            Parameters:
                    fname (char): Basis set name
                    at (char): Symbol of the element
                    install_dir_val (char) : Path to the directory containing basis set files
    """

    basisSet_fpath = fname
    start_phrase = "NewGTO   "+ at
    print(start_phrase)
    num_lines_bas = sum(1 for line_tmp1 in open(basisSet_fpath,"r"))

    for temp_num, temp_l in enumerate(open(basisSet_fpath,"r")):
        if start_phrase in temp_l.strip():
            if start_phrase == temp_l.strip():
                 bas_start_lno = temp_num+1
                 break

    #print(bas_start_lno)
    with open("input.com", "a") as new_f:

        for l1 in range(bas_start_lno,num_lines_bas):
            #print(l1)
            req_line_1 = linecache.getline(basisSet_fpath, l1)
            if "end" in req_line_1.strip():
                break
            else:
                new_f.write(req_line_1)

        new_f.write("   end\n")
# fname = value of the  basis set name from the inp file,  install_dir_val = value from initial input file by user
