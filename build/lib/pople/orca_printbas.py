import os, sys, string
import linecache, math
import numpy as np


try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources


from . import basis



def orca_printbas(fname, at):  
    """
    Prints the basis set parameters into the input file input.com

            Parameters:
                    fname (char): Basis set name
                    at (char): Symbol of the element
    """

    #bfile = pkg_resources.open_text(templates, 'GTBAS1')
    #bfile_r = pkg_resources.read_text(templates, 'GTBAS1')
    bfile_r = pkg_resources.read_text(basis, fname)
    with open("temp_bas","w") as nbfile:
        nbfile.write(bfile_r)
    #basisSet_fpath = fname
    start_phrase = "NewGTO   "+ at
    #print(start_phrase)
    num_lines_bas = sum(1 for line_tmp1 in open("temp_bas","r"))

    for temp_num, temp_l in enumerate(open("temp_bas","r")):
        if start_phrase in temp_l.strip():
            if start_phrase == temp_l.strip():
                 bas_start_lno = temp_num+1
                 break

    with open("input.com", "a") as new_f:
        linecache.clearcache()
        for l1 in range(bas_start_lno,num_lines_bas):
            req_line_1 = linecache.getline("temp_bas", l1)
            if "end" in req_line_1.strip():
                break
            else:
                new_f.write(req_line_1)

        new_f.write("   end\n")

    os.system("rm -f temp_bas")
# fname = value of the  basis set name from the inp file,  install_dir_val = value from initial input file by user
