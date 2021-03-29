import linecache, os, sys

def get_geom(setting):
    """
    Creates geometry files for pople run.
            Parameters:
                    setting (char): Type of calculation requested.

    """
    current_dir_path = os.getcwd()
    input_file = current_dir_path + "/pople.inp"
# find line no of start of geom.xyz and second geom file(is enabled)
    with open(input_file, "r") as inp_f:
        for num, line in enumerate(inp_f,1):
            if "geom_1" in line:
                geom_lno = num
            if setting == "IP" or setting == "EA" or setting == "PA" or setting == "BE":
                if "geom_2" in line:
                    sec_geom_lno = num
# creates geom.xyz
        with open("geom.xyz", "w") as g_xyz1:
            Nat = int(linecache.getline("pople.inp", geom_lno+1))
            char_mul = linecache.getline("pople.inp", geom_lno+2)
            g_xyz1.write(str(Nat) + "\n")
            g_xyz1.write(char_mul)
            for tmp_i in range(geom_lno+3, geom_lno+3+Nat):
                nl = linecache.getline("pople.inp", tmp_i)
                g_xyz1.write(nl)

#creates second geom.xyz
    if setting == "IP" or setting == "EA" or setting == "PA" or setting == "BE":
        if setting == "IP": fname2 = "geom_cation_IP.xyz"
        if setting == "EA": fname2 = "geom_anion_EA.xyz"
        if setting == "PA": fname2 = "geom_cation_PA.xyz"
        if setting == "BE": fname2 = "geom_monomer_A.xyz"
        with open(fname2, "w") as g_xyz1:
            Nat = int(linecache.getline("pople.inp", sec_geom_lno+1))
            char_mul = linecache.getline("pople.inp", sec_geom_lno+2)
            g_xyz1.write(str(Nat) + "\n")
            g_xyz1.write(char_mul)
            for tmp_i in range(sec_geom_lno+3, sec_geom_lno+3+Nat):
                nl = linecache.getline("pople.inp", tmp_i)
                g_xyz1.write(nl)
