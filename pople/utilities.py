import os
import datetime , time

def start_time():
    day = datetime.datetime.today().strftime("%A")
    date = datetime.date.today()
    cur_date = date.strftime("%B %d, %Y")
    now = datetime.datetime.now()
    cur_time = now.strftime(" %H:%M:%S")

    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *  Start time:  " + day + "  "+ cur_date + "  " + cur_time + "  " + time.tzname[0] +"  *\n\n")

def final_time():
    day = datetime.datetime.today().strftime("%A")
    date = datetime.date.today()
    cur_date = date.strftime("%B %d, %Y")
    now = datetime.datetime.now()
    cur_time = now.strftime(" %H:%M:%S")
   
    with open("Thermochemistry.out", "a") as ther_chem:
        ther_chem.write(" *  Final time:  " + day + "  "+ cur_date + "  " + cur_time + "  " + time.tzname[0] +"  *\n\n")
   
    #os.system("rm -f g* sc* fr* inp* t* l* N* e*")
    #exit()

def print_header():
    with open("Thermochemistry.out", "a") as ther_chem:
        cs1 = ( f""" 
 +--------------------------------------------------------------------+
 |                                                                    |
 |   Pople, A toolkit for ab initio thermochemistry (v21.4.1)         |
 |   https://moldis-group.github.io/pople/                            |
 |                                                                    |
 |   Standard enthalpy parameters available for the following atoms   |
 |                                                                    |
 |   ----                                                             |
 |   |H |                                                             |
 |   --------    ---------------------                                |
 |   |Li| Be|    |B  |C  |N  |O  |F  |                                |
 |   |Na| Mg|    |Al |Si |P  |S  |Cl |                                |
 |   |K | Ca|    |Ga |Ge |As |Se |Br |                                |
 |   --------    ---------------------                                |
 |                                                                    |
 |   Credits: see https://moldis-group.github.io/pople/credits.html   |
 |                                                                    |
 +--------------------------------------------------------------------+
             
""")
        ther_chem.write(cs1)

#*  Constants and conversion factors  *

#Planck's constant, h     = 6.6260700400000E-34 J s
#Boltzman's constant, kB  = 1.3806485200000E-23 J K^-1
#Rydberg's constant, Ry   = 1.0973731568508E+07 m^-1
#Elementary charge, e     = 1.6021764870000E-19 C
#Speed of light, c        = 2.9979245800000E+08 m s^-1
#Avogadro number          = 6.0221417900000E+23 mol^-1

#1 hartree                = 2.7211388291762E+01 eV
#1 hartree                = 6.2750957099203E+02 kcal/mol
#1 hartree                = 2.6255000450306E+03 kJ/mol
