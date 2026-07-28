#!/bin/python3
# This Python script generates an initial input file for ECG codes with a reasonable basis building and optimization script for the ECG codes.
#
# Usage examples:
#
# ./ecg_input_file_generator.py -basis_type=RG_0S -particles 3 -masses "1.E30 1.0 1.0" -charges "2.0 -1.0 -1.0" -symmetry "(1+P23)" -eigenvalue 1 -maxbasis 2000 -istart 50 -systemname "He_1Se"
#
# ./ecg_input_file_generator.py -basis_type=RG_0S -particles 4 -masses "1.E30 1.0 1.0 1.0" -charges "3.0 -1.0 -1.0 -1.0" -symmetry "atom-doublet" -eigenvalue 1 -maxbasis 3000 -istart 50 -systemname "Li_2Se"
#
# ./ecg_input_file_generator.py -basis_type=RG_0S -particles 5 -masses "1.E30 1.0 1.0 1.0 1.0" -charges "4.0 -1.0 -1.0 -1.0 -1.0" -symmetry "atom-triplet" -youngoperatorformat "AS" -eigenvalue 1 -maxbasis 5000 -istart 80 -systemname "Be_3Se"
#
# ./ecg_input_file_generator.py -basis_type=RG_2P -particles 7 -masses "25519.04528337 1.0 1.0 1.0 1.0 1.0 1.0" -charges "7.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0" -symmetry "(1+P23)(1+P45)(1-P24)(1-P26-P46)(1-P27-P47-P67)(1-P35)" -eigenvalue 1 -maxbasis 15000 -istart 80 -systemname "N+_3Pe"
#
# Argument basis_type is optional. If present, it will be printed in the header of the input file, where the corresponding descriptor is also optional. 
# If not present, the header will not contain the basis type descriptor.
#
# Note that the permutational symmetry can be supplied either by means of a string with an explicit form of the Young operator
# (e.g. "(1+P23)"), or can be specified using the following keywords:
# "atom-singlet", "atom-doublet", "atom-triplet", "atom-quartet", "atom-quintet", "atom-sextet","atom-septet", "atom-octet", "atom-nonet".
# In such a case the first particle is assumed to be an atomic nucleus, while all others are assumed to be identical electrons.
#
# The value of argument youngoperatorformat is only used if the permutational symmetry is defined through a keyword (e.g. "atom-quartet").
# It can only take two values - "SA" and "AS". The default/natural/textbook order is "SA", meaning that Y = S*A, where S is the symmetrizer and
# A is an antisymmetrizer. The default/natural order always yields n! terms for matrix elements. The reversed order, AS, which is implemented
# here to enable experimenting with things, may yield a somewhat smaller number of terms for high-spin states. At the same time, so far
# there has been no evidense that the reverse order can result in incorrect energies or expectation values.
#
# Argument maxbasis defines the final desired basis size
#
# Argument istart defines the size of the basis when the eigenvalue value solver option switches from 'G' to 'I' (e.g. from LAPACKS's DSYGVX to the Inverse iteration method)
#
# String passed by argument systemname (e.g. "He_1Se") is appended to the names of the input files that are saved with the SAVE_FILE script commands. Note that the eigenvalue number is also appended to the names of the saved files.
#
# In case if masses of common nuclear isotopes are needed, below are their values (in a.u., that is in electron masses). These masses
# were derived using AME 2020 atomic masses and CODATA 2022 values of fundamental constants.
#
# Masses of most common nuclear isotopes:

# M(1H)  =  1836.15267344
#             ±0.00000006
#
# M(1H)  =  1836.152673426(32)  [CODATA 2022]
#
# M(2H)  =  3670.48296764
#             ±0.00000009
#
# M(2H)  =  3670.482967655(63)  [CODATA 2022]
#
# M(3H)  =  5496.92153559
#             ±0.00000024
#
# M(3H)  =  5496.92153551(21)  [CODATA 2022]
#
# M(3He) =  5495.88527989
#             ±0.00000021
#
# M(3He) =  5495.88527984(16)  [CODATA 2022]
#
# M(4He) =  7294.29954170
#             ±0.00000042
#
# M(4He) =  7294.29954171(17)  [CODATA 2022]
#
# M(6Li) = 10961.89865320
#             ±0.00000293
#
# M(7Li) = 12786.39227775
#             ±0.00000752
#
# M(9Be) = 16424.20551744
#             ±0.00014615
#
# M(10B) = 18247.46863262
#             ±0.00002954
#
# M(11B) = 20063.73694390
#             ±0.00002411
#
# M(12C) = 21868.66385131
#             ±0.00000065
#
# M(13C) = 23697.66782831
#             ±0.00000114
#
# M(14C) = 25520.35060830
#             ±0.00000801
#
# M(14N) = 25519.04528337
#             ±0.00000142
#
# M(15N) = 27336.52871215
#             ±0.00000211
#
# M(16O) = 29148.94969878
#             ±0.00000213
#
# M(17O) = 30979.52555500
#             ±0.00000289
#
# M(18O) = 32802.46481984
#             ±0.00000292
#
# Masses of some exotic particles:
#
# M(muon)  = 206.7682827(46)  [CODATA 2022]
#
# M(tauon) = 3477.23(23)  [CODATA 2022]
#
# M(pion±) = 273.13287(47)  [Phys. Lett. B 796, 11 (2019)]

import argparse

parser = argparse.ArgumentParser(description="Generate file")
parser.add_argument("-basis_type", default='NONE', help="Basis type, optional argument. Possible values are RG_0S, RG_1P, RG_2D, RG_2P, CG_0S.", choices=["NONE", "RG_0S", "RG_1P", "RG_2D", "RG_2P", "CG_0S"])
parser.add_argument("-particles", default=3, type=int, choices=[2, 3, 4, 5, 6, 7, 8, 9], help="Number of particles, default is 3")
parser.add_argument("-maxbasis", default=5000, type=int, help="Maximum basis size, default is 5000")
parser.add_argument("-istart", default=50, type=int, help="Size of the basis after which the inverse iteration eigensolver will be used (the G option changes to I)")
parser.add_argument("-systemname", default='SYSTEMNAME', help="System name, default is SYSNAME-01")
parser.add_argument("-symmetry", default='(1)', help="Permutational symmetry, default is (1)")
parser.add_argument("-eigenvalue", default=1 ,type=int, help="Eigenvalue number")
parser.add_argument("-eigval_tolerance", default='1.0E-12', help="Eigensolver tolerance")
parser.add_argument("-masses", default='nonesupplied', help="Masses, default is 1.0 for all")
parser.add_argument("-charges", default='nonesupplied', help="Charges, default is 1.0 for all")
parser.add_argument("-youngoperatorformat", default='SA', choices=['SA', 'AS'], help="Order of symmetrizers and antisymmetryzers in the Young operator")

args = parser.parse_args()

BASIS_TYPE  = args.basis_type
PARTICLES  = args.particles
MAXBASIS = args.maxbasis
ISTART = args.istart
SYSTEMNAME = args.systemname
SYMMETRY = args.symmetry
EIGENVALUE = args.eigenvalue
EIGVAL_TOLERANCE = args.eigval_tolerance
MASSES = args.masses
CHARGES = args.charges
YOUNGOPERATORFORMAT = args.youngoperatorformat

RANDOM_TRIALS = 500  #Number of stochastic trials in BASIS_ENL
OPT_LIMIT_EXPAND = 50*(PARTICLES-1)  #Max number of energy evaluations when optimizing a newly selected basis function in BASIS_ENL
OPT_LIMIT_CYCLE = max(50*(PARTICLES-2), 10) #Max number of energy evaluations when optimizing a function in OPT_CYCLE
LIN_COEFF_THRESHOLD = 3.0 #Maximum allowed magnitude of linear coefficients

#If permutational symmetry was supplied through a keyword
#such as "atom-singlet", "atom-doublet", "atom-triplet", "atom-quartet", etc, then
#we need to generate generate a suitable Young operator

if (PARTICLES == 3):
    if (SYMMETRY == "atom-singlet"):
        S="(1+P23)"
        A=""
    if (SYMMETRY == "atom-triplet"):
        S=""
        A="(1-P23)"

if (PARTICLES == 4):
    if (SYMMETRY == "atom-doublet"):
        S="(1+P23)"
        A="(1-P24)"
    if (SYMMETRY == "atom-quartet"):
        S=""
        A="(1-P23)(1-P24-P34)"

if (PARTICLES == 5):
    if (SYMMETRY == "atom-singlet"):
        S="(1+P23)(1+P45)"
        A="(1-P24)(1-P35)"
    if (SYMMETRY == "atom-triplet"):
        S="(1+P23)"
        A="(1-P24)(1-P25-P45)"
    if (SYMMETRY == "atom-quintet"):
        S=""
        A="(1-P23)(1-P24-P34)(1-P25-P35-P45)"

if (PARTICLES == 6):
    if (SYMMETRY == "atom-doublet"):
        S="(1+P23)(1+P45)"
        A="(1-P24)(1-P26-P46)(1-P35)"
    if (SYMMETRY == "atom-quartet"):
        S="(1+P23)"
        A="(1-P24)(1-P25-P45)(1-P26-P46-P56)"
    if (SYMMETRY == "atom-sextet"):
        S=""
        A="(1-P23)(1-P24-P34)(1-P25-P35-P45)(1-P26-P36-P46-P56)"

if (PARTICLES == 7):
    if (SYMMETRY == "atom-singlet"):
        S="(1+P23)(1+P45)(1+P67)"
        A="(1-P24)(1-P26-P46)(1-P35)(1-P37-P57)"
    if (SYMMETRY == "atom-triplet"):
        S="(1+P23)(1+P45)"
        A="(1-P24)(1-P26-P46)(1-P27-P47-P67)(1-P35)"
    if (SYMMETRY == "atom-quintet"):
        S="(1+P23)"
        A="(1-P24)(1-P25-P45)(1-P26-P46-P56)(1-P27-P47-P57-P67)"
    if (SYMMETRY == "atom-septet"):
        S=""
        A="(1-P23)(1-P24-P34)(1-P25-P35-P45)(1-P26-P36-P46-P56)(1-P27-P37-P47-P57-P67)"

if (PARTICLES == 8):
    if (SYMMETRY == "atom-doublet"):
        S="(1+P23)(1+P45)(1+P67)"
        A="(1-P24)(1-P26-P46)(1-P28-P48-P68)(1-P35)(1-P37-P57)"
    if (SYMMETRY == "atom-quartet"):
        S="(1+P23)(1+P45)"
        A="(1-P24)(1-P26-P46)(1-P27-P47-P67)(1-P28-P48-P68-P78)(1-P35)"
    if (SYMMETRY == "atom-sextet"):
        S="(1+P23)"
        A="(1-P24)(1-P25-P45)(1-P26-P46-P56)(1-P27-P47-P57-P67)(1-P28-P48-P58-P68-P78)"
    if (SYMMETRY == "atom-octet"):
        S=""
        A="(1-P23)(1-P24-P34)(1-P25-P35-P45)(1-P26-P36-P46-P56)(1-P27-P37-P47-P57-P67)(1-P28-P38-P48-P58-P68-P78)"

if (PARTICLES == 9):
    if (SYMMETRY == "atom-singlet"):
        S="(1+P23)(1+P45)(1+P67)(1+P89)"
        A="(1-P24)(1-P26-P46)(1-P28-P48-P68)(1-P35)(1-P37-P57)(1-P39-P59-P79)"
    if (SYMMETRY == "atom-triplet"):
        S="(1+P23)(1+P45)(1+P67)"
        A="(1-P24)(1-P26-P46)(1-P28-P48-P68)(1-P29-P40-P69-P89)(1-P35)(1-P37-P57)"
    if (SYMMETRY == "atom-quintet"):
        S="(1+P23)(1+P45)"
        A="(1-P24)(1-P26-P46)(1-P27-P47-P67)(1-P28-P48-P68-P78)(1-P29-P49-P69-P79-P89)(1-P35)"
    if (SYMMETRY == "atom-septet"):
        S="(1+P23)"
        A="(1-P24)(1-P25-P45)(1-P26-P46-P56)(1-P27-P47-P57-P67)(1-P28-P48-P58-P68-P78)(1-P29-P49-P59-P69-P79-P89)"
    if (SYMMETRY == "atom-nonet"):
        S=""
        A="(1-P23)(1-P24-P34)(1-P25-P35-P45)(1-P26-P36-P46-P56)(1-P27-P37-P47-P57-P67)(1-P28-P38-P48-P58-P68-P78)(1-P29-P39-P49-P59-P69-P79-P89)"

if (SYMMETRY == "atom-singlet" or SYMMETRY == "atom-doublet" or SYMMETRY == "atom-triplet" or SYMMETRY == "atom-quartet" or SYMMETRY == "atom-quintet" or SYMMETRY == "atom-sextet" or SYMMETRY == "atom-septet" or SYMMETRY == "atom-octet" or SYMMETRY == "atom-nonet"):
    if (YOUNGOPERATORFORMAT == "SA"):
        #Default/natural/textbook order
        SYMMETRY = S + A
    if (YOUNGOPERATORFORMAT == "AS"):
        #Unnatural order
        SYMMETRY = A + S

#Print header first
if (BASIS_TYPE != 'NONE'):
    print(" BASIS_TYPE " + BASIS_TYPE)

print(" PARTICLES" + "{:>5}".format(PARTICLES))

print(" MASSES   ", end=" ")
if (MASSES == 'nonesupplied'):
    for _ in range(PARTICLES):
        print("  1.0E0", end=" ")
    print(" ")
else:
    print(MASSES)

print(" CHARGES  ", end=" ")
if (CHARGES == 'nonesupplied'):
    for _ in range(PARTICLES):
        print("  1.0E0", end=" ")
    print(" ")
else:
    print(CHARGES)

print(" SYMMETRY  " + SYMMETRY)
print(" BASIS_SIZE  0")
print(" CURRENT_ENERGY    1.0E30")
print(" WHICH_EIGENVALUE  " + "{:<5}".format(EIGENVALUE))
print(" EIGVAL_TOLERANCE  " + EIGVAL_TOLERANCE)
print(" INVITPARAMETER    1.000001000E0")
print(" LAST_EIGVAL_TOL   1.0E0")
print(" BEST_EIGVAL_TOL   1.0E0")
print(" WORST_EIGVAL_TOL  1.0E-30")
print(" GENERATOR_PARAM   0.7E0  1.0E0  3.0E0")
print(" ==============================")
print("    0  1.0E30      0      0      0")
print(" ==============================")

HFILE="none" #HFILE="hess.dat"
EIGSOLVER='G'

OVERLAP_THRESHOLD=0.95  #Overlap threshold
NUMCYCLES=10 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=999999 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=5 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=10 #STEP in BASIS_ENL
ilist = [3, 5, 10, 15, 20, 30]   # Basis size values for BASIS_ENL and FULL_OPT1
for i in range(len(ilist)):
    kstop = ilist[i]
    if (kstop <= MAXBASIS):    
        if (i==0):
            kstart = 1
        else:
            kstart = ilist[i-1]+1
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))        
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")

OVERLAP_THRESHOLD=0.95  #Overlap threshold
NUMCYCLES=10 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=999999 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=10 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=10 #STEP in BASIS_ENL
for i in range(30, 50, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")

OVERLAP_THRESHOLD=0.97  #Overlap threshold
NUMCYCLES=10 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=3 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=50 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=10 #STEP in BASIS_ENL
for i in range(50, 100, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")

OVERLAP_THRESHOLD=0.98  #Overlap threshold
NUMCYCLES=5 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=3 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=100 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=10 #STEP in BASIS_ENL
for i in range(100, 200, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")

OVERLAP_THRESHOLD=0.98 #Overlap threshold
NUMCYCLES=3 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=3 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=100 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=10 #STEP in BASIS_ENL
for i in range(200, 500, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")

OVERLAP_THRESHOLD=0.98  #Overlap threshold
NUMCYCLES=1 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=3 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=500 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=10 #STEP in BASIS_ENL
for i in range(500, 1000, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")    

OVERLAP_THRESHOLD=0.98  #Overlap threshold
NUMCYCLES=1 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=3 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=500 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=20 #STEP in BASIS_ENL
for i in range(1000, 2000, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")  

OVERLAP_THRESHOLD=0.98  #Overlap threshold
NUMCYCLES=1 #Number of cycles in OPT_CYCLE
OPT_LIMIT_FULLOPT=3 #Maximum allowed number of function evaluations in FULL_OPT1
SAVE_FREQ=5  #How often the file is saved during cycles in OPT_CYCLE
SFR=500 #How often (multiple of the basis size) the next inout_XXXX.txt is created
STEP=50 #STEP in BASIS_ENL
for i in range(2000, 20000, STEP):
    kstart = i+1
    kstop = i+STEP
    if (kstop <= MAXBASIS):
        if (kstart > ISTART):
            EIGSOLVER='I'
        print(" BASIS_ENL  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstart,kstop,1,RANDOM_TRIALS,OPT_LIMIT_EXPAND,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD))
        # print(" FULL_OPT1  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>10}".format(kstop,1,kstop,OPT_LIMIT_FULLOPT,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,60,60,HFILE))        
        print(" OPT_CYCLE  " + EIGSOLVER + "  {:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}".format(kstop,1,kstop,1,1,NUMCYCLES,OPT_LIMIT_CYCLE,OVERLAP_THRESHOLD,LIN_COEFF_THRESHOLD,SAVE_FREQ))
        if (kstop%SFR == 0):
            print(" SAVE_FILE  {:<6}".format(kstop) + " inout_" + SYSTEMNAME + "-" + "{:02d}".format(EIGENVALUE) + "-" + "{:05d}".format(kstop)+".txt")  

print(" ==============================")
