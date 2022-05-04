# AUTHOR : Ali Snedden
# LICENSE: MIT
# DATE   : 04/06/2022
# PURPOSE:
#   
#   
# NOTES:
#
# TODO:
#
# DEBUG:
#
import os
import sys
import time
import subprocess
from error import exit_with_error
import grp
import getpass
from classes import HEADER_V2
from classes import PARTICLE
from classes import read_metal
import struct


def print_help(Arg):
    """
    ARGS:
        Arg     : (int) exit value
    DESCRIPTION:
        Print Help. Exit with value arg
    RETURN:
        N/A
    DEBUG:
        1. Tested, it worked
    FUTURE:
    """
    sys.stdout.write(
            "\nUSAGE : python3 snapshot-reader.py snapshot_XXX\n\n"
            "       snapshot_XXX : (str) path to a Gadget snapshot file\n"
            "                         \n")
    sys.exit(Arg)



def read_gadget2_snapshot(Path=None):
    """
    ARGS:
        Path    : (str) path to gadget2 snapshot
    DESCRIPTION:
        This function reads a snapshot file that was created by gadget2-ali, which 
        includes feedback (i.e. chemical enrichment) and star formation.
    RETURN:
    DEBUG:
    FUTURE:
    """
    print("Hello World")






def main():
    """
    ARGS:
    DESCRIPTION:
    RETURN:
    DEBUG:
    FUTURE:
    """
    ### Python 3 required
    if(sys.version_info[0] != 3):
        exit_with_error("ERROR!!! Runs with python3, NOT python-{}\n\n".format(
                    sys.version_info[0]))

    ######### Get Command Line Options ##########
    if(len(sys.argv) != 2):
        print_help(1)
    elif(len(sys.argv) == 2 and (sys.argv[1] == "--help" or sys.argv[1] == "-help")):
        print_help(0)
    else:
        path=sys.argv[1]
    # Endianness matters
    endian = "little"
    if(endian == "little"):
        edn = "<"
    elif(Endian.lower() == "big"):
        edn = ">"
    else:
        exit_with_error("ERROR!!! Incorrect option ({}) for Endianess".format(Endian))

    startTime = time.time()
    floatSize = 4
    intSize   = 4
    fin = open(path, "rb")
    buf    = fin.read(4)        # buffer
    headB  = fin.read(256)
    buf    = fin.read(4)        # buffer


    ########### Header ###########
    header = HEADER_V2(Bytes=headB, Endian="little")
    ngas   = header.npartV[0]
    ndm    = header.npartV[1]
    nstar  = header.npartV[4]
    npart  = (header.npartV[0] + header.npartV[1] + header.npartV[2] + header.npartV[3] +
              header.npartV[4] + header.npartV[5])
    gasMass= header.massV[0]


    ########## Allocate Particles ########## 
    # Ugly to do in python, oh well
    pL     = []                 # List of particles
    for i in range(npart):
        pL.append(PARTICLE())


    ########### Mass ###########
    nPartWSpecMass = 0              # N_Local_Part_w_Mass in gadget2csv
    for i in range(len(header.massV)):
        if(header.massV[i] == 0):
            nPartWSpecMass += header.npartV[i]


    ########### Positions ###########
    buf    = fin.read(4)        # buffer
    i      = 0                  # Particle Index
    for n in range(6):
        for idx in range(header.npartV[n]):
            bytes = fin.read(floatSize)
            x = struct.unpack(edn+"f", bytes)[0]
            bytes = fin.read(floatSize)
            y = struct.unpack(edn+"f", bytes)[0]
            bytes = fin.read(floatSize)
            z = struct.unpack(edn+"f", bytes)[0]
            pL[i].posV[0] = x
            pL[i].posV[1] = y
            pL[i].posV[2] = z
            i += 1
    buf    = fin.read(4)        # buffer


    ########### Velocities ###########
    buf    = fin.read(4)        # buffer
    i      = 0                  # Particle Index
    for n in range(6):
        for idx in range(header.npartV[n]):
            bytes = fin.read(floatSize)
            vx = struct.unpack(edn+"f", bytes)[0]
            bytes = fin.read(floatSize)
            vy = struct.unpack(edn+"f", bytes)[0]
            bytes = fin.read(floatSize)
            vz = struct.unpack(edn+"f", bytes)[0]
            pL[i].vV[0] = x
            pL[i].vV[1] = y
            pL[i].vV[2] = z
            i += 1
    buf    = fin.read(4)        # buffer


    ########### Particle IDs ###########
    buf    = fin.read(4)        # buffer
    i      = 0                  # Particle Index
    for n in range(6):
        for idx in range(header.npartV[n]):
            bytes = fin.read(intSize)
            ident = struct.unpack(edn+"i", bytes)[0]
            pL[i].id   = ident
            pL[i].type = n
            i += 1
    buf    = fin.read(4)        # buffer


    ########### Mass ###########
    if(nPartWSpecMass > 0) :
        buf    = fin.read(4)        # buffer
    i = 0
    for n in range(6):
        for idx in range(header.npartV[n]):
            if(header.massV[n] == 0):       # Float comparison
                bytes = fin.read(floatSize)
                mass = struct.unpack(edn+"f", bytes)[0]
                pL[i].mass = mass
            else:
                pL[i].mass = header.massV[n]
            i += 1
    if(nPartWSpecMass > 0) :
        buf    = fin.read(4)        # buffer


    if(ngas > 0):
        ########### Internal Energy ###########
        buf    = fin.read(4)        # buffer
        for i in range(header.npartV[0]):
            bytes = fin.read(floatSize)
            U     = struct.unpack(edn+"f", bytes)[0]
            pL[i].u = U
        buf    = fin.read(4)        # buffer


        ########### Gas Density ###########
        buf    = fin.read(4)        # buffer
        for i in range(header.npartV[0]):
            bytes = fin.read(floatSize)
            rho   = struct.unpack(edn+"f", bytes)[0]
            pL[i].rho = rho
        buf    = fin.read(4)        # buffer
            
    
        ########### Hsml - what is it? ###########
        buf    = fin.read(4)        # buffer
        for i in range(header.npartV[0]):
            bytes = fin.read(floatSize)
            hsml = struct.unpack(edn+"f", bytes)[0]
            pL[i].hsml = hsml
        buf    = fin.read(4)        # buffer
            
    
        ########### SFR ###########
        if(header.flag_sfr == True):
            buf    = fin.read(4)        # buffer
            for i in range(header.npartV[0]):
                bytes = fin.read(floatSize)
                sfr   = struct.unpack(edn+"f", bytes)[0]
                pL[i].sfr = sfr
            buf    = fin.read(4)        # buffer
            

            ########### Stellar age ###########
            if(header.flag_stellarage == True):
                buf    = fin.read(4)        # buffer
                offset = header.npartV[0] + header.npartV[1] + header.npartV[2] + header.npartV[3]
                for i in range(header.npartV[4]):
                    bytes = fin.read(floatSize)
                    age   = struct.unpack(edn+"f", bytes)[0]
                    pL[offset + i].age = age
                buf    = fin.read(4)        # buffer
            
    
    ########### Metals Fraction ###########
    if(header.flag_metals == True):
        if(header.npartV[0] > 0 or header.npartV[4] > 0):

            ########### Carbon ###########
            if(header.flag_Carbon == True):
                read_metal(PL=pL, Header=header, Metal="Carbon", Short="Cf", Fin=fin)

            ########### Nitrogen ###########
            if(header.flag_Nitrogen == True):
                read_metal(PL=pL, Header=header, Metal="Nitrogen", Short="Nf", Fin=fin)
    
            ########### Oxygen ###########
            if(header.flag_Oxygen == True):
                read_metal(PL=pL, Header=header, Metal="Oxygen", Short="Of", Fin=fin)
    
            ########### Florine ###########
            if(header.flag_Florine== True):
                read_metal(PL=pL, Header=header, Metal="Florine", Short="Ff", Fin=fin)
    
            ########### Neon ###########
            if(header.flag_Neon == True):
                read_metal(PL=pL, Header=header, Metal="Neon", Short="Nef", Fin=fin)
    
            ########### Sodium ###########
            if(header.flag_Sodium == True):
                read_metal(PL=pL, Header=header, Metal="Sodium", Short="Naf", Fin=fin)
    
            ########### Magnesium ###########
            if(header.flag_Magnesium == True):
                read_metal(PL=pL, Header=header, Metal="Magnesium", Short="Mgf", Fin=fin)
    
            ########### Aluminum ###########
            if(header.flag_Aluminum == True):
                read_metal(PL=pL, Header=header, Metal="Aluminum", Short="Alf", Fin=fin)
    
            ########### Silicon ###########
            if(header.flag_Silicon == True):
                read_metal(PL=pL, Header=header, Metal="Silicon", Short="Sif", Fin=fin)
    
            ########### Phosphorus ###########
            if(header.flag_Phosphorus == True):
                read_metal(PL=pL, Header=header, Metal="Phosphorus", Short="Pf", Fin=fin)
    
            ########### Sulfur ###########
            if(header.flag_Sulfur == True):
                read_metal(PL=pL, Header=header, Metal="Sulfur", Short="Sf", Fin=fin)
    
            ########### Chlorine ###########
            if(header.flag_Chlorine == True):
                read_metal(PL=pL, Header=header, Metal="Chlorine", Short="Clf", Fin=fin)
    
            ########### Argon ###########
            if(header.flag_Argon == True):
                read_metal(PL=pL, Header=header, Metal="Argon", Short="Arf", Fin=fin)
    
            ########### Potassium ###########
            if(header.flag_Potassium == True):
                read_metal(PL=pL, Header=header, Metal="Potassium", Short="Kf", Fin=fin)
    
            ########### Calcium ###########
            if(header.flag_Calcium == True):
                read_metal(PL=pL, Header=header, Metal="Calcium", Short="Caf", Fin=fin)
    
            ########### Scandium ###########
            if(header.flag_Scandium == True):
                read_metal(PL=pL, Header=header, Metal="Scandium", Short="Scf", Fin=fin)
    
            ########### Titanium ###########
            if(header.flag_Titanium == True):
                read_metal(PL=pL, Header=header, Metal="Titanium", Short="Tif", Fin=fin)
    
            ########### Vanadium ###########
            if(header.flag_Vanadium == True):
                read_metal(PL=pL, Header=header, Metal="Vanadium", Short="Vf", Fin=fin)
    
            ########### Chromium ###########
            if(header.flag_Chromium == True):
                read_metal(PL=pL, Header=header, Metal="Chromium", Short="Crf", Fin=fin)
    
            ########### Manganese ###########
            if(header.flag_Manganese == True):
                read_metal(PL=pL, Header=header, Metal="Manganese", Short="Mnf", Fin=fin)
    
            ########### Iron ###########
            if(header.flag_Iron == True):
                read_metal(PL=pL, Header=header, Metal="Iron", Short="Fef", Fin=fin)
    
            ########### Cobalt ###########
            if(header.flag_Cobalt == True):
                read_metal(PL=pL, Header=header, Metal="Cobalt", Short="Cof", Fin=fin)
    
            ########### Nickel ###########
            if(header.flag_Nickel == True):
                read_metal(PL=pL, Header=header, Metal="Nickel", Short="Nif", Fin=fin)
        else:
            exit_with_error("ERROR!!! 'flag_metals' enabled but no gas or stars present!")
    
    
    
    print("\n\nEnded : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))
    sys.exit(0)

if __name__ == "__main__":
    main()
