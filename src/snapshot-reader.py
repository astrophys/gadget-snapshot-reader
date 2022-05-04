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

    ########### Get mass ###########
    nPartWSpecMass = 0              # N_Local_Part_w_Mass in gadget2csv
    for i in range(len(header.massV)):
        if(header.massV[i] == 0):
            nPartWSpecMass += header.npartV[i]

    ########### Particle Positions ###########
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


    
    print("\n\nEnded : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))
    sys.exit(0)

if __name__ == "__main__":
    main()
