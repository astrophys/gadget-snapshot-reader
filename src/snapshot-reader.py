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
import numpy as np
import grp
import getpass
import matplotlib
import matplotlib.pyplot as plt

from error import exit_with_error
from functions import read_gadget2_snapshot



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
    (pL, header) = read_gadget2_snapshot(Path=path, Endian=edn)


    ###### Let's try to plot the mass density ######
    array = np.zeros([100,100],dtype=np.float64)
    # gas, dm and stars 
    for k in range(len(pL)):
        # Check logic here
        i = int(pL[k].posV[0]/header.boxSize * array.shape[0])
        j = int(pL[k].posV[1]/header.boxSize * array.shape[1])
        #print(i,j)
        array[i,j] += pL[k].mass
    
    array = array / (header.npartV[0] + header.npartV[1] + header.npartV[2] +
                     header.npartV[3] + header.npartV[4] + header.npartV[5])
    array = np.log(array)
    fig, ax = plt.subplots()
    im = ax.imshow(array)
    plt.show()
    
    ###### Let's try to plot the Carbon mass Frac######
    array = np.zeros([100,100],dtype=np.float64)
    # gas, dm and stars 
    for k in range(len(pL)):
        if(pL[k].type == 0 or pL[k].type ==4):
            # Check logic here
            i = int(pL[k].posV[0]/header.boxSize * array.shape[0])
            j = int(pL[k].posV[1]/header.boxSize * array.shape[1])
            #print(i,j)
            array[i,j] += pL[k].Cf
    
    array = array / (header.npartV[0] + header.npartV[1] + header.npartV[2] +
                     header.npartV[3] + header.npartV[4] + header.npartV[5])
    array = np.log(array)
    fig, ax = plt.subplots()
    im = ax.imshow(array)
    plt.show()
    
    
    print("\n\nEnded : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))
    sys.exit(0)

if __name__ == "__main__":
    main()
