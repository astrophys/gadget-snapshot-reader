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
    fin = open(path, "rb")
    pL = read_gadget2_snapshot(File=fin, Endian=edn)
    ###### Ensure that all the data is read and that the EOF has been reached ######
    #   https://stackoverflow.com/a/53792531/4021436
    curPos = fin.tell()
    eof    = fin.seek(0,2)
    lastByte = fin.read(1)
    if(curPos - eof != 0 or len(lastByte) != 0):
        exit_with_error("ERROR!!! EOF ({}) not reached! Current Position = {}, lastByte = {}".
                        format(eof,curPos, lastByte))
    fin.close()
    
    
    print("\n\nEnded : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))
    sys.exit(0)

if __name__ == "__main__":
    main()
