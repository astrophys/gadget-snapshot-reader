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

    fin = open(path, "rb")
    buf    = fin.read(4)        # buffer
    headB  = fin.read(256)
    buf    = fin.read(4)        # buffer

    ##
    header = HEADER_V2(Bytes=headB, Endian="little")

    startTime = time.time()

    
    print("\n\nEnded : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))
    sys.exit(0)

if __name__ == "__main__":
    main()
