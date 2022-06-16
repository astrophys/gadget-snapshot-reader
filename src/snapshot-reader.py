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
from functions import output_array_as_matlab_int



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
            "\nUSAGE : python3 snapshot-reader.py snapshot_XXX [option] [nGrid] \n\n"
            "       snapshot_XXX : (str) path to a Gadget snapshot file\n"
            "       [option]     : (str) [plot|out_cmp_gas|out_cmp_star|out_gas_matlab] [nGrid]\n"
            "       [nGrid]      : (int) only requred when option == 'out_gas_matlab', is the\n"
            "                      number of grids along one dimension, for a total of (nGrid**3) voxels\n"
            "                         \n")
    sys.exit(Arg)




def main():
    """
    ARGS:
    DESCRIPTION:
    RETURN:
    DEBUG:
        1. Compared output to gadget2/128part-25mpc/snapshot_016/mnc2file_gas.txt
            a) Looked at mnc2file_sph code to figure out how the values were computed in 
               mnc2file_gas.txt
            #) Comparison : 
                > aa=(-header.npartV[4]-header.npartV[1]-4); np.log10(pL[aa].Cf /
                    (pL[aa].mass * 2.45E-04 * 12.011 / 1.336294e+00)) : 
                -2.3081247000154947
                $ grep "2.67651" data/gadget2/128part-25mpc/snapshot_016/mnc2file_gas.txt |
                grep 4.7546 | awk '{print $6}'
                -2.308125
            #) IDENTICAL! I have high confidence here in the values b/c I randomly checked
                          portions of the mnc2file_gas.txt. 
            #) NOTE : I had to use grep above b/c the entries are jumbled due to the fact
                      is an MPI code
            #) See p325 (10/1/14) of Daily 2014 Notes for discussion of how to compute the
               relative solar metallicity
    FUTURE:
    """
    ### Python 3 required
    if(sys.version_info[0] != 3):
        exit_with_error("ERROR!!! Runs with python3, NOT python-{}\n\n".format(
                    sys.version_info[0]))

    ######### Get Command Line Options ##########
    if(len(sys.argv) != 2 and len(sys.argv) != 3 and len(sys.argv) != 4):
        print_help(1)
    elif(len(sys.argv) == 2 and (sys.argv[1] == "--help" or sys.argv[1] == "-help")):
        print_help(0)
    elif(len(sys.argv) == 3):
        path=sys.argv[1]
        option=sys.argv[2]
    elif(len(sys.argv) == 4):
        path=sys.argv[1]
        option=sys.argv[2]
        nGrid =int(sys.argv[3])
    else:
        exit_with_error("ERROR!!! Invalid CL options : {}".format(sys.argv))
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

    ## Temporarily comment out for expediency
    ###### Let's try to plot the mass density ######
    if(option == "plot"):
        array = np.zeros([100,100],dtype=np.float64)
        # gas, dm and stars 
        for k in range(len(pL)):
            # Check logic here
            i = int(pL[k].posV[0]/header.boxSize * array.shape[0])
            j = int(pL[k].posV[1]/header.boxSize * array.shape[1])
            #print(i,j)
            array[i,j] += pL[k].mass
        
        array = array / (header.npartV[0] + header.npartV[4])
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
                array[i,j] += pL[k].metalMassD['C']
        
        array = array / (header.npartV[0] + header.npartV[4])
        array = np.log(array)
        fig, ax = plt.subplots()
        im = ax.imshow(array)
        plt.show()
      
    elif(option == "out_cmp_gas"):
        ###### Let's try to Replicate gadget2/128part-25mpc/snapshot_016/mnc2file_gas.txt ######
        print("Writing gas particles to output_gas.txt")
        fout = open("output_gas.txt", "w+")
        for p in pL:
            if(p.type == 0):
                fout.write("{:<.6e} {:<.6e} {:<.6e} {:<.6f} {:<.6f} {:<.6f} {:<.6f} {:<.6f} "
                           "{:<.6f}\n".format(p.rho, p.u, p.mass, p.metallicityD['C'],
                           p.metallicityD['O'], p.metallicityD['Ca'], p.metallicityD['Cr'],
                           p.metallicityD['Mn'], p.metallicityD['Fe']))

    elif(option == "out_cmp_star"):
        ##### Let's try to Replicate gadget2/128part-25mpc/snapshot_016/mnc2file_star.txt ######
        print("Writing gas particles to output_star.txt")
        fout = open("output_star.txt", "w+")
        for p in pL:
            if(p.type == 4):
                fout.write("{:<.6e} {:<.6e} {:<.6f} {:<.6f} {:<.6f} {:<.6f} {:<.6f} "
                           "{:<.6f}\n".format(p.mass, p.age, p.metallicityD['C'],
                           p.metallicityD['O'], p.metallicityD['Ca'], p.metallicityD['Cr'],
                           p.metallicityD['Mn'], p.metallicityD['Fe']))
    elif(option == "out_gas_matlab"):
        print("Outputting gas particles in a ({:}^3) array".format(nGrid))
        array = np.zeros([nGrid,nGrid,nGrid],dtype=np.float64)
        for p in pL:
            if(p.type == 0):
                i = int(p.posV[0] / header.boxSize * nGrid)
                j = int(p.posV[1] / header.boxSize * nGrid)
                k = int(p.posV[2] / header.boxSize * nGrid)
                array[i,j,k] = array[i,j,k] + p.mass
        #### Sanity check, array summedhere should == array in 'plot' ####
        ## When using all particles, they agree (as they should)
        # array = np.sum(array, axis=2)
        # array = np.log(array)
        # fig, ax = plt.subplots()
        # im = ax.imshow(array[:,:,0])
        # plt.show()
        #array = np.log(array+1)
        #### Plot histogram
        # fig, ax = plt.subplots()
        # (n, bins, drop) = ax.hist(np.ndarray.flatten(array), bins=30)
        # plt.show()
        #### Normalize so can be cast as int easily
        minVal = np.min(array[array>0])
        array  = (array/minVal).astype(int)


        outName = path.split("/")[-1]
        outName = outName.split(".")[0]
        outName = "{}_ngrid_{}_matlab_int.txt".format(outName,nGrid)
        output_array_as_matlab_int(Array=array, OutName=outName)
        

        
    
    print("\n\nEnded : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))
    sys.exit(0)

if __name__ == "__main__":
    main()
