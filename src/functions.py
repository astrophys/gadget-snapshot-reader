# AUTHOR : Ali Snedden
# LICENSE: MIT
# DATE   : 05/11/2022
# PURPOSE:
#   
# NOTES:
#
# TODO:
#
# DEBUG:
#
import sys
import struct
import numpy as np
from error import exit_with_error
from classes import HEADER_V2
from classes import PARTICLE
from classes import ELEMENT


def read_gadget2_snapshot(Path=None, Endian="little"):
    """
    ARGS:
        Path   : (str) path to snapshot file 
        Endian : (str) endianness
    DESCRIPTION:
        This function reads a snapshot file that was created by gadget2-ali, which 
        includes feedback (i.e. chemical enrichment) and star formation.
    RETURN:
    DEBUG:
    FUTURE:
    """
    floatSize = 4
    intSize   = 4
    edn = Endian
    fin = open(Path, "rb")
    buf    = fin.read(4)        # buffer
    headB  = fin.read(256)
    buf    = fin.read(4)        # buffer

    ########### Header ###########
    header = HEADER_V2(Bytes=headB, Endian="little")
    header.print()
    #physTime = get_phys_time_My(header)
    #print("physTime = {}".format(physTime))
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
        p = PARTICLE(Mass=0)    ### ERROR HERE!!! Seems to just be creating pointers to the same object
        pL.append(p)


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
            #print(pL[0].posV[0], i)
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
            

        ##### I think I'm missing chunks for TypeIa_SFR and TypeII_SFR, should be for star particles NOT gas.
    
        ########### Type II SN Rate ###########
        if(header.flagD['RII'] == True):
            buf    = fin.read(4)        # buffer
            offset = header.npartV[0] + header.npartV[1] + header.npartV[2] + header.npartV[3]
            for i in range(header.npartV[4]):
                bytes = fin.read(floatSize)
                rII   = struct.unpack(edn+"f", bytes)[0]
                pL[offset + i].rII = rII
            buf    = fin.read(4)        # buffer
            
        ########### Type Ia SN Rate ###########
        if(header.flagD['RIa'] == True):
            buf    = fin.read(4)        # buffer
            offset = header.npartV[0] + header.npartV[1] + header.npartV[2] + header.npartV[3]
            for i in range(header.npartV[4]):
                bytes = fin.read(floatSize)
                rIa   = struct.unpack(edn+"f", bytes)[0]
                pL[offset + i].rIa = rIa
            buf    = fin.read(4)        # buffer
            

        ########### Stellar age ###########
        if(header.flagD['sfr'] == True and header.flagD['stellarage'] == True):
            buf    = fin.read(4)        # buffer
            offset = header.npartV[0] + header.npartV[1] + header.npartV[2] + header.npartV[3]
            for i in range(header.npartV[4]):
                bytes = fin.read(floatSize)
                age   = struct.unpack(edn+"f", bytes)[0]
                pL[offset + i].age = header.timeMy - age     # How old the star is in My
            buf    = fin.read(4)        # buffer
            
    elemD, solarMassConst = read_abundances(header)
    ########### Metals Fraction ###########
    if(header.flagD['metals'] == True):
        # Get metal keys where it is True
        mKeyL = [k for k, v in header.metalMassD.items() if v == True]
        if(header.npartV[0] > 0 or header.npartV[4] > 0):
            # Only read metals where the flag in header is True
            for key in mKeyL:
                read_metal(PL=pL, Header=header, Symbol=key, ElemD=elemD,
                           SolarMassConst=solarMassConst, Fin=fin)
        else:
            exit_with_error("ERROR!!! 'flag_metals' enabled but no gas or stars present!")


    ###### Ensure that all the data is read and that the EOF has been reached ######
    #   https://stackoverflow.com/a/53792531/4021436
    curPos = fin.tell()
    eof    = fin.seek(0,2)
    lastByte = fin.read(1)
    if(curPos - eof != 0 or len(lastByte) != 0):
        exit_with_error("ERROR!!! EOF ({}) not reached! Current Position = {}, lastByte = {}".
                        format(eof,curPos, lastByte))
    fin.close()


    return(pL, header)





def read_abundances(Header=None):
    """
    ARGS:
        Header : HEADER obj, only used to use the elements we are tracking to compute
                 the SolarMassConst
    DESCRIPTION:
        This function reads a src/Solar_Abundances.txt and returns the 
        elements, atomic mass and abundance for each element.
    RETURN:
        elemD          : Dictionary of ELEMENT objects
        solarMassConst : The Sum of Atomic Mass * Number fraction_solar, used in computing
                         metallicities
    DEBUG:
    FUTURE:
    """
    path = "src/Solar_Abundances.txt"
    elemD   = dict()
    top    = 0                   # numerator
    bottom = 0                   # denominator
    fin = open(path, "r")
    for line in fin:
        # ignore comments
        if(line[0] == "#"):
            continue
        line = line.split()
        symbol = line[0]             # column 1, Periodic Table abrev. for elements
        aMass  = float(line[1])      # column 2, Atomic Mass
        abund  = float(line[2])      # column 3, (n_i/n_H)_solar = solar abundunce
        prod   = float(line[3])      # column 4, (atomic mass) * (n_i/n_H)_solar
        elemD[symbol] = ELEMENT(Symbol=symbol, AtomicMass=aMass, NziNh = abund)
        if(symbol == "H" or symbol == "He"):
            bottom += prod
        elif(symbol == "Li" or symbol == "Be" or symbol == "B"):
            continue
        elif(Header.metalMassD[symbol] == True):
            top += prod
            bottom += prod
        # Sanity check
        if(np.isclose(prod, aMass*abund, rtol=0.0001) == False):
            exit_with_error("ERROR!!! The columns don't match for {} = {} * {}"
                            "= {}".format(symbol, aMass, abund, prod))
    # compute massFracMetal     = m_metal / (m_H + m_He + m_metal)
    #for i in range(len(symbolL)):
    massFracMetal = top / bottom
    solarMassConst= bottom
    print("massFracMetal = {}".format(massFracMetal))
    print("solarMassConst= {}".format(solarMassConst))
    return(elemD, solarMassConst)






def read_metal(PL = None, Header=None, Symbol=None, FloatSize=4, ElemD=None,
               SolarMassConst=None, Fin=None, Endien="<"):
    """
    ARGS:
        PL = Particle List
        Header  = header
        Symbol  = (str) short (IUPAC) name of metal to read
        ElemD   = Dictionary of ELEMENTS for computing_metallicity()
        Fin     = Input file (already opened)
        Endien  = string, '<' == little, '>' == big
    DESCRIPTION:
        Generalizes a very repetitive chunk of my code
    RETURN:
    DEBUG:
    FUTURE:
    """
    edn = Endien
    #if(getattr(Header, "flag_" + Metal) == True):
    print("Reading {} Mass Frac".format(Symbol))
    buf    = Fin.read(4)        # buffer
    # Gas
    for i in range(Header.npartV[0]):
        bytes = Fin.read(FloatSize)
        mass  = struct.unpack(edn+"f", bytes)[0]
        #setattr(PL[i], Symbol, mass)
        PL[i].metalMassD[Symbol] =  mass
        # Regarding metallicity, we need to
        #   1. Compute our metal MASS fraction 
        #   2. Compute the metal mass fraction at solar metallicity
        #      --> This is why we need the SolarMassConst
        #   3. Divide our metal mass frac by that at solar metallicity 
        #   4. Take the log
        PL[i].metallicityD[Symbol] = np.log10((mass / PL[i].mass) * 1 / (ElemD[Symbol].nzinh * ElemD[Symbol].amass / SolarMassConst))

    # Stars
    offset = Header.npartV[0] + Header.npartV[1] + Header.npartV[2] + Header.npartV[3]
    for i in range(Header.npartV[4]):
        bytes = Fin.read(FloatSize)
        mass  = struct.unpack(edn+"f", bytes)[0]
        PL[offset + i].metalMassD[Symbol]   = mass
        # Regarding metallicity, we need to
        #   1. Compute our metal MASS fraction 
        #   2. Compute the metal mass fraction at solar metallicity
        #      --> This is why we need the SolarMassConst
        #   3. Divide our metal mass frac by that at solar metallicity 
        #   4. Take the log
        PL[offset + i].metallicityD[Symbol] = np.log10((mass / PL[offset + i].mass) * 1 / (ElemD[Symbol].nzinh * ElemD[Symbol].amass / SolarMassConst))
    buf    = Fin.read(4)        # buffer




def output_array_as_matlab_int(Array = None, OutName = None):
    """
    ARGS:
        Density_array   :   Density array,
        OutName         :   Output file name

    DESCRIPTION:
        Print out density array in the matlab int format.

    RETURN:
        N/A

    DEBUG:
        1. a) Output scaled [0, 255] value density_array (not scaling by z-axis!)
              from JD01_B_comE_16hr.lsm-Ch2-C1matlab_int.txt
           b) Output of Ray's JD01_B_comE_16hr.lsm-Ch2-C1float.vtk
           c) The output's were _identical_ (compared with diff)
           ---> I am correctly scaling the input data _and_ getting the x-y-z
                format correct!

        2. a) Output unscaled value density_array (not scaling by z-axis!)
              from JD01_B_comE_16hr.lsm-Ch2-C1matlab_int.txt
           b) Compared output to original file. Were IDENTICAL.
           ---> I am correctly reading in the original data file.

        CONCLUSION : output_array_as_matlab_int(),scale_array()
                     and read_matlab_int_file() must all be working
                     correctly, otherwise I would not have been able to
                     get the correct output.
    FUTURE:
    """
    data  = Array if Array is not None else exit_with_error("ERROR in Array!\n")
    dim   = data.shape
    fname = OutName if OutName is not None else exit_with_error("ERROR in OutName!\n")
    fout  = open(fname, "w+")
    fout.write("# vtk DataFile Version 1.0. next line is dim[]... Dx=change"
               "column, Dy=each row, Dz=every dim[0] rows\n")
    fout.write("%i %i %i\n"%(dim[0],dim[1],dim[2]))
    for k in range(0, dim[2]):
        for j in range(0, dim[1]):
            for i in range(0, dim[0]):
                fout.write("%i "%(data[i,j,k]))
            fout.write("\n")
    print("Finished outputting : %s"%(fname))


