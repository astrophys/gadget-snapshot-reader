# AUTHOR : Ali Snedden
# LICENSE: MIT
# DATE   : 4/14/22
# PURPOSE:
#
# NOTES:
# TODO:
# DEBUG:
from error import exit_with_error
import struct
import numpy as np


class HEADER_V2:
    def __init__(self, Bytes = None, Endian = "little"):
        """
        ARGS:
        DESCRIPTION:
        RETURN:
        DEBUG:
        FUTURE:
            1. Make this less ugly
        """
        h  = 0          # Header index
        dh = 4          # number of bytes
        if(Endian.lower() == "little"):
            edn = "<"
            Endian = "little"
        elif(Endian.lower() == "big"):
            edn = ">"
            Endian = "big"
        else:
            exit_with_error("ERROR!!! Incorrect option ({}) for Endianess".format(Endian))
        self.endian = Endian

        ############
        # Number of Particles - [Gas, DM, X, Star, X, X]
        self.npartV     = np.zeros(6, dtype=np.int64)           # int
        for i in range(6):
            # I worry about type casting here
            self.npartV[i] = int.from_bytes(Bytes[h:h+dh], Endian)
            #self.npartV[i] = struct.unpack(Bytes[h:h+dh],"little")
            h += dh

        ############
        # Mass of each group  - [Gas, DM, X, Star, X, X]
        self.massV      = np.zeros(6, dtype=np.float64)           # doubles
        dh = 8
        for i in range(6):
            # If 0.0, then mass stored explicitlyin mass block
            self.massV[i] = struct.unpack(edn+'d', Bytes[h:h+dh])[0]
            h += dh

        ############
        # a 
        self.time       = struct.unpack(edn+'d', Bytes[h:h+dh])[0] # double, I think a
        h += dh

        ############
        # redshift
        self.z          = struct.unpack(edn+'d', Bytes[h:h+dh])[0] # double, redshift
        h += dh

        ############
        # SFR flag 
        self.flagD  = dict()
        dh = 4
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        h += dh
        if(flag == 1):
            self.flagD['sfr'] = True          # bool, Star formation 
        elif(flag == 0):
            self.flagD['sfr'] = False         # bool, Star formation 
            
        ############
        # Feedback flag 
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        h += dh
        if(flag == 1):
            self.flagD['feed'] = True          # bool, Feedback
        elif(flag == 0):
            self.flagD['feed'] = False         # bool, Feedback

        ############
        # int, Num particles of each type in this snapshot
        self.npartTotal = np.zeros(6, dtype=np.int64)
        for i in range(6):
            self.npartTotal[i] = struct.unpack(edn+'I', Bytes[h:h+dh])[0]
            h += dh

        ############
        # Cooling flag 
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        h += dh
        if(flag == 1):
            self.flagD['cooling'] = True        # bool, Cooling included
        elif(flag == 0):
            self.flagD['cooling'] = False       # bool, Cooling included

        ############
        # N files in multi-file snapshot
        self.num_files  = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # int
        h += dh

        ############
        # Box-size of simulation
        dh = 8
        self.boxSize    = struct.unpack(edn+'d', Bytes[h:h+dh])[0]             # double
        h += dh

        ############
        # matter density in units of critical density
        self.Omega0     = struct.unpack(edn+'d', Bytes[h:h+dh])[0]             # double 
        h += dh

        ############
        # cosmological constant parameter
        self.OmegaLambda= struct.unpack(edn+'d', Bytes[h:h+dh])[0]             # double 
        h += dh

        ############
        # Hubble parameter in units of 100 km/sec/Mpc
        self.HubbleParam= struct.unpack(edn+'d', Bytes[h:h+dh])[0]             # double 
        h += dh

        ############
        # File contains formation times of star particles
        dh = 4
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flagD['stellarage'] = True         # bool
        elif(flag == 0):
            self.flagD['stellarage'] = False        # bool
        h += dh

        ############
        # File contains metallicity for gas and star particles
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flagD['metals'] = True              # bool
        elif(flag == 0):
            self.flagD['metals'] = False             # bool
        h += dh

        ############
        # High word of the total number of particles of each type
        self.npartTotalHighWord = np.zeros(6, dtype=np.int64)
        for i in range(6):
            self.npartTotalHighWord[i] = struct.unpack(edn+'I', Bytes[h:h+dh])[0]
            h += dh

        ############
        # flags that IC-file contains entropy instead of u 
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flagD['entropy_not_u'] = True # 
        elif(flag == 0):
            self.flagD['entropy_not_u'] = False # 
        h += dh
        

        ############
        # flags Type II SN were present
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flagD['RII'] = True # 
        elif(flag == 0):
            self.flagD['RII'] = False # 
        h += dh


        ############
        # flags Type Ia SN were present
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flagD['RIa'] = True # 
        elif(flag == 0):
            self.flagD['RIa'] = False # 
        h += dh

        ############
        # Recall that the words used by C are 8bytes in this case, looking at
        # gadget2-ali/src/allvars.h (simulated by src/struct-size.c), the size is 208bytes
        # where there is a 4byte buffer b/c of flag_R1a. This means that h=208. If you 
        # include the 4byte buffer below, you get the 'correct' value.  HOWEVER 
        # on the CRC, ~/Lab/phillips/pscratch/me/128/25Mpc/Makefile indicates that
        # Carbon, Oxygen, Calcium, Chromium, Manganese, Iron). Those flags ONLY work
        # if h=204 and the below offset is excluded. Why this is the case, I have no friggin' 
        # idea, other than there must be a bug somewhere.
        #h += 4
        self.metalD = dict()

        ############
        # Carbon flag
        dh = 1
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['C'] = True # 
        elif(flag == 0):
            self.metalD['C'] = False # 
        h += dh

        ############
        # Nitrogen flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['N'] = True # 
        elif(flag == 0):
            self.metalD['N'] = False # 
        h += dh

        ############
        # Oxygen flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['O'] = True # 
        elif(flag == 0):
            self.metalD['O'] = False # 
        h += dh

        ############
        # Florine flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['F'] = True # 
        elif(flag == 0):
            self.metalD['F'] = False # 
        h += dh

        ############
        # Neon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Ne'] = True # 
        elif(flag == 0):
            self.metalD['Ne'] = False # 
        h += dh

        ############
        # Sodium Flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Na'] = True # 
        elif(flag == 0):
            self.metalD['Na'] = False # 
        h += dh

        ############
        # Magnesium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Mg'] = True # 
        elif(flag == 0):
            self.metalD['Mg'] = False # 
        h += dh

        ############
        # Aluminum flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Al'] = True # 
        elif(flag == 0):
            self.metalD['Al'] = False # 
        h += dh

        ############
        # Silicon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Si'] = True # 
        elif(flag == 0):
            self.metalD['Si'] = False # 
        h += dh

        ############
        # Phosphorus flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['P'] = True # 
        elif(flag == 0):
            self.metalD['P'] = False # 
        h += dh

        ############
        # Sulfur flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['S'] = True # 
        elif(flag == 0):
            self.metalD['S'] = False # 
        h += dh

        ############
        # Chlorine flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Cl'] = True # 
        elif(flag == 0):
            self.metalD['Cl'] = False # 
        h += dh

        ############
        # Argon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Ar'] = True # 
        elif(flag == 0):
            self.metalD['Ar'] = False # 
        h += dh

        ############
        # Potassium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['K'] = True # 
        elif(flag == 0):
            self.metalD['K'] = False # 
        h += dh

        ############
        # Calcium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Ca'] = True # 
        elif(flag == 0):
            self.metalD['Ca'] = False # 
        h += dh

        ############
        # Scandium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Sc'] = True # 
        elif(flag == 0):
            self.metalD['Sc'] = False # 
        h += dh

        ############
        # Titanium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Ti'] = True # 
        elif(flag == 0):
            self.metalD['Ti'] = False # 
        h += dh

        ############
        # Vanadium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['V'] = True # 
        elif(flag == 0):
            self.metalD['V'] = False # 
        h += dh

        ############
        # Chromium flag 
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Cr'] = True # 
        elif(flag == 0):
            self.metalD['Cr'] = False # 
        h += dh

        ############
        # Manganese flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Mn'] = True # 
        elif(flag == 0):
            self.metalD['Mn'] = False # 
        h += dh

        ############
        # Iron flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Fe'] = True # 
        elif(flag == 0):
            self.metalD['Fe'] = False # 
        h += dh

        ############
        # Cobalt flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Co'] = True # 
        elif(flag == 0):
            self.metalD['Co'] = False # 
        h += dh

        ############
        # Nickle flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Ni'] = True # 
        elif(flag == 0):
            self.metalD['Ni'] = False # 
        h += dh

        ############
        # Copper flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Cu'] = True # 
        elif(flag == 0):
            self.metalD['Cu'] = False # 
        h += dh

        ############
        # Zinc flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalD['Zn'] = True # 
        elif(flag == 0):
            self.metalD['Zn'] = False # 
        h += dh

        #char fill[27];               // fills to 256 Bytes
        #print(vars(self))


    def print(self):
        """
        ARGS:
        DESCRIPTION:
            This pretty print of the header
        RETURN:
        DEBUG:
        FUTURE:
        """
        print("###########################################################################")
        print("############################ Snapshot Header ##############################")
        print("###########################################################################")
        # Endian
        print("{:>22} : {}".format("Endian", self.endian))

        # Number of Particles
        print("{:>22} : ".format("Particles"))
        print("{:>22}   {} gas".format("", self.npartV[0]))
        print("{:>22}   {} dm".format("", self.npartV[1]))
        print("{:>22}   {} ".format("", self.npartV[2]))
        print("{:>22}   {} ".format("", self.npartV[3]))
        print("{:>22}   {} stars".format("", self.npartV[4]))
        print("{:>22}   {} ".format("", self.npartV[5]))


        # Mass
        print("{:>22} : ".format("Mass (10^10 M_sun)"))
        print("{:>22}   {} gas".format("", self.massV[0]))
        print("{:>22}   {} dm".format("", self.massV[1]))
        print("{:>22}   {} ".format("", self.massV[2]))
        print("{:>22}   {} ".format("", self.massV[3]))
        print("{:>22}   {} ".format("", self.massV[4]))
        print("{:>22}   {} ".format("", self.massV[5]))

        # Total number of particles 
        print("{:>22} : {} ".format("npartTotal", self.npartTotal))

        # Number of files in multi-file snapshot
        print("{:>22} : {} ".format("N files/snap", self.num_files))

        # High word of the total number of particles of each type
        print("{:>22} : ".format("Total High Word"))
        print("{:>22} : {} gas".format("", self.npartTotalHighWord[0]))
        print("{:>22} : {} dm".format("", self.npartTotalHighWord[1]))
        print("{:>22} : {} ".format("", self.npartTotalHighWord[2]))
        print("{:>22} : {} ".format("", self.npartTotalHighWord[3]))
        print("{:>22} : {} ".format("", self.npartTotalHighWord[4]))
        print("{:>22} : {} ".format("", self.npartTotalHighWord[5]))

        # flags that IC-file contains entropy instead of U
        print("{:>22} : {} ".format("", self.flagD['entropy_not_u']))
            
        # Cosmological Parameters
        print("\n\n############ Cosmological Parameters ############")
        print("{:>22} : {} ".format("Omega0", self.Omega0))
        print("{:>22} : {} ".format("OmegaLambda", self.OmegaLambda))
        print("{:>22} : {} ".format("HubbleParam", self.HubbleParam))
        # Scale factor - a
        print("{:>22} : {}".format("a", self.time))
        # Redshift - z
        print("{:>22} : {}".format("redshift", self.z))
        # Box-size
        print("{:>22} : {} Mpc/h".format("Box-size", self.boxSize))

        print("\n\n############ Feedback Parameters ############")
        # Cooling
        print("{:>22} : {} ".format("cooling", self.flagD['cooling']))
        # Star formation flags
        print("{:>22} : {}".format("flag_sfr", self.flagD['sfr']))
        print("{:>22} : {}".format("flag_feed", self.flagD['feed']))

        # Stellar Age
        print("{:>22} : {} ".format("StellarAge", self.flagD['stellarage']))

        # Metals
        print("{:>22} : {} \n".format("flag_metals", self.flagD['metals']))

        # Type II SN
        print("{:>22} : {} ".format("Type II SNR", self.flagD['RII']))

        # Type Ia SN
        print("{:>22} : {} ".format("Type Ia SNR", self.flagD['RIa']))

        # Metals
        if(self.metalD['C']):
            print("{:>22} : {} ".format("metalD['C']", self.metalD['C']))
        if(self.metalD['N']):
            print("{:>22} : {} ".format("metalD['N']", self.metalD['N']))
        if(self.metalD['O']):
            print("{:>22} : {} ".format("metalD['O']", self.metalD['O']))
        if(self.metalD['F']):
            print("{:>22} : {} ".format("metalD['F']", self.metalD['F']))
        if(self.metalD['Ne']):
            print("{:>22} : {} ".format("metalD['Ne']", self.metalD['Ne']))
        if(self.metalD['Na']):
            print("{:>22} : {} ".format("metalD['Na']", self.metalD['Na']))
        if(self.metalD['Mg']):
            print("{:>22} : {} ".format("metalD['Mg']", self.metalD['Mg']))
        if(self.metalD['Al']):
            print("{:>22} : {} ".format("metalD['Al']", self.metalD['Al']))
        if(self.metalD['Si']):
            print("{:>22} : {} ".format("metalD['Si']", self.metalD['Si']))
        if(self.metalD['P']):
            print("{:>22} : {} ".format("metalD['P']", self.metalD['P']))
        if(self.metalD['S']):
            print("{:>22} : {} ".format("metalD['S']", self.metalD['S']))
        if(self.metalD['Cl']):
            print("{:>22} : {} ".format("metalD['Cl']", self.metalD['Cl']))
        if(self.metalD['Ar']):
            print("{:>22} : {} ".format("metalD['Ar']", self.metalD['Ar']))
        if(self.metalD['K']):
            print("{:>22} : {} ".format("metalD['K']", self.metalD['K']))
        if(self.metalD['Ca']):
            print("{:>22} : {} ".format("metalD['Ca']", self.metalD['Ca']))
        if(self.metalD['Sc']):
            print("{:>22} : {} ".format("metalD['Sc']", self.metalD['Sc']))
        if(self.metalD['Ti']):
            print("{:>22} : {} ".format("metalD['Ti']", self.metalD['Ti']))
        if(self.metalD['V']):
            print("{:>22} : {} ".format("metalD['V']", self.metalD['V']))
        if(self.metalD['Cr']):
            print("{:>22} : {} ".format("metalD['Cr']", self.metalD['Cr']))
        if(self.metalD['Mn']):
            print("{:>22} : {} ".format("metalD['Mn']", self.metalD['Mn']))
        if(self.metalD['Fe']):
            print("{:>22} : {} ".format("metalD['Fe']", self.metalD['Fe']))
        if(self.metalD['Co']):
            print("{:>22} : {} ".format("metalD['Co']", self.metalD['Co']))
        if(self.metalD['Ni']):
            print("{:>22} : {} ".format("metalD['Ni']", self.metalD['Ni']))
        if(self.metalD['Cu']):
            print("{:>22} : {} ".format("metalD['Cu']", self.metalD['Cu']))
        if(self.metalD['Zn']):
            print("{:>22} : {} ".format("metalD['Zn']", self.metalD['Zn']))







class PARTICLE:
    def __init__(self, PosV=None, Mass=None, Type=None, VelV=None, ID=None,
                 U=None, Rho=None, Hsml=None):
        """
        ARGS:
            PosL = 
            Mass = 
        DESCRIPTION:
        RETURN:
        DEBUG:
        FUTURE:
            1. Make this less ugly
        """
        if(PosV is None):
            self.posV = np.zeros(3)
        else:
            self.posV = PosV
        self.mass = Mass
        self.type = Type
        if(VelV is None):
            self.vV = np.zeros(3)
        else:
            self.vV = VelV
        self.id   = ID
        self.u    = U
        self.rho  = Rho
        self.hsml = Hsml
        self.rII  = -1
        self.rIa  = -1
        self.u    = VelV
        self.age  = -1
        self.metalD = dict()
        # These are masses, e.g. self.C is the mass of carbon in the part
        self.metalD['C']   = -1   
        self.metalD['N']   = -1
        self.metalD['O']   = -1
        self.metalD['F']   = -1
        self.metalD['Ne']  = -1
        self.metalD['Na']  = -1
        self.metalD['Ma']  = -1
        self.metalD['Al']  = -1
        self.metalD['Si'] = -1
        self.metalD['P']  = -1
        self.metalD['S']  = -1
        self.metalD['Cl'] = -1
        self.metalD['Ar'] = -1
        self.metalD['K']  = -1
        self.metalD['Ca'] = -1
        self.metalD['Sc'] = -1
        self.metalD['Ti'] = -1
        self.metalD['V']  = -1
        self.metalD['Cr']  = -1
        self.metalD['Mn']  = -1
        self.metalD['Fe']  = -1
        self.metalD['Co']  = -1
        self.metalD['Ni']  = -1
        self.metalD['Cu']  = -1
        self.metalD['Zn']  = -1
        # Recall that metallicity is defined as log10[n_zi / n_H] STOPPED HERE!!
        #self.C   = -1   
        #self.N   = -1
        #self.O   = -1
        #self.F   = -1
        #self.Ne  = -1
        #self.Na  = -1
        #self.Ma  = -1
        #self.Al  = -1
        #self.Si  = -1
        #self.P   = -1
        #self.S   = -1
        #self.Cl  = -1
        #self.Ar  = -1
        #self.K   = -1
        #self.Ca  = -1
        #self.Sc  = -1
        #self.Ti  = -1
        #self.V   = -1
        #self.Cr  = -1
        #self.Mn  = -1
        #self.Fe  = -1
        #self.Co  = -1
        #self.Ni  = -1
        #self.Cu  = -1
        #self.Zn  = -1


class ELEMENT:
    def __init__(self, AtomMass=None, SolarMetalNormConstant=1.3765e-02):
        """
        ARGS:
            PosL = 
            Mass = 
        DESCRIPTION:
        RETURN:
        DEBUG:
        FUTURE:
            1. Make this less ugly
        """



def read_metal(PL = None, Header=None, Short=None, FloatSize=4, Fin=None,
               Endien="<"):
    """
    ARGS:
        PL = Particle List
        Header = header
        Short  = (str) short (IUPAC) name of metal to read
        Fin    = Input file (already opened)
        Endien = string, '<' == little, '>' == big
    DESCRIPTION:
        Generalizes a very repetitive chunk of my code
    RETURN:
    DEBUG:
    FUTURE:
    """
    edn = Endien
    #if(getattr(Header, "flag_" + Metal) == True):
    print("Reading {} Mass Frac".format(Short))
    buf    = Fin.read(4)        # buffer
    # Gas 
    for i in range(Header.npartV[0]):
        bytes = Fin.read(FloatSize)
        mass  = struct.unpack(edn+"f", bytes)[0]
        setattr(PL[i], Short, mass)

    # Stars
    offset = Header.npartV[0] + Header.npartV[1] + Header.npartV[2] + Header.npartV[3]
    for i in range(Header.npartV[4]):
        bytes = Fin.read(FloatSize)
        mass  = struct.unpack(edn+"f", bytes)[0]
        setattr(PL[offset + i], Short, mass)
    buf    = Fin.read(4)        # buffer
