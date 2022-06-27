# AUTHOR : Ali Snedden
# LICENSE: MIT
# DATE   : 4/14/22
# PURPOSE:
#
# NOTES:
# TODO:
# DEBUG:
from error import exit_with_error
#from functions import get_phys_time_My
from scipy.integrate import quad
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
        self.metalMassD = dict()

        ############
        # Carbon flag
        dh = 1
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['C'] = True # 
        elif(flag == 0):
            self.metalMassD['C'] = False # 
        h += dh

        ############
        # Nitrogen flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['N'] = True # 
        elif(flag == 0):
            self.metalMassD['N'] = False # 
        h += dh

        ############
        # Oxygen flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['O'] = True # 
        elif(flag == 0):
            self.metalMassD['O'] = False # 
        h += dh

        ############
        # Florine flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['F'] = True # 
        elif(flag == 0):
            self.metalMassD['F'] = False # 
        h += dh

        ############
        # Neon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Ne'] = True # 
        elif(flag == 0):
            self.metalMassD['Ne'] = False # 
        h += dh

        ############
        # Sodium Flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Na'] = True # 
        elif(flag == 0):
            self.metalMassD['Na'] = False # 
        h += dh

        ############
        # Magnesium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Mg'] = True # 
        elif(flag == 0):
            self.metalMassD['Mg'] = False # 
        h += dh

        ############
        # Aluminum flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Al'] = True # 
        elif(flag == 0):
            self.metalMassD['Al'] = False # 
        h += dh

        ############
        # Silicon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Si'] = True # 
        elif(flag == 0):
            self.metalMassD['Si'] = False # 
        h += dh

        ############
        # Phosphorus flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['P'] = True # 
        elif(flag == 0):
            self.metalMassD['P'] = False # 
        h += dh

        ############
        # Sulfur flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['S'] = True # 
        elif(flag == 0):
            self.metalMassD['S'] = False # 
        h += dh

        ############
        # Chlorine flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Cl'] = True # 
        elif(flag == 0):
            self.metalMassD['Cl'] = False # 
        h += dh

        ############
        # Argon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Ar'] = True # 
        elif(flag == 0):
            self.metalMassD['Ar'] = False # 
        h += dh

        ############
        # Potassium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['K'] = True # 
        elif(flag == 0):
            self.metalMassD['K'] = False # 
        h += dh

        ############
        # Calcium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Ca'] = True # 
        elif(flag == 0):
            self.metalMassD['Ca'] = False # 
        h += dh

        ############
        # Scandium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Sc'] = True # 
        elif(flag == 0):
            self.metalMassD['Sc'] = False # 
        h += dh

        ############
        # Titanium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Ti'] = True # 
        elif(flag == 0):
            self.metalMassD['Ti'] = False # 
        h += dh

        ############
        # Vanadium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['V'] = True # 
        elif(flag == 0):
            self.metalMassD['V'] = False # 
        h += dh

        ############
        # Chromium flag 
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Cr'] = True # 
        elif(flag == 0):
            self.metalMassD['Cr'] = False # 
        h += dh

        ############
        # Manganese flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Mn'] = True # 
        elif(flag == 0):
            self.metalMassD['Mn'] = False # 
        h += dh

        ############
        # Iron flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Fe'] = True # 
        elif(flag == 0):
            self.metalMassD['Fe'] = False # 
        h += dh

        ############
        # Cobalt flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Co'] = True # 
        elif(flag == 0):
            self.metalMassD['Co'] = False # 
        h += dh

        ############
        # Nickle flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Ni'] = True # 
        elif(flag == 0):
            self.metalMassD['Ni'] = False # 
        h += dh

        ############
        # Copper flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Cu'] = True # 
        elif(flag == 0):
            self.metalMassD['Cu'] = False # 
        h += dh

        ############
        # Zinc flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.metalMassD['Zn'] = True # 
        elif(flag == 0):
            self.metalMassD['Zn'] = False # 
        h += dh

        ### Do relevant computation
        self.timeMy     = get_phys_time_My(self)                   # Physical time in My
        #char fill[27];               // fills to 256 Bytes
        #print(vars(self))


    def print_header(self):
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
        print("{:>22} : {:<.3f}My".format("physTime", self.timeMy))
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
        if(self.metalMassD['C']):
            print("{:>22} : {} ".format("metalMassD['C']", self.metalMassD['C']))
        if(self.metalMassD['N']):
            print("{:>22} : {} ".format("metalMassD['N']", self.metalMassD['N']))
        if(self.metalMassD['O']):
            print("{:>22} : {} ".format("metalMassD['O']", self.metalMassD['O']))
        if(self.metalMassD['F']):
            print("{:>22} : {} ".format("metalMassD['F']", self.metalMassD['F']))
        if(self.metalMassD['Ne']):
            print("{:>22} : {} ".format("metalMassD['Ne']", self.metalMassD['Ne']))
        if(self.metalMassD['Na']):
            print("{:>22} : {} ".format("metalMassD['Na']", self.metalMassD['Na']))
        if(self.metalMassD['Mg']):
            print("{:>22} : {} ".format("metalMassD['Mg']", self.metalMassD['Mg']))
        if(self.metalMassD['Al']):
            print("{:>22} : {} ".format("metalMassD['Al']", self.metalMassD['Al']))
        if(self.metalMassD['Si']):
            print("{:>22} : {} ".format("metalMassD['Si']", self.metalMassD['Si']))
        if(self.metalMassD['P']):
            print("{:>22} : {} ".format("metalMassD['P']", self.metalMassD['P']))
        if(self.metalMassD['S']):
            print("{:>22} : {} ".format("metalMassD['S']", self.metalMassD['S']))
        if(self.metalMassD['Cl']):
            print("{:>22} : {} ".format("metalMassD['Cl']", self.metalMassD['Cl']))
        if(self.metalMassD['Ar']):
            print("{:>22} : {} ".format("metalMassD['Ar']", self.metalMassD['Ar']))
        if(self.metalMassD['K']):
            print("{:>22} : {} ".format("metalMassD['K']", self.metalMassD['K']))
        if(self.metalMassD['Ca']):
            print("{:>22} : {} ".format("metalMassD['Ca']", self.metalMassD['Ca']))
        if(self.metalMassD['Sc']):
            print("{:>22} : {} ".format("metalMassD['Sc']", self.metalMassD['Sc']))
        if(self.metalMassD['Ti']):
            print("{:>22} : {} ".format("metalMassD['Ti']", self.metalMassD['Ti']))
        if(self.metalMassD['V']):
            print("{:>22} : {} ".format("metalMassD['V']", self.metalMassD['V']))
        if(self.metalMassD['Cr']):
            print("{:>22} : {} ".format("metalMassD['Cr']", self.metalMassD['Cr']))
        if(self.metalMassD['Mn']):
            print("{:>22} : {} ".format("metalMassD['Mn']", self.metalMassD['Mn']))
        if(self.metalMassD['Fe']):
            print("{:>22} : {} ".format("metalMassD['Fe']", self.metalMassD['Fe']))
        if(self.metalMassD['Co']):
            print("{:>22} : {} ".format("metalMassD['Co']", self.metalMassD['Co']))
        if(self.metalMassD['Ni']):
            print("{:>22} : {} ".format("metalMassD['Ni']", self.metalMassD['Ni']))
        if(self.metalMassD['Cu']):
            print("{:>22} : {} ".format("metalMassD['Cu']", self.metalMassD['Cu']))
        if(self.metalMassD['Zn']):
            print("{:>22} : {} ".format("metalMassD['Zn']", self.metalMassD['Zn']))







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
        self.metalMassD = dict()
        self.metallicityD= dict()
        # These are masses, e.g. self.C is the mass of carbon in the part
        self.metalMassD['C']   = -1   
        self.metalMassD['N']   = -1
        self.metalMassD['O']   = -1
        self.metalMassD['F']   = -1
        self.metalMassD['Ne']  = -1
        self.metalMassD['Na']  = -1
        self.metalMassD['Ma']  = -1
        self.metalMassD['Al']  = -1
        self.metalMassD['Si'] = -1
        self.metalMassD['P']  = -1
        self.metalMassD['S']  = -1
        self.metalMassD['Cl'] = -1
        self.metalMassD['Ar'] = -1
        self.metalMassD['K']  = -1
        self.metalMassD['Ca'] = -1
        self.metalMassD['Sc'] = -1
        self.metalMassD['Ti'] = -1
        self.metalMassD['V']  = -1
        self.metalMassD['Cr']  = -1
        self.metalMassD['Mn']  = -1
        self.metalMassD['Fe']  = -1
        self.metalMassD['Co']  = -1
        self.metalMassD['Ni']  = -1
        self.metalMassD['Cu']  = -1
        self.metalMassD['Zn']  = -1
        # Recall that metallicity is defined as log10[n_zi / n_H] STOPPED HERE!!
        self.metallicityD['C']   = -1   
        self.metallicityD['N']   = -1
        self.metallicityD['O']   = -1
        self.metallicityD['F']   = -1
        self.metallicityD['Ne']  = -1
        self.metallicityD['Na']  = -1
        self.metallicityD['Ma']  = -1
        self.metallicityD['Al']  = -1
        self.metallicityD['Si'] = -1
        self.metallicityD['P']  = -1
        self.metallicityD['S']  = -1
        self.metallicityD['Cl'] = -1
        self.metallicityD['Ar'] = -1
        self.metallicityD['K']  = -1
        self.metallicityD['Ca'] = -1
        self.metallicityD['Sc'] = -1
        self.metallicityD['Ti'] = -1
        self.metallicityD['V']  = -1
        self.metallicityD['Cr']  = -1
        self.metallicityD['Mn']  = -1
        self.metallicityD['Fe']  = -1
        self.metallicityD['Co']  = -1
        self.metallicityD['Ni']  = -1
        self.metallicityD['Cu']  = -1
        self.metallicityD['Zn']  = -1



class ELEMENT:
    def __init__(self, Symbol=None, AtomicMass=None, NziNh = None):
        """
        ARGS:
            Symbol      = (str) IUPAC element symbol
            AtomicMass  = (float) Atomic Mass
            NziNh       = (float) Number Fraction of element 'Z' per 'H' at solar
                                  metallicity
        DESCRIPTION:
        RETURN:
        DEBUG:
        FUTURE:
            1. Make this less ugly
        """
        self.symbol = Symbol
        self.amass  = float(AtomicMass)
        self.nzinh  = float(NziNh)
        



def get_phys_time_My(Header=None):
    """
    ARGS:
        Header : HEADER obj, only used to use the elements we are tracking to compute
                 the SolarMassConst
    DESCRIPTION:
        This function integrates Friedman's eqn where we know H(a) = adot/a.
        Separate and integrate

                da/dt = a * H(a)  ---->  dt = da / (a * H(a)).

        We take H(a) from Carrol & Ostlie (2nd) eqn 29.128. Recall Omega0 = Omega_matter.
    RETURN:
        Time in MegaYears since the Big Bang
    NOTES:
        1. Omega_k = 1 - All.Omega0 - All.OmegaLambda, Omega_r = 0.
        2. HUBBLE is in (h/sec), All.HubbleParam is little 'h'
        3. Taken almost exclusively from mnc2file_sph/src/functions.c:phys_time_integrand()
    DEBUG:
    FUTURE:
    """
    def integrand(a, Header):
        """
        ARGS:
            a : (float) scale factor
            aF: (float) Final (current) value of scale factor
            Header : HEADER obj, only used to use the elements we are tracking to compute
                     the SolarMassConst
        DESCRIPTION:
        RETURN:
        DEBUG:
        FUTURE:
        """
        h = 0.702               # 'little h'
        H = 3.2407789e-18       # (h/s) mncfile_sph/src/allvars.h : HUBBLE
        return 1.0 / (H * h) * pow(Header.Omega0 / a + Header.OmegaLambda * pow(a, 2.0) + (1 - Header.Omega0 - Header.OmegaLambda), -0.5);
        #return a

    secPerMy = 3.1557e13
    # Integrate from a=0 (z~inf) to current Header.time (i.e. scale factor)
    time = quad(integrand, 0, Header.time, args=(Header))
    return(time[0]/secPerMy)
    #return(quad(integrand, 0, 3, args=(Header)))
