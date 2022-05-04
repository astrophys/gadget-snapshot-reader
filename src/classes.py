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
        dh = 4
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        h += dh
        if(flag == 1):
            self.flag_sfr   = True          # bool, Star formation 
        elif(flag == 0):
            self.flag_sfr   = False         # bool, Star formation 
            
        ############
        # Feedback flag 
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        h += dh
        if(flag == 1):
            self.flag_feed  = True          # bool, Feedback
        elif(flag == 0):
            self.flag_feed  = False         # bool, Feedback

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
            self.flag_cooling = True        # bool, Cooling included
        elif(flag == 0):
            self.flag_cooling = False       # bool, Cooling included

        ############
        # N files in multi-file snapshot
        self.num_files  = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # int
        h += dh

        ############
        # Box-size of simulation
        dh = 8
        self.BoxSize    = struct.unpack(edn+'d', Bytes[h:h+dh])[0]             # double
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
            self.flag_stellarage = True         # bool
        elif(flag == 0):
            self.flag_stellarage = False        # bool
        h += dh

        ############
        # File contains metallicity for gas and star particles
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flag_metals= True              # bool
        elif(flag == 0):
            self.flag_metals= False             # bool
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
            self.flag_entropy_not_u = True # 
        elif(flag == 0):
            self.flag_entropy_not_u =False # 
        h += dh
        

        ############
        # flags Type II SN were present
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flag_RII = True # 
        elif(flag == 0):
            self.flag_RII =False # 
        h += dh


        ############
        # flags Type Ia SN were present
        flag = struct.unpack(edn+'i', Bytes[h:h+dh])[0] # double, redshift
        if(flag == 1):
            self.flag_RIa= True # 
        elif(flag == 0):
            self.flag_RIa=False # 
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

        ############
        # Carbon flag
        dh = 1
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Carbon = True # 
        elif(flag == 0):
            self.flag_Carbon =False # 
        h += dh

        ############
        # Nitrogen flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Nitrogen = True # 
        elif(flag == 0):
            self.flag_Nitrogen =False # 
        h += dh

        ############
        # Oxygen flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Oxygen    = True # 
        elif(flag == 0):
            self.flag_Oxygen    =False # 
        h += dh

        ############
        # Florine flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Florine       = True # 
        elif(flag == 0):
            self.flag_Florine       =False # 
        h += dh

        ############
        # Neon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Neon = True # 
        elif(flag == 0):
            self.flag_Neon =False # 
        h += dh

        ############
        # Sodium Flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Sodium = True # 
        elif(flag == 0):
            self.flag_Sodium =False # 
        h += dh

        ############
        # Magnesium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Magnesium = True # 
        elif(flag == 0):
            self.flag_Magnesium =False # 
        h += dh

        ############
        # Aluminum flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Aluminum  = True # 
        elif(flag == 0):
            self.flag_Aluminum  =False # 
        h += dh

        ############
        # Silicon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Silicon       = True # 
        elif(flag == 0):
            self.flag_Silicon       =False # 
        h += dh

        ############
        # Phosphorus flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Phosphorus    = True # 
        elif(flag == 0):
            self.flag_Phosphorus    =False # 
        h += dh

        ############
        # Sulfur flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Sulfur        = True # 
        elif(flag == 0):
            self.flag_Sulfur        =False # 
        h += dh

        ############
        # Chlorine flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Chlorine      = True # 
        elif(flag == 0):
            self.flag_Chlorine      =False # 
        h += dh

        ############
        # Argon flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Argon         = True # 
        elif(flag == 0):
            self.flag_Argon         =False # 
        h += dh

        ############
        # Potassium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Potassium     = True # 
        elif(flag == 0):
            self.flag_Potassium     =False # 
        h += dh

        ############
        # Calcium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Calcium       = True # 
        elif(flag == 0):
            self.flag_Calcium       =False # 
        h += dh

        ############
        # Scandium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Scandium      = True # 
        elif(flag == 0):
            self.flag_Scandium      =False # 
        h += dh

        ############
        # Titanium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Titanium      = True # 
        elif(flag == 0):
            self.flag_Titanium      =False # 
        h += dh

        ############
        # Vanadium flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Vanadium      = True # 
        elif(flag == 0):
            self.flag_Vanadium      =False # 
        h += dh

        ############
        # Chromium flag 
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Chromium      = True # 
        elif(flag == 0):
            self.flag_Chromium      =False # 
        h += dh

        ############
        # Manganese flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Manganese     = True # 
        elif(flag == 0):
            self.flag_Manganese     =False # 
        h += dh

        ############
        # Iron flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Iron          = True # 
        elif(flag == 0):
            self.flag_Iron          =False # 
        h += dh

        ############
        # Cobalt flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Cobalt        = True # 
        elif(flag == 0):
            self.flag_Cobalt        =False # 
        h += dh

        ############
        # Nickle flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Nickel        = True # 
        elif(flag == 0):
            self.flag_Nickel        =False # 
        h += dh

        ############
        # Copper flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Copper        = True # 
        elif(flag == 0):
            self.flag_Copper        =False # 
        h += dh

        ############
        # Zinc flag
        flag = int(struct.unpack(edn+'c', Bytes[h:h+dh])[0]) # 
        if(flag == 1):
            self.flag_Zinc          = True # 
        elif(flag == 0):
            self.flag_Zinc          =False # 
        h += dh

        #char fill[27];               // fills to 256 Bytes
        print(vars(self))



class PARTICLE:
    def __init__(self, PosV = np.zeros(3), Mass=None, Type=None, VelV=np.zeros(3), ID=None,
                 U=None, Rho=None, Hsml=None, SFR=None):
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
        self.posV = PosV
        self.mass = Mass
        self.type = Type
        self.vV   = VelV
        self.id   = ID
        self.u    = U
        self.rho  = Rho
        self.hsml = Hsml
        self.sfr  = SFR
        self.u    = VelV
        self.age  = -1
        self.Cf   = -1
        self.Nf   = -1
        self.Of   = -1
        self.Ff   = -1
        self.Nef  = -1
        self.Naf  = -1
        self.Maf  = -1
        self.Alf  = -1
        self.Sif  = -1
        self.Pf   = -1
        self.Sf   = -1
        self.Clf  = -1
        self.Arf  = -1
        self.Kf   = -1
        self.Caf  = -1
        self.Scf  = -1
        self.Tif  = -1
        self.Vf   = -1
        self.Crf  = -1
        self.Mnf  = -1
        self.Fef  = -1
        self.Cof  = -1
        self.Nif  = -1
        self.Cuf  = -1
        self.Znf  = -1



def read_metal(PL = None, Header=None, Metal=None, Short=None, FloatSize=4, Fin=None,
               Endien="<"):
    """
    ARGS:
        PL = Particle List
        Header = header
        Metal  = (str) full name of metal to read
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
    if(getattr(Header, "flag_" + Metal) == True):
        print("Reading {} Mass Frac".format(Metal))
        buf    = Fin.read(4)        # buffer
        # Gas 
        for i in range(Header.npartV[0]):
            bytes = Fin.read(FloatSize)
            frac  = struct.unpack(edn+"f", bytes)[0]
            setattr(PL[i], Short, frac)

        # Stars
        offset = Header.npartV[0] + Header.npartV[1] + Header.npartV[2] + Header.npartV[3]
        for i in range(Header.npartV[4]):
            bytes = Fin.read(FloatSize)
            frac  = struct.unpack(edn+"f", bytes)[0]
            setattr(PL[offset + i], Short, frac)
        buf    = Fin.read(4)        # buffer
