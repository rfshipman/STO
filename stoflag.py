###########
# StoFlag
###########
import numpy as np


class StoFlag:
    
    def __init__(self, bitOffset, flagname):
        '''
        Define flags for STO-2.  Create a flag with bitOffset and flagname.
        
            bitOffset:  Index (int) of bit:  Example, bitOffset =3 means the flag is set for the 4th bit.
            flagname:  String name describing the flag.
            
            Usage:
            flagvalue=StoFlag(bit offset ,'Flag Name')
            
            Example:
            MASTER=StoFlag(0,"Master")
            print(MASTER.getFlag)
            1
            
            To set a flag at for a channel or range of channels in x
            x[range]=MASTER.setFlag(x[range])
        '''
      
        try:
            bitOffset < 0 or bitOffset >= 32
        except:
            print ("Bit offset values shoud be between 0 and 32")
            
        self.index=bitOffset
        self.bit=2**bitOffset
        self.name=flagname
        self.check=2**bitOffset
        
        
    def getName(self):
        return self.name
    
    def getBitOffset(self):
        return self.index
    
    def getFlag(self):
        return self.bit
    
    def getIndex(self):
        return self.index
    
    def setFlag(self,mask):
        flag=self.bit
        #is mask an array?
        if hasattr(mask,"__len__"):
            flagv=[]
            for m in mask:
                flagv.append(flag | m)
            flagv=np.array(flagv)
            return flagv
        else:
            return flag | mask
    
    def isSet(self, mask): 
        #
        dtype=isinstance(mask, int)
        if dtype:
            test = (self.check & mask) != 0
            return test
        elif hasattr(mask,"__len__"):
            #
            #b=np.zeros(mask.__len__(),dtype=np.bool)
            b=[]
            for m in mask:
                 test=(self.check & m) != 0
                 b.append(test)
            b=np.array(b)
            return b
        else:
            print("Cannot Set the Flag of this object")
            

        
        
        
    