import numpy as np
import atomic_mass

class Atom:
    
    # Static variables #
    one_letter_elements = {"H","B","C","N","O","F","P","K","S","V","I","U","Y","W","I"}     # one letter elements
    two_letter_elements = {"He","Li","Be","Ne",   
                           "Na","Mg","Al","Si","Cl","Ar","Ca",   # the 2nd letter is in the lower-case.
                           "Sc","Ti","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                           "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Zr",
                           "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                           "Sb","Te","Xe","Cs","Ba","La","Ce","Pr","Nd",
                           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                           "Lu","Hf","Ta","Re","Os","Ir","Pt","Au","Hg",
                           "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
                           "Pa","Np","Pu","Am","Cm","Bk","Cf"} #
    def __init__(self,line_ext):
        self.pos = np.array([0,0,0],dtype='f')
        line_ext=line_ext.replace("\n","") # getting rid of \n and \r in the line
        line_ext=line_ext.replace("\r","") # getting rid of \n and \r in the line
        self.line = "%-80s\n" % line_ext
        self.pos[0] = float(self.line[30:38])
        self.pos[1] = float(self.line[38:46])
        self.pos[2] = float(self.line[46:54])
        self.name   = self.line[12:16]
        self.resname = self.line[17:21]
        self.index  = 1
        self.element = self.ExtractElement()
        self.bonded_atoms=[]
        
    def ApplyTransform(self,trans): # trans is a vtk.vtkTransfom()
        trans.TransformPoint(self.pos,self.pos)
        self.UpdateCrd()
    
    def SetPosition(self,ex_pos):
        self.pos = ex_pos
        self.UpdateCrd()
    
    def UpdateCrd(self):
        self.line = self.line[:30] + "%8.3f" % self.pos[0] + self.line[38:]
        self.line = self.line[:38] + "%8.3f" % self.pos[1] + self.line[46:]
        self.line = self.line[:46] + "%8.3f" % self.pos[2] + self.line[54:]
        ## fix element ##
        if not self.element:  
            pass
        
    def GetMolNameGMX(self):
        """Returns the molecules name with Gromacs standard.
        """
        return self.line[17:21].strip()
    
    def GetStr(self,atom_sq_number=1,res_sq_number=1,atom_index=1):
        self.UpdateCrd()
        self.line  = self.line[:6]   + '%5d' % (atom_sq_number % 100000)  + self.line[11:]
        self.line  = self.line[:22] + '%4d' % (res_sq_number % 10000)    + self.line[26:]
        self.line  = self.line[:76] + '%2s' % self.element()[-2:] + self.line[78:]
        self.index = atom_index
        return self.line
    
    def TakeToOrigin(self):
        self.pos[0] = 0.0
        self.pos[1] = 0.0
        self.pos[2] = 0.0
        
    def ExtractElement(self):
        """ Try to extract the atom element from the atom's name.
        """
        element_str = self.line[76:78].strip()
        name_strip = self.name.strip()
        if element_str !='' :
            return element_str
        if name_strip[0:2] in Atom.two_letter_elements: # two letter elements
            out = name_strip[0].upper() + name_strip[1].lower()              # return element in upper case.
            return out.strip()
        else:
            if name_strip[0] in Atom.one_letter_elements: # one letter elements
                return name_strip[0].upper()              # return element in upper case.
        return ""
                    
                
