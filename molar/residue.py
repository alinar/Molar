'''
Created on Mar 14, 2014

@author: alinar
'''
from . import atom
import numpy as np

class Residue:

    def __init__(self,name_ext=""):
        self.atoms = []
        self.name = name_ext
        
    def AddAtom(self,line_str):
        self.atoms.append(atom.Atom(line_str))
        return self.atoms[-1]
        
    def AddAtomFromCoordinates(self,coor=[0.0,0.0,0.0],name="C"):
        self.atoms.append( atom.Atom.FromCoordinates(coor,name) )
        return self.atoms[-1]
    
    def GetStr(self,atom_number=1,res_number=1,atom_index=1):
        out_str = str()
        for atom_itr in self.atoms:
            out_str = out_str + atom_itr.GetStr(atom_number,res_number,atom_index)
            atom_number = atom_number + 1
            atom_index += 1
        return out_str
    
    def Center(self):
        center = np.array([0,0,0],dtype='f')
        num_atom = 0
        for atom_itr in self.atoms:
            center = center + atom_itr.pos;
            num_atom = num_atom +1
        center = center / float(num_atom)
        return center
    
#################################
######## iterator ###############
    def __iter__(self):
        self._i = [0,0,0,0]
        return self

    def next(self):
        """returning the atoms one by one in a for loop"""
        if self._i[3] == len(self.molecules[self._i[0]].chains[self._i[1]].residues[self._i[2]].atoms):
            self._i[3] = 0
            self._i[2] = self._i[2] + 1
        if self._i[2] == len(self.molecules[self._i[0]].chains[self._i[1]].residues):
            self._i[3] = 0
            self._i[2] = 0
            self._i[1] = self._i[1] + 1
        if self._i[1] == len(self.molecules[self._i[0]].chains):
            self._i[3] = 0
            self._i[2] = 0
            self._i[1] = 0
            self._i[0] = self._i[0] + 1
        if self._i[0] == len(self.molecules):
            raise StopIteration
        out_atom=self.molecules[self._i[0]].chains[self._i[1]].residues[self._i[2]].atoms[self._i[3]]
        self._i[3] = self._i[3] + 1
        return out_atom
    #####################################
    
###################################################
###################################################
