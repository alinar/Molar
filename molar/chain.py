'''
Created on Mar 14, 2014

@author: alinar
'''
from . import residue
import numpy as np

class Chain:
    
    def __init__(self,id_ext):
        self.id = id_ext
        self.residues = []
        
    def AddResidue(self,name_str):
        self.residues.append(residue.Residue(name_str))
        
    def GetStr(self,atom_index=1,resid=1,make_TER=True):
        out_str = str()
        atom_number = 1
        res_number  = resid

        for res_itr in self.residues:
            out_str = out_str + res_itr.GetStr(atom_number,res_number,atom_index)
            atom_number = atom_number + len(res_itr.atoms)
            res_number  = res_number  + 1
            atom_index += len(res_itr.atoms)
        # make TER
        if make_TER:
            ter = 'TER' + ' ' * 24 + '\n'
        else:
            ter = ""
        out_str = out_str + ter
        return out_str
    
    def NumOfAtoms(self):
        output = 0
        for res_itr in self.residues:
            output += len(res_itr.atoms)
        return output
    
    def Center(self):
        center = np.array([0,0,0],dtype='f')
        num_res = 0
        for res_itr in self.residues:
            center = center + res_itr.Center();
            num_res = num_res +1
        center = center / num_res
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
