'''
Created on Mar 14, 2014

@author: alinar
'''
from . import chain
import numpy as np
import random

class Molecule:

    def __init__(self,name_=""):
        self.chains = []
        self.name   = name_ #if needed!#
        self.bonds_established = False
        
    def AddChain(self,chain_id_str=" "):
        new_chain = chain.Chain(chain_id_str)
        self.chains.append(new_chain)
        return self.chains[-1]
        
    def GetStr(self,atom_index = 1, resid=1 ,make_TER=False,update_name=False):
        """name_res writes the molecule name on residue name domain on the pdb file.
           (residue name has 4 characters to be compatible with GMX.)
        """
        out_str = str()
        if not update_name:
            for chain_itr in self.chains:
                out_str = out_str + chain_itr.GetStr(atom_index,make_TER)
                atom_index += chain_itr.NumOfAtoms()
        else: # update_name ==True:
            for chain_itr in self.chains:
                for line in chain_itr.GetStr( atom_index=atom_index , resid=resid , make_TER=make_TER ).splitlines():
                    line = line [:17] + "%4.4s"%self.name + line[21:] + "\n"
                    out_str = out_str + line
                    atom_index += chain_itr.NumOfAtoms()
                resid += len(chain_itr.residues)
        return out_str
    
    def NumOfAtoms(self):
        output = 0
        for chain_itr in self.chains:
            output += chain_itr.NumOfAtoms()
        return output
    
    def Center(self):
        center = np.array([0,0,0],dtype='f')
        num_ch = 0
        for chain_itr in self.chains:
            center = center + chain_itr.Center();
            num_ch = num_ch + 1
        center = center / num_ch
        return center
    
    def AddAtom(self,line):
        """Adding atom by the provided string(line) to the molecule.
           Acting similar to PdbBasic.ReadFile()
        """
        q = 0
        s = -1
        if line[0:4]  == 'ATOM':
            id1 = line[21]
            q   = int(line[22:26])
            if len(self.chains)==0:
                id2 = "!"
                s   = -1
            else:
                id2 = self.chains[-1].residues[-1].atoms[-1].line[21]
                s   = int(self.chains[-1].residues[-1].atoms[-1].line[22:26])

            if id1 != id2 :
                self.AddChain(id1)
            if q != s:
                self.chains[-1].AddResidue(line[17:20])
            self.chains[-1].residues[-1].AddAtom(line)
            
    def ApplyTransform(self,trans): # trans is vtk.vtkTransform
        for atom in self:
            atom.ApplyTransform(trans)

    def GetAtomsByName(self,name_str):
        """Atom name has maximum four chars.
        """
        output = list()
        for atom in self:
            n = atom.name.strip()
            if n == name_str:
                output.append(atom)
        return output
            
    def Bounds(self):
        min_x = 1e9 
        max_x = -1e9
        min_y = 1e9
        max_y = -1e9
        min_z = 1e9
        max_z = -1e9
        for atom in self:
            if min_x > atom.pos[0]:
                min_x = atom.pos[0]
            if min_y > atom.pos[1]:
                min_y = atom.pos[1]
            if min_z > atom.pos[2]:
                min_z = atom.pos[2]
            if max_x < atom.pos[0]:
                max_x = atom.pos[0]
            if max_y < atom.pos[1]:
                max_y = atom.pos[1]
            if max_z < atom.pos[2]:
                max_z = atom.pos[2]
        output = np.array([[min_x,max_x],[min_y,max_y],[min_z,max_z]])
        return output
            
    def EstablishBonds(self,max_distance = 1.9):
        """
        max_distance is the maximum distance to be accepted as a bond. 
        """
        if self.bonds_established:
            return
        else:
            self.bonds_established = True
        atom_list = []
        for atom in self:  # atom_list should've been made since the iterator does not work for nested for loops.
            atom_list.append(atom)
            ## also reset any existing bond links.
            atom.bonded_atoms = []
        for atom1 in atom_list:
            if atom1.ExtractElement().strip() == "H": # not calculated if the first atom is hydrogen
                continue
            for atom2 in atom_list:
                if atom1 is not atom2:
                    d = Distance(atom1,atom2)
                    if d < max_distance:
                        if atom2.ExtractElement().strip() == "H":
                            if len(atom2.bonded_atoms) > 0: #if atom2 already has a bond link. also hasattr(atom2, 'distance_aux'):
                                atom2.distance_aux = min (atom2.distance_aux , d)
                                if atom2.distance_aux == d:
                                    # remove the previous link between atom2 and another atom. 
                                    atom2.bonded_atoms[0].bonded_atoms.remove(atom2)
                                    # make a new link between atom1 and atom2.
                                    atom2.bonded_atoms=[ atom1 ]
                                    atom1.bonded_atoms.append(atom2)
                            else:
                                atom2.distance_aux = d
                                atom2.bonded_atoms=[ atom1 ]
                                atom1.bonded_atoms.append(atom2)
                        else:
                            atom1.bonded_atoms.append(atom2)
                            # atom2 will find and append atom1 automatically since none are hydrogens
    
    def AddHydrogen(self):
        """
        works if atoms have num_hydrogen set. The positions of hydrogens are NOT accurate.
        """
        # add hydrogen:
        if True:
            atom_list = []
            for atom in self:  # atom_list should've been made since the iterator does not work for nested for loops.
                atom_list.append(atom)
            self.MakeBackwardLinks()
            for atom in atom_list:  # atom_list should be used instead of "self" because atoms are added.
                 for i in range(atom.num_hydrogen ):
                     d = np.array([0.0,0.0,1.1])
                     d = RotateX(d,360*random.random())
                     d = RotateY(d,360*random.random())
                     d = RotateZ(d,360*random.random())
                     H_atom = atom.residue.AddAtomFromCoordinates(atom.pos+d,"H",next_to_atom=atom)
                     H_atom.bonded_atoms = [atom]
                     atom.bonded_atoms.append(H_atom)
                 pass
            
                    
    def MakeBackwardLinks(self):
        """ Every component will know its parent component.
        """
        for chain in self.chains:
            chain.molecule = self
            for residue in chain.residues:
                residue.chain = chain
                for atom in residue.atoms:
                    atom.residue = residue
        pass        


    #################################
    ######## iterator ###############
    def __iter__(self):
        self._i = [0,0,0,0] # (_i[0] is useless)
        return self

    def next(self):
        """returning the atoms one by one in a for loop"""
        if self._i[3] == len(self.chains[self._i[1]].residues[self._i[2]].atoms):
            self._i[3] = 0
            self._i[2] = self._i[2] + 1
        if self._i[2] == len(self.chains[self._i[1]].residues):
            self._i[3] = 0
            self._i[2] = 0
            self._i[1] = self._i[1] + 1
        if self._i[1] == len(self.chains):
            raise StopIteration
        out_atom=self.chains[self._i[1]].residues[self._i[2]].atoms[self._i[3]]
        self._i[3] = self._i[3] + 1
        return out_atom
    #####################################
    
###################################################
###################################################

def Distance(atom1,atom2):
    return np.linalg.norm(  [ atom1.pos[0] - atom2.pos[0] , atom1.pos[1] - atom2.pos[1] , atom1.pos[2] - atom2.pos[2] ] )
    pass
    
def RotateX(np_array, angle):
    theta = np.deg2rad(angle)
    rotation_matrix = np.array ( [  [1. ,0. , 0.] , [ 0. , np.cos(theta) ,np.sin(theta) ] , [ 0. , -np.sin(theta) ,  np.cos(theta)]  ] ).transpose()
    return np.matmul (rotation_matrix , np_array  )
    
def RotateY(np_array, angle):
    theta = np.deg2rad(angle)
    rotation_matrix = np.array ( [  [np.cos(theta)  ,0. , np.sin(theta)] , [ 0. , 1., 0. ] , [ -np.sin(theta) , 0 ,  np.cos(theta)]  ] ).transpose()
    return(np.matmul (rotation_matrix , np_array  ) )

def RotateZ(np_array, angle):
    theta = np.deg2rad(angle)
    rotation_matrix = np.array ( [  [np.cos(theta)  , -np.sin(theta) , 0.] , [ np.sin(theta), np.cos(theta) , 0. ] , [ 0. , 0.  ,  1.  ] ] ).transpose()
    return(np.matmul (rotation_matrix , np_array  ) )
    