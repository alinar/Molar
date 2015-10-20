
'''
Created on Mar 14, 2014

@author: alinar
'''
import molecule
import atom
import numpy as np

class PdbBasic:

    def __init__(self):
        self.molecules = []
        self.index     = dict()     # make index file to be used with Gromacs.
        self.include_HETATM = True  # change it to False if you do not want to include HETATM.

    def AddMolecule(self,ex_mol=None):
        if ex_mol:
            self.molecules.append(ex_mol)
        else:           
            mol = molecule.Molecule();
            self.molecules.append(mol)
    
    def ReadFile(self,file_name_str):
        pdb_file   = open(file_name_str)
        self.lines = pdb_file.readlines()
        pdb_file.close()
        termination_reached = True # to add the first molecule
        for line in self.lines:
            if line[0:4] =='ATOM' or (line[0:6]=='HETATM' and self.include_HETATM) :
                ## add next molecule
                if termination_reached:
                    self.AddMolecule()
                    self.molecules[-1].name = line[17:21]
                    termination_reached = False
                    id1 = ""
                    id2 = ""
                    res_num_1 = 0
                    res_num_2 = -1
                    res1 = str("a")
                    res2 = str("b")
                    
                id1 = line[21]
                res_num_1   = int(line[22:26])
                res1        = line[17:20]
                if id1 != id2 :
                    ## add chain
                    self.molecules[-1].AddChain(id1)
                if (res_num_1 != res_num_2) or (res1 != res2):
                    ## add residue
                    self.molecules[-1].chains[-1].AddResidue(line[17:20])
                ## add atom
                self.molecules[-1].chains[-1].residues[-1].AddAtom(line)
                id2 = line[21] 
                res_num_2   = int(line[22:26])
                res2        = line[17:20]
            elif line[0:3]=='TER' : ## prepare to add a molecule.
                termination_reached = True
    
    def PrintInfo(self,mol=0,chain=0):
        """printing information about the distributions.
        """
        mol_set = dict()
        for mol in self.molecules:
            if mol.name in mol_set:
                mol_set[mol.name] += 1
            else:
                mol_set[mol.name] = 1
        print "molecule list and quantity:\n" , mol_set
        
    def WriteOnFile(self,file_name_str):
        file_str = str()
        for molecule in self.molecules:
            file_str = file_str + molecule.GetStr()
        new_file = file(file_name_str , 'w')
        new_file.write(file_str)
        new_file.close()
    
    def WriteOnFileGMX(self,file_name_str,write_index=False,index_file="index.ndx"):
        """Write files to be used in GROMACS: 
        same type molecules are written consecutively.
        The type of a molecule is its "name" property.
        """
        mol_set = dict()
        for mol in self.molecules:
            if mol.name in mol_set:
                mol_set[mol.name] += 1
            else:
                mol_set[mol.name] = 1
        file_str = str()
        atom_index = 1
        for type in mol_set:
            print "writing molecule type: ",type
            for molecule in self.molecules:
                if type == molecule.name:
                    file_str = file_str + molecule.GetStr(atom_index)
                    atom_index += molecule.NumOfAtoms()
        new_file = file(file_name_str , 'w')
        new_file.write(file_str)
        new_file.close()
        
        ############### index #######################
        ## make System index including all the atoms.
        if write_index:
            ## System group with all of the atoms.
            self.index["System"]=list()
            for atom in self:
                self.index["System"].append(atom)
            ## 
            ndx_string = str()
            for key,value in self.index.iteritems():
                ndx_string += "[ %s ]\n" % (key)
                line_watch = 0
                for atom in value:
                    if line_watch == 10:
                        ndx_string += "\n"
                        line_watch = 0
                    ndx_string += "%9d " % (atom.index)
                    line_watch += 1
                ndx_string += "\n"
            new_file = file(index_file , 'w')
            new_file.write(ndx_string)
            new_file.close()
    
    def Center(self):
        center = np.array([0,0,0],dtype='f')
        num_mol = 0
        for atom in self:
            center = center + atom.pos
            num_mol = num_mol +1
        center = center / num_mol
        return center
    
    def NumOfResidues(self):
        num = 0
        for mol in self.molecules:
            for chain in mol.chains:
                num = num + len(chain.residues)
        return num
    
    def BoundAtoms(self):
        min_x = 1e9 
        max_x = -1e9
        min_y = 1e9
        max_y = -1e9
        min_z = 1e9
        max_z = -1e9
        for atom in self:
            if min_x > atom.pos[0]:
                atom_min_x = atom
            if min_y > atom.pos[1]:
                atom_min_y = atom
            if min_z > atom.pos[2]:
                atom_min_z = atom
            if max_x < atom.pos[0]:
                atom_max_x = atom
            if max_y < atom.pos[1]:
                atom_max_y = atom
            if max_z < atom.pos[2]:
                atom_max_z = atom
        output = [[atom_min_x,atom_max_x],[atom_min_y,atom_max_y],[atom_min_z,atom_max_z]]
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

    def CopyFrom(self,pdb_ex,res1,res2):
        """copies residues of pdb_ex from res1 until (including) res2 into the pdb"""
        if len(self.molecules) == 0:
            self.AddMolecule()
        if len(self.molecules[-1].chains) == 0:
            self.molecules[-1].AddChain("A")
        for i in range(res1,res2+1):
            self.molecules[-1].chains[-1].AddResidue((pdb_ex.molecules[-1].chains[-1].residues[i].name))
            for atom in pdb_ex.molecules[-1].chains[-1].residues[i].atoms:
                self.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)
        
    def AddFrom(self,pdb_ex,res1,res2):
        """adds the residues of pdb_ex from  res1 up to res2. it uses the pointers and does not copy 
        the residues"""
        if len(self.molecules) == 0:
            self.AddMolecule()
        if len(self.molecules[-1].chains) == 0:
            self.molecules[-1].AddChain("A")
        for i in range(res1,res2+1):
            self.molecules[-1].chains[-1].residues.append(pdb_ex.molecules[-1].chains[-1].residues[i])
            
    def Prune(self):
        """Remove empty components.
        """
        # Residues
        for mol in self.molecules:
            for ch in mol.chains:
                for re in ch.residues[:]: # with [:] a copy is made from the list.
                                          # CAUTION: do not loop over a list and change the list inside the loop.
                    if len(re.atoms)==0:
                         ch.residues.remove(re)                            
        # Chains
        for mol in self.molecules[:]:
            for ch in mol.chains[:]:
                if len(ch.residues)==0:
                    mol.chains.remove(ch)
        # Molecules
        for mol in self.molecules[:]:
            if len(mol.chains)==0:
                self.molecules.remove(mol)
                
    def GetAtomsByName(self,name_str):
        """Atom name has maximum four chars.
        """
        output = list()
        for atom in self:
            n = atom.name.replace(" ", "")
            if n == name_str:
                output.append(atom)
        return output

    def GetAtomByName(self,name_str):
        """Returns the first match
        Atom name has maximum four chars.
        """
        output = list()
        for atom in self:
            n = atom.name.replace(" ", "")
            if n == name_str:
                return atom 
        pass

    def RemoveAtomByFunc(self,func):
        """removes atoms that make func return anything else than False or None
        """
        for mol in self.molecules:
            for ch in mol.chains:
                for re in ch.residues:
                    for a in re.atoms[:]:
                        if func(a):
                            re.atoms.remove(a)

    def RemoveAtom(self,atom_in):
        for mol in self.molecules:
            for ch in mol.chains:
                for re in ch.residues:
                    for a in re.atoms[:]:
                        if a is atom_in:
                            re.atoms.remove(a)
                            return
                            
    def Concatenate(self,pdb_ext):
        """concatenates pdb_ext to slef."""
        for molecule in pdb_ext.molecules:
            self.AddMolecule(molecule)
    
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