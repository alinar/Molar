'''
Created on Mar 14, 2014

@author: alinar
'''
from . import molecule
from . import atom
import numpy as np
import sys
from collections import OrderedDict

class PdbBasic:

    def __init__(self):
        self.molecules = []
        self.index     = OrderedDict()     # make index file to be used with Gromacs. OrderedDict is like Dict but it keeps the order.
        self.include_HETATM = True  # change it to False if you do not want to include HETATM.
        self.mol_set = False
        self.empty   = True         # some methods use it.
        self.cryst   = False
        self.cryst_string = False
 
    def Update(self):
        """to be overridden by subclasses if needed. This method is called before any file writing.
        """
        pass
 
    def AddMolecule(self,ex_mol=None,name=""):
        if ex_mol:
            self.molecules.append(ex_mol)
        else:           
            mol = molecule.Molecule(name);
            self.molecules.append(mol)
        return self.molecules[-1]
    
    def ReadFile(self,file_name_str):
        pdb_file   = open(file_name_str)
        self.lines = pdb_file.readlines()
        pdb_file.close()
        termination_reached = True # to add the first molecule
        for line in self.lines:
            if line[0:4] =='ATOM' or (line[0:6]=='HETATM' and self.include_HETATM) :
                ## add next molecule
                if termination_reached:
                    self.AddMolecule( name = line[17:21].strip() )
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
                if (res_num_1 != res_num_2) or (res1 != res2) or len(self.molecules[-1].chains[-1].residues)==0:
                    ## add residue
                    self.molecules[-1].chains[-1].AddResidue(line[17:20])
                ## add atom
                self.molecules[-1].chains[-1].residues[-1].AddAtom(line)
                id2 = id1
                res_num_2   = int(line[22:26])
                res2        = line[17:20]
            elif line[0:3]=='TER' : ## prepare to add a molecule.
                termination_reached = True
            elif line[0:6]=='CRYST1': ## read crystal line
                self.cryst_string = line
                
    def ReadFileGMX(self,file_name_str):
        """add new molecule for every residue names. 
           ignoring chains.
           residue name field has one more character on the right side.
           ([17:21] instead of [17:20])
        """
        pdb_file   = open(file_name_str)
        self.lines = pdb_file.readlines()
        pdb_file.close()
        res_num_1 = 0
        res_num_2 = -1
        res1 = str("a")
        res2 = str("b")
        termination_reached = True # to add the first molecule
        for line in self.lines:
            if line[0:4] =='ATOM' or (line[0:6]=='HETATM' and self.include_HETATM) :
                res_num_1   = int(line[22:26])
                res1        = line[17:21]
                ## add next molecule
                if termination_reached or (res_num_1 != res_num_2) or (res1 != res2):
                    ##
                    self.AddMolecule( name=line[17:21].strip() ) # not allow any space
                    self.molecules[-1].AddChain(" ")
                    self.molecules[-1].chains[-1].AddResidue(line[17:21])
                    ##
                    res_num_1   = int(line[22:26])
                    res1        = line[17:21]
                    termination_reached = False
                ## add atom
                self.molecules[-1].chains[-1].residues[-1].AddAtom(line)
                res_num_2   = int(line[22:26])
                res2        = line[17:21]
            elif line[0:3]=='TER' : ## prepare to add a molecule.
                termination_reached = True
            elif line[0:6]=="CRYST1" :  ## read the unit-cell's dimensions
                cx = float(line[6:15])
                cy = float(line[15:24])
                cz = float(line[24:33])
                self.cryst=[cx,cy,cz]

    def Concatenate(self,pdb_ext):
        """Transforms a hard copy of pdb_ext and concatenates it to self."""
        for molecule in pdb_ext.molecules:
            self.AddMolecule(name=molecule.name)
            for chain in molecule.chains:
                self.molecules[-1].AddChain(chain.id)
                for residue in chain.residues:
                    self.molecules[-1].chains[-1].AddResidue(residue.name)
                    for atom in residue.atoms:
                        atom.UpdateCrd()
                        self.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)

    def MakeMolList(self):
        self.mol_set = OrderedDict() # instead of dict() to preserve the order 
        for mol in self.molecules:
            if mol.name in self.mol_set:
                self.mol_set[mol.name] += 1
            else:
                self.mol_set[mol.name] = 1
        
    def PrintInfo(self,mol=0,chain=0):
        """printing information about the distributions.
        """
        self.MakeMolList()
        n = self.NumOfResidues()
        for mol_name in self.mol_set:
            print "Molecule Type: %s,   Molecule Quantity: %d ,  ratio: %.1f%%" % (mol_name , self.mol_set[mol_name],  100*float(self.mol_set[mol_name])/n )
        #print "molecule list and quantity:\n" , self.mol_set
        
    def DefCryst(self,a=100,b=100,c=100,alpha=90,beta=90,gamma=90,sGroup="P 1", zvalue=1):
        """This method is to be called by the user.
        """
        self.cryst_string = "CRYST1% 9.3f% 9.3f% 9.3f% 7.2f% 7.2f% 7.2f %-11s%4d\n" % (a,b,c,alpha,beta,gamma,sGroup,zvalue)
        
    def DefBox(self , d = [0.0,0.0,0.0]):
        """define periodic box around the system added the distance d.
        """
        b = self.Bounds()
        l0=b[0,1]-b[0,0]
        l1=b[1,1]-b[1,0]
        l2=b[2,1]-b[2,0]
        self.DefCryst ( a=l0 + 2*d[0] , b=l1 + 2*d[1] ,c=l2 + 2*d[2] )
        
    def WriteOnFile(self,file_name_str,make_TER=False,include_CONECT=False):
        self.Update()
        file_str = str()
        atom_index = 1
        if self.cryst_string:
            file_str = file_str + self.cryst_string
        for molecule in self.molecules:
            file_str = file_str + molecule.GetStr(atom_index = 1,make_TER=False,update_name=False)
            atom_index         += molecule.NumOfAtoms()
        if include_CONECT:
            file_str+=self.GetConectStr()
        new_file = file(file_name_str , 'w')
        print "Writing to the file: ",file_name_str
        new_file.write(file_str)
        new_file.close()
    
    def WriteOnFileGMX(self,file_name_str,make_TER=False,include_CONECT=False):
        """Write files to be used in GROMACS: 
        same-type molecules are written consecutively.
        The type of a molecule is its "name" property.
        """
        self.Update()
        self.MakeMolList()
        file_str = str()
        if self.cryst_string:
            file_str = file_str + self.cryst_string
        atom_index_ = 1
        resid_ = 1
        for type in self.mol_set:
            for molecule in self.molecules:
                if molecule.name == type:
                    file_str = file_str + molecule.GetStr(atom_index=atom_index_ , resid=resid_, make_TER=make_TER , update_name=True )
                    atom_index_         += molecule.NumOfAtoms()
                    resid_ +=1
            print "writing molecule type: ",type
        
        if include_CONECT:
            file_str+=self.GetConectStr()    
        new_file = file(file_name_str , 'w')
        new_file.write(file_str)
        new_file.close()
        
    def GetConectStr(self):
        """ Returning the CONECT string. It should be called after or inside WriteOnFile or WriteOnFileGMX
        because the atom.index should be established first.
        """
        for mol in self.molecules:
            mol.EstablishBonds()
        out_str = ""
        for atom in self:
            if len(atom.bonded_atoms) ==1:
                line="CONECT%5d%5d" % (atom.index,atom.bonded_atoms[0].index)
            elif len(atom.bonded_atoms) ==2:
                line="CONECT%5d%5d%5d" % (atom.index,atom.bonded_atoms[0].index , atom.bonded_atoms[1].index)
            elif len(atom.bonded_atoms) ==3:
                line="CONECT%5d%5d%5d%5d" % (atom.index,atom.bonded_atoms[0].index , atom.bonded_atoms[1].index ,  atom.bonded_atoms[2].index)
            elif len(atom.bonded_atoms) ==4:
                line="CONECT%5d%5d%5d%5d%5d" % (atom.index,atom.bonded_atoms[0].index , atom.bonded_atoms[1].index ,  atom.bonded_atoms[2].index , atom.bonded_atoms[3].index)
            out_str+=line
            out_str+="\n"
        return out_str
        
    def WriteIndexGMX(self,index_file="index.ndx",include_system=False):
        ############### index ###################
        ## make System index including all the atoms.
        """Must be called after WriteOnFileGMX()
        """
        self.Update()
        if include_system:
            ## System group with all of the atoms. Other index groups should be made.
            self.index["System"]=list()
            for atom in self:
                self.index["System"].append(atom)
            ##
        if True: 
            ndx_string = str()
            for key,value in self.index.iteritems():
                ndx_string += "[ %s ]\n" % (key)
                line_watch = 0
                for atom in value:
                    #print key ,  atom.name , atom.index , atom.residue.chain.molecule.name
                    if line_watch == 10:
                        ndx_string += "\n"
                        line_watch = 0
                    ndx_string += "%10d " % (atom.index)
                    line_watch += 1
                ndx_string += "\n"
            new_file = file(index_file , 'w')
            new_file.write(ndx_string)
            new_file.close()
        print "Index file is written to : ",index_file
            
    def WriteTopologyGMX(self,topol_file="topol_list.top",name_map_dic=False):
        ############### topology ###################
        self.Update()
        if not self.mol_set:
            self.MakeMolList()
        if True:
            top_string=str()
            for type,value in self.mol_set.iteritems():
                if (not name_map_dic) or (type not in name_map_dic) :
                    top_string += "%-15s%d\n" % (type,value)
                else:
                    top_string += "%-15s%d\n" % (name_map_dic[type],value)
            new_file = file(topol_file , 'w')
            new_file.write(top_string)
            new_file.close()
        print "Topology file is written to : ",topol_file

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
    
    def NumOfAtoms(self):
        i=0
        for atom in self:
            i+=1
        return i
    
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
                                          # CAUTION: do not loop over a list and remove 
                                          #          elements from the same list inside the loop.
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
                
    def GetMoleculesByName(self,name_str):
        """Atom name has maximum four chars.
        """
        output = list()
        for mol in self.molecules:
            n = mol.name.strip()
            if n == name_str:
                output.append(mol)
        return output
                
    def GetAtomsByName(self,name_str):
        """Atom name has maximum four chars.
        """
        output = list()
        for atom in self:
            n = atom.name.strip()
            if n == name_str:
                output.append(atom)
        return output

    def GetAtomByName(self,name_str):
        """Returns the first match
        Atom name has maximum four chars.
        """
        output = list()
        for atom in self:
            n = atom.name.strip()
            if n == name_str:
                return atom 
        return False

    def RemoveAtomByFunc(self,func):
        """removes atoms that make func return anything other than False or None
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
            self.AddMolecule(molecule,name=molecule.name)
    
    def MakeBackwardLinks(self):
        """ Every component will know its parent component.
        """
        for molecule in self.molecules:
            for chain in molecule.chains:
                chain.molecule = molecule
                for residue in chain.residues:
                    residue.chain = chain
                    for atom in residue.atoms:
                        atom.residue = residue
        pass
    
    def BoxMask(self,box_dimensions,verbose=False):
        """Mask self into the cuboid with the box_dimensions.
        self should be centered first [self.BringToCenter()]
        IMPORTANT: If the structure is too large, consider using molar.mask_box instead.
        """
        self.Update()
        half_dim = [ 0.5*box_dimensions[0] , 0.5*box_dimensions[1] , 0.5*box_dimensions[2] ]
        pdb_out = PdbBasic()
        i=0
        if verbose:
            N=self.NumOfAtoms()
        for molecule in self.molecules:
            pdb_out.AddMolecule(name=molecule.name)
            no_atom_in_mol = True
            for chain in molecule.chains:
                pdb_out.molecules[-1].AddChain(chain.id)
                for res in chain.residues:
                    pdb_out.molecules[-1].chains[-1].AddResidue(res.name)
                    for atom in res.atoms:
                        atom.UpdateCrd()
                        i+=1
                        if abs(atom.pos[0]) < half_dim[0] and abs(atom.pos[1]) < half_dim[1] and abs(atom.pos[2]) < half_dim[2]:
                            pdb_out.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)
                            no_atom_in_mol = False
                        if verbose:
                            update_progress(float(i)/N)
            if no_atom_in_mol:
                # remove the last molecule if no atom is added. to avoid empty molecule.
                pdb_out.molecules.pop(-1)
        ###
        del(self.molecules) # to make sure although it is garbage-collected in the next line
        self.molecules = pdb_out.molecules
    
    #################################
    ######## iterator ##################
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

def update_progress(progress):
    """update_progress() : Displays or updates a console progress bar
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%
    """
    barLength = 30 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rProgress: [%s%s] %3.1f %% %s"  %  ("|"*block , "-"*(barLength-block) , progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
