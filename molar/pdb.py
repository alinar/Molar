import pdb_basic
import atom
import pdb_viewer
import vtk
import random as rand
import numpy as np
import math
import copy
import time, sys
import itertools

class Pdb(pdb_basic.PdbBasic):
    """Stores the pdb_basic file for making transformations."""
    def __init__(self,input_file=False):
        pdb_basic.PdbBasic.__init__(self)
        self.hinge=[0,1]
        if input_file:
            self.ReadFile(input_file)
    
    def ApplyTransform(self,trans):
        """Transforms the whole pdb_basic"""
        for molecule in self.molecules:
            for chain in molecule.chains:
                for residue in chain.residues:
                    for atom in residue.atoms:
                        atom.ApplyTransform(trans)
                        
    def CatTransformed(self,pdb_ext,trans):
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
                        self_atom = self.molecules[-1].chains[-1].residues[-1].atoms[-1]
                        self_atom.ApplyTransform(trans)
                        
    def SetHinge(self,res_num1=0,res_num2=1):
        """Sets the hinge to be between two residues"""
        self.hinge = [res_num1,res_num2]
    
    def HingeTrans(self,trans,mode_str='up'):
        """transforms only one of the parts separated from the hinge."""
        translation=vtk.vtkTransform()
        if mode_str == 'down':
            pos=self.molecules[-1].chains[-1].residues[self.hinge[1]].atoms[0].pos.copy()
            translation.Translate(-pos[0], -pos[1], -pos[2])
            self.ApplyTransform(translation)
            for residue in self.molecules[-1].chains[-1].residues[:self.hinge[1]]:
                for atom in residue.atoms:
                    atom.ApplyTransform(trans)
            
        elif mode_str == 'up':
            pos=self.molecules[-1].chains[-1].residues[self.hinge[0]].atoms[-1].pos.copy()
            translation.Translate(-pos[0], -pos[1], -pos[2])
            self.ApplyTransform(translation)
            for residue in self.molecules[-1].chains[-1].residues[self.hinge[1]:]:
                for atom in residue.atoms:
                    atom.ApplyTransform(trans)
        else:
            print ("Enter a valid mode.")
            return
        translation.Identity()
        translation.Translate(pos[0], pos[1], pos[2])
        self.ApplyTransform(translation)

    def BreakHinge(self,mode_str='up'):
        """Breaks the pdb_basic form the hinge and returns one of the parts
        depending on the mode_str"""
        pdb_out = Pdb()
        pdb_out.AddMolecule()
        pdb_out.molecules[-1].AddChain('A')
        if mode_str == 'up':
            for residue in self.molecules[-1].chains[-1].residues[self.hinge[1]:]:
                pdb_out.molecules[-1].chains[-1].AddResidue(residue.name)
                for atom in residue.atoms:
                    pdb_out.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)
        elif mode_str == 'down':
            for residue in self.molecules[-1].chains[-1].residues[:self.hinge[1]]:
                pdb_out.molecules[-1].chains[-1].AddResidue(residue.name)
                for atom in residue.atoms:
                    pdb_out.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)
        return pdb_out
    
    def MakeCopy(self):
        """makes a hard copy and returns it."""
        self.Update()
        out = Pdb()
        for molecule in self.molecules:
            out.AddMolecule(name=molecule.name)
            for chain in molecule.chains:
                out.molecules[-1].AddChain(chain.id)
                for residue in chain.residues:
                    out.molecules[-1].chains[-1].AddResidue(residue.name)
                    for atom in residue.atoms:
                        atom.UpdateCrd()
                        out.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)
        return out
        
    def Clone(self):
        self.Update()
        pdb_out = Pdb()
        for molecule in self.molecules:
            pdb_out.AddMolecule(name=molecule.name)
            for chain in molecule.chains:
                pdb_out.molecules[-1].AddChain(chain.id)
                for res in chain.residues:
                    pdb_out.molecules[-1].chains[-1].AddResidue(res.name)
                    for atom in res.atoms:
                        atom.UpdateCrd()
                        pdb_out.molecules[-1].chains[-1].residues[-1].AddAtom(atom.line)
        return pdb_out
    
    def Show(self,mode="dot",ext_actors=None,axis=True,only_ext=False):
        """
           Visualize the system.
           ext_actors is a list of external vtkActor(s).
        """
        self.Update()
        viewer = pdb_viewer.PdbViewer(ext_actors) 
        viewer.SetPdb(self)
        if axis ==True:
            viewer.AxesOn()
        viewer.Show(mode,only_ext)

    def RetrunRenderer(self,mode="dot",background=[0.1,0.2,0.3]):
        self.Update()
        viewer = pdb_viewer.PdbViewer()
        viewer.SetPdb(self)
        renderer = viewer.RetrunRenderer(mode,False,background)
        return renderer
    
    def PointReflect(self):
        for atom in self:
            atom.SetPosition(-1 * atom.pos)
    
    def BringToCenter(self):
        c = self.Center()
        trans = vtk.vtkTransform()
        trans.PostMultiply()
        trans.Translate(-1 * c)
        self.ApplyTransform(trans)
        
    def BringToPositiveY(self):
        self.BringToCenter()
        trans = vtk.vtkTransform()
        trans.Translate(0.0,-self.Bounds()[1,0],0.0)
        self.ApplyTransform(trans)
        
    def BringToNegativeY(self):
        self.BringToCenter()
        trans = vtk.vtkTransform()
        trans.Translate(0.0,-self.Bounds()[1,1],0.0)
        self.ApplyTransform(trans)

    def BringToPositiveCorner(self):
        bound = self.Bounds()
        trans = vtk.vtkTransform()
        trans.Translate([-bound[0,0],-bound[1,0],-bound[2,0]])
        self.ApplyTransform(trans)

    def SelectBond(self,atom1,atom2):
        self.bond_atoms = [atom1,atom2]
                    
    def RotateBond(self,angle,mode="up"):
        for mol in self.molecules:
            mol.EstablishBonds()
        trans = vtk.vtkTransform()
        trans.PostMultiply()
        pos = self.bond_atoms[1].pos.copy()
        v   = self.bond_atoms[0].pos - self.bond_atoms[1].pos
        self.Translate(-1*pos)
        trans.RotateWXYZ(angle, v[0], v[1], v[2])
        ## nested recursive functions ##
        def Recursive_trans(in_atom,in_trans,in_prev_atoms_list):
            """
            """
            in_atom.ApplyTransform(in_trans)
            in_prev_atoms_list.append(in_atom)
            for bonded_atom in in_atom.bonded_atoms:
                if bonded_atom not in in_prev_atoms_list:
                    Recursive_trans(bonded_atom,in_trans,in_prev_atoms_list)
            return
        ###
        prev_atoms_list_ = [ self.bond_atoms[0] ]
        Recursive_trans( self.bond_atoms[1] , trans , prev_atoms_list_)
        
        self.Translate(pos)
                    
    def RandomizeLocation(self,vec):
        """Randomiztion of location respecting to the vec. 
        """
        trans = vtk.vtkTransform() 
        random_vector = [vec[0]*rand.random(),vec[1]*rand.random(),vec[2]*rand.random()]
        trans.Translate(random_vector)
        self.ApplyTransform(trans)
        
    def Translate(self,vec):
        trans = vtk.vtkTransform()
        trans.Translate(vec)
        self.ApplyTransform(trans)
        
    def RotateZ(self,angle):
        trans = vtk.vtkTransform()
        trans.RotateZ(angle)
        self.ApplyTransform(trans)
        
    def RotateX(self,angle):
        trans = vtk.vtkTransform()
        trans.RotateX(angle)
        self.ApplyTransform(trans)
    
    def RotateY(self,angle):
        trans = vtk.vtkTransform()
        trans.RotateY(angle)
        self.ApplyTransform(trans)
        
    def RotateWXYZ(self,angle,vec):
        trans = vtk.vtkTransform()
        trans.RotateWXYZ(angle,vec)
        self.ApplyTransform(trans)
        
    def Flip(self):
        trans = vtk.vtkTransform()
        trans.RotateZ(180)
        self.ApplyTransform(trans)
        
    def CatTransformedCautious(self,pdb_ext,trans, pointlocator, cutoff = 1.2):
        """Transforms pdb_ext and concatenates it to self. checks if the new molecule is not within the
        distance of cutoff from all the points in pointlocator.dataset .
        pointlocator :  vtkPointLocator()
        """
        pdb_aux = pdb_ext.MakeCopy()
        pdb_aux.ApplyTransform(trans)
        dataset = pointlocator.GetDataSet()
        points = dataset.GetPoints()
        distance = vtk.mutable(0.0) # FindClosestPointWithinRadius writes on distance and its useless here.
        if not points:
            points = vtk.vtkPoints()
            dataset.SetPoints(points)
        else:
            pointlocator.Update() # not sure what it does! but points should be defined before!
        if pointlocator.GetDataSet().GetNumberOfPoints()!=0:      # do not check if empty
            for atom_itr in pdb_aux:
                #print pointlocator.GetDataSet().GetNumberOfPoints()
                pointid = pointlocator.FindClosestPointWithinRadius(cutoff , atom_itr.pos, distance)
                if pointid != -1: # There have been an atom within "cutoff" distance.
                    return False
        ## There have been no atom within the cutoff distance. 
        self.Concatenate(pdb_aux)
        ## add the newly added atom's position and its repeats to the pointlocator.
        pos = [0,0,0]
        for atom_itr in pdb_aux:
            points.InsertNextPoint(atom_itr.pos)
            pointlocator.InsertNextPoint(atom_itr.pos)
            ## if atom position is close to any border of cell, it should be repeated in the locator.
            close_coords = []  
            for coord in [0,1,2]:
                if (self.unit_cell_size - atom_itr.pos[coord]) < (2*cutoff)  or (atom_itr.pos[coord] > self.unit_cell_size):
                    close_coords.append(-1*coord)
                elif abs(atom_itr.pos[coord]) < (2*cutoff) or (atom_itr.pos[coord] < 0) :
                    close_coords.append(coord)

            for i in [1,2,3]: # extract subsets of close_coords with 1, 2 or 3 members
                close_coords_set_of_subsets = set( itertools.combinations(close_coords,i) )
                for close_coords_subset in close_coords_set_of_subsets:
                    atom_pos_aux = np.copy(atom_itr.pos)
                    for close_coords_member in close_coords_subset:
                        atom_pos_aux[abs(close_coords_member)] += np.sign(close_coords_member) * self.unit_cell_size
                    points.InsertNextPoint(atom_pos_aux)
                    pointlocator.InsertNextPoint(atom_pos_aux)
        return True
         
###################################################
###################################################

def MergePdb(pdb_list):
    out = Pdb()
    for p in pdb_list:
        for m in p.molecules:
            out.AddMolecule(m)
    return out
    
def RotationToY(vec_ext):
    """Returns a transform that rotates vec to be parallel to Y axis."""
    vec = vec_ext / np.linalg.norm(vec_ext)
    trans = vtk.vtkTransform()
    trans.PostMultiply()
    phi = np.degrees(math.acos(vec[1]))
    if np.abs(vec[0]) > 0.000001:
        theta = math.degrees(math.atan(vec[2]/vec[0]))
    else:
        if np.abs(vec[2]) > 0.000001:
            theta = np.degrees(np.arctan(vec[2] * np.inf))
        else:
            theta = 0.0
    if vec[0] < 0:
        theta = theta + 180
    trans.RotateY(theta)
    trans.RotateZ(phi)
    return trans

def RotationToZ(vec_ext):
    """Returns a transform that rotates vec to be parallel to Z axis."""
    vec = vec_ext / np.linalg.norm(vec_ext)
    trans = vtk.vtkTransform()
    trans.PostMultiply()
    phi = math.degrees(math.acos(vec[2]))
    
    if np.abs(vec[1]) > 0.000001:
        theta = math.degrees(math.atan(vec[0]/vec[1]))
    else:
        if np.abs(vec[0]) > 0.000001:
            theta = np.degrees(np.arctan(vec[0] * np.inf))
        else:
            theta = 0.0
    if vec[1] < 0:
        theta = theta + 180
    trans.RotateZ(theta)
    trans.RotateX(phi)
    return trans

def RotateToParallel(v_const,v_2):
    """ returns a transform that rotates v_2 to be parallel to v_const.
    """
    trans1 = RotationToY(v_const)
    trans2 = RotationToY(v_2)
    trans1.Inverse()
    trans2.PostMultiply()
    trans2.Concatenate(trans1)
    return trans2

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

def FindClosestPoint(point_in,point_set):
    """returns the ID if the closest point to point_in in point_set. (rigorous algorithm)
    """ 
    min_distance = 1e8
    id           = -1
    for i,p in enumerate(point_set):
        d = (point_in[0] - p[0])**2 + (point_in[1] - p[1])**2 + (point_in[2] - p[2])**2
        if d<min_distance:
            min_distance = d
            id = i
    return id