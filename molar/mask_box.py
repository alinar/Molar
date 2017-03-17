import random as rand
import numpy as np
import vtk
import os
import sys
import molar,molar.pdb
import vtk

PATH     =  os.environ['MOLECULES']

class Box():
    def __init__(self,x=753,y=753,z=300):
        #self.vertices=np.array(np.empty([8,3],dtype=float))
        x=0.5*x
        y=0.5*y
        z=0.5*z
        self.vertices=np.array([\
                                            [-x,-y,-z],\
                                            [x,-y,-z],\
                                            [x,y,-z],\
                                            [-x,y,-z],\
                                            [-x,-y,z],\
                                            [x,-y,z],\
                                            [x,y,z],\
                                            [-x,y,z],\
                                            ])
        # self.polygon_faces[i], ([0]-[1]) X ([2]-[1]) 's direction always points towards outside of the box.
        self.polygon_faces=[\
                                        [0,1,2],\
                                        [0,2,3],\
                                        [1,5,2],\
                                        [5,6,2],\
                                        [2,6,3],\
                                        [6,7,3],\
                                        [3,4,0],\
                                        [7,4,3],\
                                        [4,5,0],\
                                        [5,1,0],\
                                        [6,5,4],\
                                        [4,7,6],\
                                        ]
        pass
    
    def WriteTclFile(self,file_name="cuboid_box.tcl"):
        ### make the string
        str = ""
        for triangle in self.polygon_faces:
            tr_str="draw triangle "
            for vertex_indx in triangle: # each vertex has three values for its coordinates.
                tr_str += " {%.3f %.3f %.3f} " % tuple(self.vertices[vertex_indx])
            tr_str += "\n"
            str+=tr_str
        new_file = file(file_name , 'w')
        new_file.write(str)
        new_file.close()
        print "file %s was written." % file_name
        
    def FixBoxRotation(self,rotation_angles=[0.0,0.0,0.0]):
        trans=vtk.vtkTransform()
        trans.RotateZ(rotation_angles[0]) 
        trans.RotateX(rotation_angles[1])
        trans.RotateZ(rotation_angles[2])
        for v in self.vertices:
            trans.TransformPoint(v,v) # in-place transformation is tested in this case and works fine. 
        pass
        
    def MaskPdb(self,pdb_in_file,out_pdb_filename="masked_pdb.pdb",specimen_rotation_angles=False): 
        """ The box is stationary and the specimen rotates inside. the specimen's center will be aligned to the center of the box which is the origin.
        """
        inward_vectors = []
        for triangle in self.polygon_faces:
            v1 = self.vertices[triangle[1]] - self.vertices[triangle[0]] 
            v2 = self.vertices[triangle[1]] - self.vertices[triangle[2]]
            inward_vectors.append(np.cross(v2,v1))
        triangle_vector = zip(self.polygon_faces,inward_vectors)
        #print triangle_vector
        new_pdb_str = ""
        pdb_file   = open(pdb_in_file)
        lines = pdb_file.readlines()
        pdb_file.close()
        pos = np.array([0.0,0.0,0.0],dtype='f')
        t=len(lines)
        ## find the center
        n=0
        center = np.array([0.0,0.0,0.0],dtype='f')
        print "finding the center ..."
        for line in lines:
            molar.pdb.update_progress(float(n)/t)
            if line[0:4] =='ATOM':
                ## extract pos
                pos[0] = float(line[30:38])
                pos[1] = float(line[38:46])
                pos[2] = float(line[46:54])
                ##
                center+=pos
                n+=1
        center /= n
        print "center is at: ", center
        ## Transformation  ##
        trans=vtk.vtkTransform()
        trans.PostMultiply()
        trans.Identity()
        trans.Translate(-1*center)
        if specimen_rotation_angles:
            trans.RotateZ(specimen_rotation_angles[0]) 
            trans.RotateX(specimen_rotation_angles[1])
            trans.RotateZ(specimen_rotation_angles[2])
        ##
        i=0
        print "\nRotating and masking ..."
        for line in lines:
            i+=1
            molar.pdb.update_progress(float(i)/t)
            if line[0:4] =='ATOM':
                ## extract pos
                pos[0] = float(line[30:38])
                pos[1] = float(line[38:46])
                pos[2] = float(line[46:54])
                ##
                trans.TransformPoint(pos,pos)
                atom_is_in = True
                for triangle , vector in triangle_vector:
                    value = np.dot( ( pos - self.vertices[triangle[1]]) , vector )
                    if value < 0.0:
                        atom_is_in = False
                        break
                if atom_is_in:
                    ## write to new_pdb_str
                    line = line[:30] + "%8.3f" % pos[0] + line[38:]
                    line = line[:38] + "%8.3f" % pos[1] + line[46:]
                    line = line[:46] + "%8.3f" % pos[2] + line[54:]
                    new_pdb_str+=line
            else:
                ##NOT ATOM
                #new_pdb_str+=line
                pass
        
        new_file = file(out_pdb_filename , 'w')
        new_file.write(new_pdb_str)
        new_file.close()
        print "file %s was written." % out_pdb_filename

        ####
       
############################################
############################################
 
def Make_test(dim=[400,400,400],distance = 20,out_file="test.pdb"):
    test_atom=molar.pdb.Pdb(os.path.join(PATH , "test.pdb"))
    print test_atom.molecules[0].chains[0].residues[0].atoms[0].pos
    #test_atom.Show("sphere")
    out = molar.pdb.Pdb()
    trans=vtk.vtkTransform()
    for x in range (0,dim[0],distance):
        for y in range (0,dim[1],distance):
            for z in range (0,dim[2],distance):
                trans.Identity()
                trans.Translate([x,y,z])
                out.CatTransformed(test_atom,trans)
    out.WriteOnFile(out_file)
        