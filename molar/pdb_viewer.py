from . import pdb_basic
import vtk
import math

class PdbInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self,parent=None):
        #self.AddObserver("LeftButtonPressEvent",self.middleButtonPressEvent)
        #self.AddObserver("LeftButtonReleaseEvent",self.middleButtonReleaseEvent)
        #self.AddObserver("KeyPressEvent",self.MyOnKeyPress)
        pass
    
    def OnKeyPress(self,a,b):
        print "key pressed"
        super(PdbInteractorStyle,self).OnKeyPress(self,a, b)
        
    def middleButtonPressEvent(self,obj,event):
        print "Middle Button pressed"
        self.OnMiddleButtonDown()
        return
 
    def middleButtonReleaseEvent(self,obj,event):
        print "Middle Button released"
        self.OnMiddleButtonUp()
        return
    

class PdbViewer:
    axis = False
    
    def __init__(self,ext_actors=None): #ext_actors is a list of any external vtkActors.
        #initializations:
        self.renderer  = vtk.vtkRenderer()
        self.window    = vtk.vtkRenderWindow()
        self.window.SetSize(1000,1000)
        self.mapper    = vtk.vtkPolyDataMapper()
        self.points    = vtk.vtkPoints()
        self.poly_data = vtk.vtkPolyData()
        self.glyph3d   = vtk.vtkGlyph3D()
        self.actor     = vtk.vtkActor()
        self.point_s   = vtk.vtkPointSource()
        self.sphere    = vtk.vtkSphereSource() 
        self.interactor= vtk.vtkRenderWindowInteractor()
        self.inter_sty = PdbInteractorStyle()
        self.axes_actor= vtk.vtkAxesActor()
        #configurations:
        self.point_s.SetNumberOfPoints(1)
        self.sphere.SetRadius(1.0)
        self.interactor.SetInteractorStyle(self.inter_sty)
        self.interactor.SetRenderWindow(self.window)
        self.axes_actor.SetTotalLength(100,100,100)
        if ext_actors:
            self.ex_actors = ext_actors
        else:
            self.ex_actors=[]
            
    def SetPdb(self,pdb_ext):
        self.pdb = pdb_ext
        atom_numbers = 0
        for atom in self.pdb:
            self.points.InsertNextPoint(atom.pos[0],atom.pos[1],atom.pos[2])
            atom_numbers = atom_numbers + 1
        #configurations:
        self.poly_data.SetPoints(self.points)
        resolution = int(math.sqrt(300000.0/atom_numbers))
        if resolution > 20:
            resolution = 20
        if resolution < 4:
            resolution = 4
        
        self.sphere.SetThetaResolution(resolution)
        self.sphere.SetPhiResolution(resolution)
         
    def RetrunRenderer(self,mode="dot",only_ext=False,background=[0.1,0.2,0.3]):
        if hasattr(self,'pdb'):
            if mode == "dot":
                self.glyph3d.SetSourceConnection(self.point_s.GetOutputPort())
            elif mode == "sphere":
                self.glyph3d.SetSourceConnection(self.sphere.GetOutputPort())
            self.glyph3d.SetInputData(self.poly_data)
            self.mapper.SetInputConnection(self.glyph3d.GetOutputPort())
            self.actor.SetMapper(self.mapper)
            if not only_ext:
                self.renderer.AddActor(self.actor)
            ## add external actors ##
            for act in self.ex_actors:
                self.renderer.AddActor(act)
            ##
            self.renderer.SetBackground(background)
            if self.axis:
                self.renderer.AddActor(self.axes_actor)
            self.renderer.ResetCamera()
        return self.renderer

    def Show(self,mode="dot",only_ext=False):
        if hasattr(self,'pdb'):
            if mode == "dot":
                self.glyph3d.SetSourceConnection(self.point_s.GetOutputPort())
            elif mode == "sphere":
                self.glyph3d.SetSourceConnection(self.sphere.GetOutputPort())
            self.glyph3d.SetInputData(self.poly_data)
            self.mapper.SetInputConnection(self.glyph3d.GetOutputPort())
            self.actor.SetMapper(self.mapper)
            self.window.AddRenderer(self.renderer)
            if not only_ext:
                self.renderer.AddActor(self.actor)
            ## add external actors ##
            for act in self.ex_actors:
                self.renderer.AddActor(act)
            ##
            self.renderer.SetBackground(0.1,0.2,0.3)
            if self.axis:
                self.renderer.AddActor(self.axes_actor)
            self.renderer.ResetCamera()
            self.window.Render()
            self.interactor.Start()
        else:
            print ('Set a pdb via SetPdb()')
            
    def AxesOn(self):
        self.axis = True

    def AxesOff(self):
        self.axis = False
        
    def ExportObj(self,file_name):
        exporter = vtk.vtkOBJExporter()
