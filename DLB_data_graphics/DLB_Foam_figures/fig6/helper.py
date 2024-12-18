from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk

def generate_target_grid(x,z,xspacing,yspacing):
    xi = np.arange(np.min(x),np.max(x),xspacing)
    zi = np.arange(np.min(z),np.max(z),yspacing)
    xi,zi = np.meshgrid(xi,zi)
    return xi,zi
def interpolate_grid(x,z,field,xi,zi):
    data = griddata((x,z),field,(xi,zi),method='linear')
    return data

def readVTKFile(path):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(path)
        reader.Update()
        output = reader.GetOutput()
        return output
def readScalarField(path,fieldName):
        vtk_f = readVTKFile(path)
        s1 = vtk_f.GetPointData().GetArray(fieldName)
        return vtk_to_numpy(s1)

