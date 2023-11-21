import vtk as vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt
import sys
from vtk.numpy_interface import dataset_adapter as dsa
from collections import OrderedDict
import scipy.interpolate as si
import files_pvtu as fvtu
import os
import scipy.signal as sp
import scipy as sc


if not os.path.exists("BL"):
    os.makedirs("BL") 
outfolder = "BL/"


############################################################
# FieldF - plane slice obtained from IC3  
Wall = "Folder_containing_SurfPVTU/"
FieldF = "Folder_containing_slicePVTU/"

Uinf = None                                         #Free stream velocity
Rhoinf = None                                       #Free stream density
intrp = None                                        #High res line for BL; specify int value
height = 0.02                                       #Height for the profiles 
norm_surf = True                                    #Surface normal profiles
pnts = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]        #Points for the profiles
zero_correction = True                              #Correction for the surface height
# Figure parameter
Figure = True                                       #For a default latex style figure
fig_delta = True                                    #To include d99 profile in the default figure
###########################################################




wall_file = fvtu.pvtu_files(Wall)[-1]
wall_reader = vtk.vtkXMLPUnstructuredGridReader()
wall_reader.SetFileName(Wall+wall_file)
wall_reader.Update()
wall = wall_reader.GetOutput()

field_file = fvtu.pvtu_files(FieldF)[-1]
field_reader = vtk.vtkXMLPUnstructuredGridReader()
field_reader.SetFileName(FieldF+field_file)
field_reader.Update()
field = field_reader.GetOutput()

# Set the origin for cut plane
center = field_reader.GetOutput().GetCenter()
plane = vtk.vtkPlane()
plane.SetOrigin(center)
plane.SetNormal(0,0,1) # set norm here

#set wall mapper
wall_map = vtk.vtkCellDataToPointData()
wall_map.AddInputData(wall)
wall_map.Update()
field_map = vtk.vtkCellDataToPointData()
field_map.AddInputData(field)
field_map.Update()


wall_slice = vtk.vtkCutter()
wall_slice.SetCutFunction(plane)
wall_slice.SetInputConnection(wall_map.GetOutputPort())
wall_slice.Update()
field_slice = vtk.vtkCutter()
field_slice.SetCutFunction(plane)
field_slice.SetInputConnection(field_map.GetOutputPort())
field_slice.Update()

#data aquisition
wall_data = wall_slice.GetOutput().GetPoints().GetData()
wall_coords = np.array([wall_data.GetTuple3(x)
                    for x in range(wall_data.GetNumberOfTuples())])                   
# sort to x
xsort = np.argsort(wall_coords[:,0])

#wall coordinates
xwall = wall_coords[xsort,0]
ywall = wall_coords[xsort,1]
zwall = wall_coords[xsort,2]


###################
xbl,ybl,uxbl=[],[],[]
step = 0

for i in pnts:

    xpt = i*np.max(xwall)


    def find_nearest(array, value):
        n = [abs(i-value) for i in array]
        idx = n.index(min(n))
        return idx
    xidx = find_nearest(xwall,xpt)

    if norm_surf == True:
        norms = [(xwall[xidx-1]-xwall[xidx+1]),(ywall[xidx-1]-ywall[xidx+1]),0]
    else:
        norms = [1,0,0]

    BLplane = vtk.vtkPlane()
    BLplane.SetOrigin(xpt,0,0)
    BLplane.SetNormal(*norms)

    BLslice = vtk.vtkCutter()
    BLslice.SetCutFunction(BLplane)
    BLslice.SetInputConnection(field_slice.GetOutputPort())
    BLslice.Update()
    slice_data = BLslice.GetOutput().GetPoints().GetData()
    BLcoord = np.array([slice_data.GetTuple3(x)
                    for x in range(slice_data.GetNumberOfTuples())])
    u = vtk_to_numpy(BLslice.GetOutput().GetPointData().GetArray("U_AVG"))
    rho = vtk_to_numpy(BLslice.GetOutput().GetPointData().GetArray("RHO_AVG"))   
    # sort to x
    ysort = np.argsort(BLcoord[:,1])            
    xb =  BLcoord[ysort,0]
    yb =  BLcoord[ysort,1]
    ux = u[ysort,0]
    rho = rho[ysort]
    
    
    # Referece size of the array to have same shape throughout
    if (step==0):
        idx = np.where(np.logical_and(yb>=ywall[xidx], yb<=height+ywall[xidx]))
        ref_size = np.size(idx)
    else:
        ldx = np.array(np.where(yb>=ywall[xidx]))
        idx = ldx[0][0:ref_size]

    x = xb[idx]
    y = yb[idx]

    if Uinf == None:
        Uinf = np.mean(ux[-20:-1])
    if Rhoinf == None:
        Rhoinf = np.mean(rho[-20:-1])

    ux = ux[idx]/Uinf 
    
    # Interpolate the data using a cubic spline to len intrp 
    if isinstance(intrp, int)==True:
        new_x = np.linspace(x.min(), x.max(), intrp)
        new_y = np.linspace(y.min(), y.max(), intrp)
        new_ux = sc.interpolate.interp1d(y, ux, kind='cubic')(new_y)
    else:
        new_x = x
        new_y = y
        new_ux = ux
        
    if zero_correction == True:
        new_y = new_y-ywall[xidx]

    if(step == 0):
        xbl  = np.append(xbl,  new_x, axis=0)
        ybl  = np.append(ybl,  new_y, axis=0)
        uxbl = np.append(uxbl, new_ux, axis=0)
    else:
        xbl = np.vstack([xbl, new_x])
        ybl = np.vstack([ybl, new_y])
        uxbl = np.vstack([uxbl, new_ux])
    step = step+1



import matplotlib as mpl
font = {'family':'serif','weight':'medium','size':18}
mpl.rc('font',**font)
mpl.rc('text',usetex=True)
mpl.rcParams['axes.linewidth'] = 1.0
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpcol
import matplotlib.font_manager as fm





    

# High res BL stats
resolution = 10 # for high res stats
delta, deltas, theta = [],[],[]
step = 0
res = np.linspace(xwall.min()+0.05,0.99,resolution)

for i in res:
    xpt = i*np.max(xwall)

    def find_nearest(array, value):
        n = [abs(i-value) for i in array]
        idx = n.index(min(n))
        return idx
    
    xidx = find_nearest(xwall,xpt)

    if norm_surf == True:
        norms = [(xwall[xidx-1]-xwall[xidx+1]),(ywall[xidx-1]-ywall[xidx+1]),0]
    else:
        norms = [1,0,0]

    BLplane = vtk.vtkPlane()
    BLplane.SetOrigin(xpt,0,0)
    BLplane.SetNormal(*norms)

    BLslice = vtk.vtkCutter()
    BLslice.SetCutFunction(BLplane)
    BLslice.SetInputConnection(field_slice.GetOutputPort())
    BLslice.Update()
    slice_data = BLslice.GetOutput().GetPoints().GetData()
    BLcoord = np.array([slice_data.GetTuple3(x)
                    for x in range(slice_data.GetNumberOfTuples())])
    u = vtk_to_numpy(BLslice.GetOutput().GetPointData().GetArray("U_AVG"))
    rho = vtk_to_numpy(BLslice.GetOutput().GetPointData().GetArray("RHO_AVG"))   
    # sort to x
    ysort = np.argsort(BLcoord[:,1])            
    xb =  BLcoord[ysort,0]
    yb =  BLcoord[ysort,1]
    ux = u[ysort,0]
    rho = rho[ysort]

        # Referece size of the array to have same shape throughout
    if (step==0):
        idx = np.where(np.logical_and(yb>=ywall[xidx], yb<=height+ywall[xidx]))
        ref_size = np.size(idx)
    else:
        ldx = np.array(np.where(yb>=ywall[xidx]))
        idx = ldx[0][0:ref_size]
        
    x = xb[idx]
    y = yb[idx]
    ux = ux[idx]
    rho = rho[idx]



    if Uinf == None:
        Uinf = np.mean(ux[-20:-1])
    if Rhoinf == None:
        Rhoinf = np.mean(rho[-20:-1])
    
    

    
    # point for integration
    for j in range(1,np.shape(ux)[0]):
        if (ux[j]>=0.99*Uinf):
            blindx = j
            break
    
    # Delta    
    dell = (y[blindx-1]+y[blindx]+y[blindx+1])/3.
    
    #Delta_S
    ds1 = np.trapz((1-(rho[:blindx-1]*ux[:blindx-1])/(Uinf*Rhoinf)),y[:blindx-1])
    ds2 = np.trapz((1-(rho[:blindx]*ux[:blindx])/(Uinf*Rhoinf)),y[:blindx])
    ds3 = np.trapz((1-(rho[:blindx+1]*ux[:blindx+1])/(Uinf*Rhoinf)),y[:blindx+1])
    del_s = (ds1+ds2+ds3)/3
    
    
    # theta
    th1 = np.trapz((rho[:blindx-1]*ux[:blindx-1])/(Rhoinf*Uinf)*(1-ux[:blindx-1]/Uinf),y[:blindx-1])
    th2 = np.trapz((rho[:blindx]*ux[:blindx])/(Rhoinf*Uinf)*(1-ux[:blindx]/Uinf),y[:blindx])
    th3 = np.trapz((rho[:blindx+1]*ux[:blindx+1])/(Rhoinf*Uinf)*(1-ux[:blindx+1]/Uinf),y[:blindx+1])
    tht = (th1+th2+th3)/3
        
    if zero_correction==True:
        delta.append(dell-ywall[xidx])
        deltas.append(del_s)
        theta.append(tht)
    else:
        delta.append(dell)
    step = step+1


#Saves
np.save(outfolder+"y",ybl)
np.save(outfolder+"x",xbl)
np.save(outfolder+"U_profile",uxbl)

params = np.vstack((delta,deltas,theta))
np.save(outfolder+"parameters",params)

if Figure==True:
    fig, ax = plt.subplots(figsize=(8,3))
    offset=0
    for n in range(np.shape(pnts)[0]):
        plt.plot(uxbl[n,:]*0.1+offset,ybl[n,:],'k--')
        offset = offset+0.1
    plt.ylim(0,height)
    plt.xlabel('$x/C$')
    plt.ylabel('$y/C$')
    plt.ylim(0,height)
    plt.xlim(0,1)
    
    plt.plot(res,deltas,'r-',label='$\delta^*$')
    plt.plot(res,theta,'g:',label='$\theta$')
    if fig_delta ==True:
        plt.plot(res,delta,'k-',label="$\delta$")
    if zero_correction==False:
        plt.plot(xwall,ywall)
    
    plt.tight_layout()
    plt.savefig("demo.png")
   
    plt.show()    
    






















