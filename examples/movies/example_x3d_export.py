#example to export x3d file that contains essentially a contour file with a color mapping (here wall normal velocity)
#note that this file probably don't work as is as I kep the essential of a bigger file without testing it
try: paraview.simple
except: from paraview.simple import *
import sys,glob,os
import numpy

exporters=servermanager.createModule("exporters")

paraview.simple._DisableFirstRenderCameraReset()
dirname=sys.argv[1]

filein = XMLUnstructuredGridReader(FileName="dummy" )
cont = Contour( Input=filein,PointMergeMethod="Uniform Binning" )
cont.PointMergeMethod = "Uniform Binning"
cont.ContourBy = ['POINTS', 'Q_CRIT']
qval=float(sys.argv[2]) #e-3, + units

w1=XMLPolyDataWriter(Input=cont)
x3dExporter1=exporters.X3DExporter(FileName="dummy.x3d")

cont.Isosurfaces = [qval/pow(mu,2.0)*pow(rhoub/ub,2.0)*pow(utau,4)*1.0e-3]
  
filein.FileName=file
file1='%s/contourQ-%07d-q=%.1f.vtp'%(dirname,framenum,qval)
w1.FileName=file1
if (not os.path.isfile(file1)) or (force_recompute) :  w1.UpdatePipeline() #save the vtp file
cont.UpdatePipeline()
src=SetActiveSource(cont) 
view1 = GetRenderView()
proxy=GetActiveSource()
dr = Show(proxy=proxy,view=view1)
dr.ColorArrayName = ('POINT_DATA', 'U')
#this has to be extracted from Paraview
lut = GetLookupTableForArray( "U", 3, RGBPoints=[-0.2, 0.278431, 0.278431, 0.858824, -0.14136892128125444, 0.0, 0.0, 0.360784, -0.08315231671139994, 0.0, 1.0, 1.0, -0.02411074569148189, 0.0, 0.501961, 0.0, 0.03410911675377873, 1.0, 1.0, 0.0, 0.0927363946178838, 1.0, 0.380392, 0.0, 0.15136729234355134, 0.419608, 0.0, 0.0, 0.20999999999999996, 0.878431, 0.301961, 0.301961] )
lut.VectorMode = "Component" 
lut.VectorComponent=1
dr.LookupTable = lut
render=Render()

x3dExporter1.FileName=file4
x3dExporter1.SetView(view1) # <===== NEW LINE
x3dExporter1.Write()
