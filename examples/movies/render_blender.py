import bpy,math,sys
import os

def DeleteObject(name_obj):
    bpy.ops.object.select_all(action='DESELECT')
    if name_obj in bpy.data.objects :
        bpy.data.objects[name_obj].select=True
        bpy.ops.object.delete()

#essentially parameters to find the right file
index=int(sys.argv[-1])
framenum=int(sys.argv[-2])
qval=float(sys.argv[-3]) #e-3, + units
dirname=sys.argv[-4]

#scenes to render
scenes=['scene1-nodata.blend','scene2-nodata.blend','scene3-nodata.blend','scene4-nodata.blend']
scenes=['scene1-nodata.blend']

#x3d file to import
f='%s/contourQ-%07d-q=%.1f.x3d'%(dirname,framenum,qval)

print('processing file '+ f)
print('frame %03d'%framenum)

#default name of the imported object in blender

name_obj="ShapeIndexedFaceSet"

for s in scenes:
    print(s)
    out_dir='%s/%s-%.1f'%(dirname,s.split('.')[0],qval)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    #open the scene
    bpy.ops.wm.open_mainfile(filepath=s)
    #import the x3d
    bpy.ops.import_scene.x3d(filepath=f)
    #remove lights coming from paraview
    for obj in bpy.data.objects:
      if obj.name[:4]=="TODO":
        print("deactivating light:, ",obj.name)
        bpy.data.objects[obj.name].hide_render=True
    
    contour = bpy.data.objects[name_obj]

    #put the object at the right place
    contour.location=(0.,-3.5,1.03)

    scn = bpy.context.scene

    for face in contour.data.polygons:
        face.use_smooth = True

    scn.frame_current=index
    scn.render.resolution_percentage=50

    #here 'MatQ' has been saved in the scene file: trick -> assign it to whatever objects that won't be rendered
    # in order to make sure it will stay in the file .blend
    contour.data.materials.append(bpy.data.materials['MatQ'])
    contour.data.materials['MatQ'].use_vertex_color_paint=True


    #exemple of possible change in the scripts
    #contour1.data.materials['MatQ'].diffuse_intensity=1.0 #0.7
    #contour1.data.materials['MatQ'].alpha=0.8
    #contour1.data.materials['MatQ'].diffuse_color=(0,0,0)
    #contour1.data.materials['MatQ'].specular_color=(1,1,1)
    #contour1.data.materials['MatQ'].use_transparency=False
    #contour1.data.materials['MatQ'].raytrace_mirror.use=True
    #contour1.data.materials['MatQ'].raytrace_mirror.fresnel=0.5


    #lamp=bpy.data.objects['Lamp']
    #lamp.data.shadow_color=(0.7,0.7,0.7) 
    #lamp.data.shadow_method='NOSHADOW'

    #set path and image name
    scn.render.filepath=out_dir+'/%07d.png'%framenum
    
    #render
    bpy.ops.render.render(write_still=True)
