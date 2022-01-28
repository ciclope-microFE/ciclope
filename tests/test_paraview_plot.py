#!/usr/bin/pvpython
from paraview.simple import *
import time

#read a vtp
reader = LegacyVTKReader(FileNames="/home/gianthk/PycharmProjects/CT2FE/test_data/steel_foam/B_matrix_tetraFE_Nlgeom.10.vtk")

# #position camera
view = GetActiveView()
if not view:
    # When using the ParaView UI, the View will be present, not otherwise.
    view = CreateRenderView()
view.CameraViewUp = [0, 0, 1]
view.CameraFocalPoint = [0, 0, 0]
view.CameraViewAngle = 45
view.CameraPosition = [5,0,0]

# slicer = Slice(Input=reader, SliceType="Plane")
# slicer.SliceType.Origin = [0, 0, 0]
# slicer.SliceType.Normal = [0, 0, 1]
#
# # To render the result, do this:
# Show(slicer)
# Render()

#draw the object
Show()

#set the background color
view.Background = [1,1,1]  #white

#set image size
view.ViewSize = [800, 800] #[width, height]

dp = GetDisplayProperties()

#set point color
dp.AmbientColor = [1, 0, 0] #red

#set surface color
dp.DiffuseColor = [0, 1, 0] #blue

#set point size
dp.PointSize = 2

#set representation
dp.Representation = "Surface"

Render()

#save screenshot
WriteImage("pippo.png")