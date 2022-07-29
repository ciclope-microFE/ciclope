#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ciclope - Postprocessing Module
"""

import sys
# sys.path.append('/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib')
sys.path.append('/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib/python3.8/site-packages')

from paraview.simple import *

def paraview_plot(filein, fileout=None, slicenormal="X", RepresentationType="SurfaceWithEdges", crinkle=False, colorby="S_Mises", Roll=0, ImageResolution=[1280, 960], TransparentBackground=True, colormap='Viridis (matplotlib)'):
    """Save plot of field data using Paraview.

    Parameters
    ----------
    image
        Image data.
    threshold_value (optional)
        Threshold value. If empty an Otsu threshold is calculated.

    Returns
    -------
    """

    # read vtk file
    data = LegacyVTKReader(FileNames=filein)

    slice = Slice(registrationName='Slice1', Input=data)
    slice.SliceType = 'Plane'
    slice.HyperTreeGridSlicer = 'Plane'

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # draw the object
    # show data in view
    slice1Display = Show(slice, renderView1, 'GeometryRepresentation')

    # get the object bounds
    bounds = data.GetDataInformation().GetBounds()
    center = [(bounds[1] - bounds[0]) / 2, (bounds[3] - bounds[2]) / 2, (bounds[5] - bounds[4]) / 2]

    slice.SliceOffsetValues = [0.0]

    # slice.SliceType.Origin = [2., 2., 1.4]
    slice.SliceType.Origin = center

    colorby_string = ''.join(colorby)

    if fileout is None:
        import os
        filename_out_base, ext_out = os.path.splitext(filein)

    if 'x' in slicenormal.lower():
        slice.SliceType.Normal = [1.0, 0.0, 0.0]

        if crinkle:
            slice.Crinkleslice = 1

        # reset view to fit data
        renderView1.ResetCamera()

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [10000, center[1], center[2]]
        renderView1.CameraFocalPoint = center

        plot_slice(renderView1, slice1Display, filename_out_base+'_'+colorby_string+'_YZ.png', RepresentationType, colorby, Roll, ImageResolution, TransparentBackground, colormap)

    if 'y' in slicenormal.lower():
        slice.SliceType.Normal = [0.0, 1.0, 0.0]

        if crinkle:
            slice.Crinkleslice = 1

        # reset view to fit data
        renderView1.ResetCamera()

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [center[0], 10000, center[2]]
        renderView1.CameraFocalPoint = center

        plot_slice(renderView1, slice1Display, filename_out_base + '_' + colorby_string + '_XZ.png', RepresentationType, colorby, (-90+Roll), ImageResolution, TransparentBackground, colormap)

    if 'z' in slicenormal.lower():
        slice.SliceType.Normal = [0.0, 0.0, 1.0]

        if crinkle:
            slice.Crinkleslice = 1

        # reset view to fit data
        renderView1.ResetCamera()

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [center[0], center[1], 10000]
        renderView1.CameraFocalPoint = center

        plot_slice(renderView1, slice1Display, filename_out_base + '_' + colorby_string + '_XY.png', RepresentationType, colorby, Roll, ImageResolution, TransparentBackground, colormap)

    if not ('x' in slicenormal.lower()) | ('y' in slicenormal.lower()) | ('z' in slicenormal.lower()):
        raise IOError('Invalid Slice Normal.')

    return

def plot_slice(renderView1, slice1Display, fileout, RepresentationType, colorby, Roll, ImageResolution, TransparentBackground, colormap):
    """Save plot of field data using Paraview.

    Parameters
    ----------
    image
        Image data.
    threshold_value (optional)
        Threshold value. If empty an Otsu threshold is calculated.

    Returns
    -------
    """

    import time

    # change representation type
    # slice1Display.SetRepresentationType('Surface')
    slice1Display.SetRepresentationType(RepresentationType)

    # set scalar coloring
    if len(colorby) == 2:
        ColorBy(slice1Display, ('POINTS', colorby[0], colorby[1]))
    else:
        ColorBy(slice1Display, ('POINTS', colorby, colorby))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map
    if len(colorby) == 2:
        LUT = GetColorTransferFunction(colorby[0])
    else:
        LUT = GetColorTransferFunction(colorby)

    # get opacity transfer function/opacity map
    if len(colorby) == 2:
        PWF = GetOpacityTransferFunction(colorby[0])
    else:
        PWF = GetOpacityTransferFunction(colorby)

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    LUT.ApplyPreset(colormap, True)

    # get color legend/bar for LUT in view renderView1
    LUTColorBar = GetScalarBar(LUT, renderView1)

    # change scalar bar placement
    # LUTColorBar.Orientation = 'Horizontal'
    LUTColorBar.Orientation = 'Vertical'
    LUTColorBar.WindowLocation = 'UpperRightCorner'  # 'AnyLocation', 'LowerRightCorner', 'LowerLeftCorner', 'LowerCenter', 'UpperLeftCorner', 'UpperRightCorner', 'UpperCenter'
    LUTColorBar.TitleColor = [0, 0, 0]  # switch to black
    LUTColorBar.LabelColor = [0, 0, 0]  # switch to black
    LUTColorBar.TitleFontSize = 10
    LUTColorBar.LabelFontSize = 10
    LUTColorBar.ScalarBarThickness = 8

    # Data Axis visibility
    slice1Display.DataAxesGrid.GridAxesVisibility = 1

    # Data Axis Font
    slice1Display.DataAxesGrid.XTitleFontSize = 40
    slice1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]  # switch to black
    slice1Display.DataAxesGrid.XLabelFontSize = 40
    slice1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
    slice1Display.DataAxesGrid.YTitleFontSize = 40
    slice1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]  # switch to black
    slice1Display.DataAxesGrid.YLabelFontSize = 40
    slice1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
    slice1Display.DataAxesGrid.ZTitleFontSize = 40
    slice1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]  # switch to black
    slice1Display.DataAxesGrid.ZLabelFontSize = 40
    slice1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

    slice1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]  # grid color to black

    # update the view to ensure updated data information
    renderView1.Update()

    # dp = GetDisplayProperties()
    #
    # #set point color
    # dp.AmbientColor = [1, 0, 0] #red
    #
    # #set surface color
    # dp.DiffuseColor = [0, 1, 0] #blue
    #
    # #set point size
    # dp.PointSize = 2

    # #set representation
    # dp.Representation = "Surface"

    camera = GetActiveCamera()
    # camera.Elevation(45)
    camera.Roll(Roll)
    Render()

    # save screenshot
    if TransparentBackground:
        SaveScreenshot(fileout, ImageResolution=ImageResolution, TransparentBackground=1)
    else:
        SaveScreenshot(fileout, ImageResolution=ImageResolution)

    time.sleep(1)

    return