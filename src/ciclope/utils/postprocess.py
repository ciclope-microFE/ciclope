#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope postprocessing module
"""

import numpy as np
import math
import h5py

try:
    from paraview.simple import *
except ImportError:
    import warnings
    warnings.warn("ParaView.simple is required!", RuntimeWarning)

def paraview_plot(filein, fileout=None, slicenormal="XYZ", RepresentationType="Surface", Crinkle=False, ColorBy="S_Mises", Roll=0, ImageResolution=[1280, 960], TransparentBackground=True, ColorMap='Viridis (matplotlib)'):
    """Plot field data using ParaView.
    ParaView must be installed and a link to its python library must be added to your system path.

    Parameters
    ----------
    filein : str
        Input data (VTK or other 3D ParaView file).
    fileout : str
        Output image file name.
    slicenormal : str
        Any combination of 'X', 'Y', and 'Z'. Default='XYZ'.
    RepresentationType : str
        'Surface', 'SurfaceWithEdges', 'Volume', 'Points', 'Feature Edges', or '3D Glyphs'. Default='Surface'.
    Crinkle : bool
        Crinkle the slice. Default=False.
    ColorBy : str
        Field name for coloring. Default='S_Mises'
    Roll : int
        View roll angle. Default=0.
    ImageResolution : int
        Output image resolution [X, Y] in pixels. Default=[1280, 960].
    TransparentBackground : bool
        Transparent background. Default=True.
    Colormap : str
        Default = Viridis.
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

    colorby_string = ''.join(ColorBy)

    if fileout is None:
        import os
        filename_out_base, ext_out = os.path.splitext(filein)

    if 'x' in slicenormal.lower():
        slice.SliceType.Normal = [1.0, 0.0, 0.0]

        if Crinkle:
            slice.Crinkleslice = 1

        # reset view to fit data
        renderView1.ResetCamera()

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [10000, center[1], center[2]]
        renderView1.CameraFocalPoint = center

        plot_slice(renderView1, slice1Display, filename_out_base +'_' + colorby_string +'_YZ.png', RepresentationType, ColorBy, Roll, ImageResolution, TransparentBackground, ColorMap)

    if 'y' in slicenormal.lower():
        slice.SliceType.Normal = [0.0, 1.0, 0.0]

        if Crinkle:
            slice.Crinkleslice = 1

        # reset view to fit data
        renderView1.ResetCamera()

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [center[0], 10000, center[2]]
        renderView1.CameraFocalPoint = center

        plot_slice(renderView1, slice1Display, filename_out_base + '_' + colorby_string + '_XZ.png', RepresentationType, ColorBy, (-90 + Roll), ImageResolution, TransparentBackground, ColorMap)

    if 'z' in slicenormal.lower():
        slice.SliceType.Normal = [0.0, 0.0, 1.0]

        if Crinkle:
            slice.Crinkleslice = 1

        # reset view to fit data
        renderView1.ResetCamera()

        # changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [center[0], center[1], 10000]
        renderView1.CameraFocalPoint = center

        plot_slice(renderView1, slice1Display, filename_out_base + '_' + colorby_string + '_XY.png', RepresentationType, ColorBy, Roll, ImageResolution, TransparentBackground, ColorMap)

    if not ('x' in slicenormal.lower()) | ('y' in slicenormal.lower()) | ('z' in slicenormal.lower()):
        raise IOError('Invalid Slice Normal.')

    return

def plot_slice(renderView1, slice1Display, fileout, RepresentationType, colorby, Roll, ImageResolution, TransparentBackground, ColorMap):
    """Save plots using Paraview.

    Parameters
    ----------
    renderView1
    slice1Display
    fileout
    RepresentationType
    colorby
    Roll
    ImageResolution
    TransparentBackground
    ColorMap
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
    LUT.ApplyPreset(ColorMap, True)

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

def calculate_total_force(filename_dat):
    """
    Calculate the total force from a data file.

    Parameters:
    - filename_dat (str): The path to the data file containing force components.

    Returns:
    - float or None: The total force calculated from the data file, or None if an error occurs.
    """
    try:
        with open(filename_dat, 'r') as file:
            lines = file.readlines()

        # Check if there are enough lines
        if len(lines) >= 4:
            # Extract the data line (fourth line) and split it
            data_line = lines[3].split()

            # Check if data_line contains enough elements
            if len(data_line) >= 3:
                # Extract the values of the three components
                fx = float(data_line[0])
                fy = float(data_line[1])
                fz = float(data_line[2])

                # Calculate the total force
                total_force = math.sqrt(fx**2 + fy**2 + fz**2)

                return total_force
            else:
                raise ValueError("Not enough data in the lines of force component values.")
        else:
            raise ValueError("The file does not contain enough lines.")

    except FileNotFoundError:
        print("The specified file does not exist.")
        return None
    
def circular_masks_BVTV(L, diameter, pixel_spacing_mm):
    """
    Create circular masks and calculate BVTV.

    Parameters
    ----------
    L : list of numpy.ndarray
        List of input slices.

    diameter : float
        Diameter of the circle.

    pixel_spacing_mm : float
        Pixel spacing in millimeters.

    Returns
    -------
    circular_masks : list of numpy.ndarray
        List of circular masks.

    BVTV : float
        Bone volume-to-total volume ratio.
    """
    # Calculate the circle radius in pixels
    circle_radius_pixel = diameter / (2 * pixel_spacing_mm)

    # Initialize a list for circular masks
    circular_masks = []

    # Calculate a circular mask for each slice in L
    for slice_mask in L:
        # Find the center of the slice (assuming it is at the center)
        center_x, center_y = slice_mask.shape[1] // 2, slice_mask.shape[0] // 2

        # Create a grid of coordinates
        y, x = np.ogrid[:slice_mask.shape[0], :slice_mask.shape[1]]

        # Calculate the circular mask for this slice
        circular_mask = ((x - center_x) ** 2 + (y - center_y) ** 2 <= circle_radius_pixel ** 2)

        # Add the circular mask to the list
        circular_masks.append(circular_mask)

    # Calculate the BVTV of the entire model
    num_pixel_osso_totale = 0
    num_pixel_osseo_vuoto_totale = 0

    for slice_mask, circular_mask in zip(L, circular_masks):
        num_pixel_osso = np.sum(np.logical_and(slice_mask, circular_mask))
        num_pixel_osseo_vuoto = np.sum(np.logical_or(slice_mask, circular_mask))

        num_pixel_osso_totale += num_pixel_osso
        num_pixel_osseo_vuoto_totale += num_pixel_osseo_vuoto

    BVTV = num_pixel_osso_totale / (num_pixel_osseo_vuoto_totale)

    return circular_masks, BVTV

 def reaction_forces(file_path, vs):
    """
    Calculate total reaction force and Z value from an HDF5 file.

    Parameters:
    file_path (str): Path to the HDF5 file.
    vs (float): Voxel size.

    Returns:
    Z_value (float): The calculated Z value.
    total_force (numpy.ndarray): The total force (fx, fy, fz).
    F_tot (float): Magnitude of the total force.
    """
    with h5py.File(file_path, 'r') as file:
        # Access the Fixed_Displacement_Coordinates dataset within the Image_Data group
        fixed_disp_coords = file['Image_Data/Fixed_Displacement_Coordinates'][()]
        
        # Find the slice ID with the highest numerical value
        max_slice_id = np.max(fixed_disp_coords[:, 0])
        
        # Calculate Z_value
        Z_value = max_slice_id * vs
        
        print(f"Z_value: {Z_value}")
        print()
        
        # Read the node coordinates from the Mesh group
        coordinates = file['Mesh/Coordinates'][()]
        
        # Identify nodes with the specified Z value
        z_indices = np.where(coordinates[:, 2] == Z_value)[0]
        
        # Assuming we can access the nodal force data in some manner
        nodal_forces = file['Solution/Nodal forces'][()]
        
        # Assuming nodal_forces has a one-to-one mapping with coordinates
        # Extract forces corresponding to Z indices
        forces_filtered = nodal_forces[z_indices, :]
        
        # Calculate the sum of forces for fx, fy, and fz
        total_force = np.sum(forces_filtered, axis=0)
        
        print(f"Total force (fx, fy, fz) for set NODES_Z0 and time  0.1000000E+01: {total_force}")
        print()
        
        # Apply Pythagoras' theorem to calculate the total force
        F_tot = np.sqrt(np.sum(total_force**2))
        
        print(f"F_tot: {F_tot:.2f} N")
        
        return Z_value, total_force, F_tot
    
def sample_height(input_folder, vs):
    """
    Calculate the height of the sample in millimeters based on the number of .tif slices.

    Parameters
    ----------
    input_folder : str
        The path to the folder containing the .tif slices.
    voxel_size : float
        The size of a voxel in millimeters.

    Returns
    -------
    float
        The height of the sample in millimeters.
    """
    # Count the number of .tif files in the specified folder
    num_slices = len([f for f in os.listdir(input_folder) if f.endswith('.tif')])

    # Calculate the height of the sample
    height_mm = num_slices * vs

    return height_mm