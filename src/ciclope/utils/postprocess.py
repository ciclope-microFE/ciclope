#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope postprocessing module
"""

import numpy as np
import math
import h5py
import os

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
    num_pixel_total_bone = 0
    num_pixel_total_bone-empty = 0

    for slice_mask, circular_mask in zip(L, circular_masks):
        num_pixel_bone = np.sum(np.logical_and(slice_mask, circular_mask))
        num_pixel_bone-empty = np.sum(np.logical_or(slice_mask, circular_mask))

        num_pixel_total_bone += num_pixel_bone
        num_pixel_total_bone-empty += num_pixel_bone-empty

    BVTV = num_pixel_total_bone / (num_pixel_total_bone-empty)

    return circular_masks, BVTV
    
def cyl_binary_mask2bvtv(mask: np.ndarray, voxel_size: float, radius: float, height: float) -> float:
    """
    Calculates the BV/TV (Bone Volume/Total Volume) ratio of a trabecular bone sample.
    
    Parameters:
        mask (np.ndarray): Binarized (3D) mask with 1 = bone, 0 = empty.
        voxel_size (float): Size of the voxel in mm.
        radius (float): Radius of the sample cylinder in mm.
        height (float): Height of the sample in mm.
    
    returns:
        float: BV/TV expressed as a percentage.
    """
    # Compute the number of bone voxels
    bone_pixel_number = np.sum(mask)
    
    # Compute the bone volume
    voxel_volume = voxel_size ** 3
    bone_volume = bone_pixel_number * voxel_volume
    
    # Compute of the geometric volume of the ideal cylinder
    ideal_bone_volume = np.pi * (radius ** 2) * height
    
    # Compute the BV/TV percentage
    bvtv = (bone_volume / ideal_bone_volume) * 100
    print(f'BV/TV = {bvtv:.2f}%')
    
    return bvtv

def count_fixed_displacements(file_path, slice_level):
    """
    Count the number of nodes fixed (bottom region) and nodes
    fixed with top_displacement (top region), based on voxel coordinates only.

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file.
    slice_level: int
        Number of planes locked at the top (just the nodes on the base of 
        the voxels, not all nodes of the voxel planes)

    Returns
    -------
    tuple of int
        nodes_z0_count, nodes_z1_count
    """
    nodes_z0_count, nodes_z1_count = 0, 0

    with h5py.File(file_path, 'r') as file:
        vs = file['/Image_Data/Voxelsize'][0]
        fixed_disp_coords = file['/Image_Data/Fixed_Displacement_Coordinates'][()]
        z_slices = fixed_disp_coords[:, 0]
        max_slice_id = np.max(z_slices)

        # Nodes at bottom (Z=0) --> slices from max_id - slice_level +1 to max_id
        z0_mask = (z_slices >= (max_slice_id - slice_level + 1)) & (z_slices <= max_slice_id)
        z0_raw = np.sum(z0_mask)
        nodes_z0_count = int(z0_raw / 3)  # each DOF counted

        # Nodes at top (-0.04) --> slices 0 to 10 inclusive
        z1_mask = (z_slices >= 0) & (z_slices <= 10)
        z1_raw = np.sum(z1_mask)
        nodes_z1_count = int(z1_raw / 3)

    return nodes_z0_count, nodes_z1_count

def reaction_forces(file_path, slice_level):
    """
    Compute reaction forces and mesh information from an HDF5 file,
    considering a region of nodes within a specified range of slices
    defined by `slice_level`.

    The function identifies the nodes whose Z-coordinate falls within
    the range corresponding to the last `slice_level` slices at the
    top of the voxel grid, sums their nodal forces, and returns related
    quantities.

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file containing mesh and solution data.
    slice_level : int
        Number of voxel levels (just the nodes on the base of the voxels,
        not all nodes of the voxel planes) used to define the locked region 
        for boundary conditions. The reaction forces are computed for all
        nodes within this region.

    Returns
    -------
    dict
        Dictionary containing:
        - 'Z_min' : float
            Lower Z limit of the locked region.
        - 'Z_max' : float
            Upper Z limit of the locked region.
        - 'total_force' : ndarray of shape (3,)
            Total reaction force vector [Fx, Fy, Fz].
        - 'F_tot' : float
            Norm of the total reaction force.
        - 'num_nodes' : int
            Total number of nodes in the mesh.
        - 'num_elements' : int
            Total number of elements in the mesh.
        - 'vs' : float
            Voxel size.
        - 'nodes_z0_count' : int
            Number of nodes constrained at Z = 0.
        - 'nodes_z1_count' : int
            Number of nodes constrained at Z = Z_max.
    """
    with h5py.File(file_path, 'r') as file:
        vs = file['/Image_Data/Voxelsize'][0]

        fixed_disp_coords = file['Image_Data/Fixed_Displacement_Coordinates'][()]
        max_slice_id = np.max(fixed_disp_coords[:, 0])

        Z_min = (max_slice_id - slice_level + 1) * vs
        Z_max = max_slice_id * vs

        coordinates = file['Mesh/Coordinates'][()]
        num_nodes = coordinates.shape[0]

        z_indices = np.where(
            (coordinates[:, 2] >= Z_min) & (coordinates[:, 2] <= Z_max)
        )[0]

        nodal_forces = file['Solution/Nodal forces'][()]
        forces_filtered = nodal_forces[z_indices, :]
        total_force = np.sum(forces_filtered, axis=0)
        F_tot = np.linalg.norm(total_force)

        elements = file['Mesh/Elements'][()]
        num_elements = elements.shape[0]

        nodes_z0_count, nodes_z1_count = count_fixed_displacements(file_path, slice_level)

        return {
            'Z_min': Z_min,
            'Z_max': Z_max,
            'total_force': total_force,
            'F_tot': F_tot,
            'num_nodes': num_nodes,
            'num_elements': num_elements,
            'vs': vs,
            'nodes_z0_count': nodes_z0_count,
            'nodes_z1_count': nodes_z1_count
        }
    
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