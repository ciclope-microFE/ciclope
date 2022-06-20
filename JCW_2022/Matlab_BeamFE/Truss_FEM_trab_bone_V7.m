% Script: "Truss_FEM_trabecular_bone_V6.m" - Matlab script
% Created by Martino pani - martino.pani@port.ac.uk
% Last change: 10th Februrary 2020
% -------------------------------------------------------------------------
%
% This cript builds a Truss FEM model from a stack of microCT images.
%
% INPUT:
%
% - microCT dataset (8bit raw file)
%
%
% OUTPUT:
%
% - CalculiX ".inp" file for the static analyisis (uniaxial compression, 
%   cinematic BCs) and the related linear Buckling analysis (loadr from the
%   nodal reation retrieved from the static analysis)
%
% A set if graphical outputs are available. By default the graphics is
% disabled and it is recommended to enable it only for SMALL DATASETS.
%
%
