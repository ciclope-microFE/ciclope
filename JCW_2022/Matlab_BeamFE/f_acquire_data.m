function [dataset,FEMsolverPATH] = f_acquire_data(varargin)

% including required libraries
addpath ([pwd,'\Libraries_trussFEM_bone\skeleton3d-matlab-master']);
addpath ([pwd,'\Libraries_trussFEM_bone\skel2graph3d-matlab-master']);

% folder where the FEM solver (i.e. the file "ccx216.exe") is located:
FEMsolverPATH = [pwd,'\Libraries_trussFEM_bone\CL34-win64\bin\ccx'];

if isempty(varargin)
    [filename,foldername]=uigetfile('*.raw','File RAW to open','*.raw');
    dataset.info.fileraw = sprintf('%s%s',foldername,filename);
    
    % dataset info
    prompt = {'R:','C:','S:','Voxel Size [mm]:','Threshold [1-255]:'};
    dlgtitle = 'Dataset Details';
    dims = [1 25];
    definput = {' ',' ',' ','0.0195','102'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    dataset.info.R = str2num(answer{1});
    dataset.info.C = str2num(answer{2});
    dataset.info.S = str2num(answer{3});
    dataset.info.sp = str2num(answer{4});
    dataset.info.Threshold = str2num(answer{5});
    
    % dataset cropping
    prompt = {'Ri:','Rf:','Ci:','Cf:','Si:','Sf:'};
    dims = [1 25];
    definput = {'1',    num2str(dataset.info.R),'1',...
        num2str(dataset.info.C),'1',num2str(dataset.info.S)};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    dataset.crop.Ri = str2num(answer{1});
    dataset.crop.Rf = str2num(answer{2});
    dataset.crop.Ci = str2num(answer{3});
    dataset.crop.Cf = str2num(answer{4});
    dataset.crop.Si = str2num(answer{5});
    dataset.crop.Sf = str2num(answer{6});
    
    % modelling setting
    dataset.FEM.dia_type = 1;
    dataset.FEM.alphaDisp = 1/100;
else
    
end
