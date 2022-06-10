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
% REQUIRED LIBRARIES:
%
% - Skeleton 3D from https://github.com/phi-max/skeleton3d-matlab based on [1]
% - Skel2Graph 3D from https://github.com/phi-max/skel2graph3d-matlab based
%   on [1]. Line 188 in the function “Skel2Graph3D.m” has been changed from
%   the original version (i.e. the matrix is then initialised as sparse
%   instead of full), to avoid error in memory allocation.

% - FEM solver Calculix / Abaqus
% 
% References: 
%
% [1] Kerschnitzki, Kollmannsberger et al.,
% "Architecture of the osteocyte network correlates with bone material quality."
% Journal of Bone and Mineral Research, 28(8):1837-1845, 2013.
%
% about beam element in calculix
% http://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node53.html
% -------------------------------------------------------------------------
% 
% Procedure:
%
% 1 -   Read the RAW file (8bit) of the microCT dataset
%
% 2 -   Binarise the dataset. A fixed global threshold is used by default,
%       but different approaches are available.
%
% 3 -   Identify the skeleton associated to the binarised dataset (by using
%       the functio "Skeleton3D" of the "Skeleton 3D library" and calculating
%       its distance from the bony surface.
%
% 4 -   retrieve the graph associated to the skeleton (the function
%       "Skel2Graph3D" is use) and select the biggest graph (i.e. the 
%       brances associated to the isolated cluster of bone voxels are removed).
%
% 5 -   crate a truss element of each branch of the graph. Both matrixes for
%       nodal coordinates and nodal connectivity are created 
%
% 6 -   set the .inp file for the truss model to be solved in
%       Abaqus/calculix; firt the static analysis associated to kinematic
%       BCs on the nodes of the top surfaces is solved; nodal reations are
%       extracted, then applied as loads for the associated buckling analysis.
%       Each model is solved by calling the executable script and the buckling
%       proportionality factor is extracted. (Calculix is assumed)
%
% -------------------------------------------------------------------------
clc; clear all; close all; more off;

% % including required libraries
% addpath ([pwd,'\Skeleton_buckling_bone\skeleton3d-matlab-master']);
% addpath ([pwd,'\Skeleton_buckling_bone\skel2graph3d-matlab-master']);
% 
% % folder where the FEM solver is located
% %FEMsolverPATH = [pwd,'\CL34-win64\bin\ccx'];

% -------------------------------------------------------------------------
% 1 -   Read the RAW file (8bit) of the microCT dataset
% -------------------------------------------------------------------------
[dataset,FEMsolverPATH] = f_acquire_data;
%  file ".raw" to open
fileraw = dataset.info.fileraw;

% Details about the dataset (spacing and size)
sp = dataset.info.sp;    % isotropic voxel size (spacing)
C = dataset.info.C;
R = dataset.info.R;
S = dataset.info.S;

% parameters for cropping
Ri = dataset.crop.Ri;
Rf = dataset.crop.Rf;
Ci = dataset.crop.Ci;
Cf = dataset.crop.Cf;
Si = dataset.crop.Si;
Sf = dataset.crop.Sf;

% fixed global threshold for binarisation
T = dataset.info.Threshold;

% shrinkage percentage (of the model length)
alphaDisp = dataset.FEM.alphaDisp;
% Truss diameter type: 1(min), 2(max), 3(mean) value of trabecular
% thickness
dia_type = dataset.FEM.dia_type;

% reading data
fid = fopen(fileraw,'rb');
data = fread(fid,R*C*S,'uint8');
fclose(fid);

% PERCENTAGE OF CONSTRAIEND LENGRTH
clp = 5/100;

% -------------------------------------------------------------------------
% 2 -   Binarise the dataset.
% -------------------------------------------------------------------------
% we do want bone = 1 and background = 0;
data = -1.*(data-255);
I = reshape(data,R,C,S);
% cropping - if appropriate
I = I(Ri:Rf,Ci:Cf,Si:Sf);

% T = fixed global threshold

% Compute a global threshold using the histogram counts.
% counts = imhist(I, 256);
% T = otsuthresh(counts)*256;
% another method
%T = graythresh(I./255)*255;

% Binarize image using computed threshold.
BW = imbinarize(I,T);
%BW = imbinarize(I,'global');%,'ForegroundPolarity','dark','Sensitivity',0.4);

% -------------------------------------------------------------------------
% 3 - SKELETONISATION AND TRABECULAR THICKNESS ASSESSMENT
% -------------------------------------------------------------------------
skel = Skeleton3D(BW);

% DISTANCE MAPPING (min dist between skeleton ans bone surface) -----------
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% coordinated of the voxels forming the skeleton
IDskel = find(skel); % position of the skeleton on the matrix
[x,y,z] = ind2sub([w,l,h],IDskel);

% nearest distance beteen each point of the skeleton and the surface
C = -1.*(BW-1); % inverting the bins-> bone = 0 and background = 1
D = bwdist(C);  % default = Eucledean distance

skelD = zeros(size(skel));
% associating each elements of the skeleton its distance from the bone surface
skelD(IDskel) = D(IDskel); 

% -------------------------------------------------------------------------
% 4 - MAKING THE GRAPH
% -------------------------------------------------------------------------
% calculate the adjacent matrix for the skeleton
[A,node,link] = Skel2Graph3D(skel,0);
% making the Graph form the adjacent matrix
G = graph(A);
nlink = size(link,2);
Truss_diameters = zeros(nlink,3);
for k = 1:nlink
    diameters = 2*double(D(link(k).point))*sp;
    Truss_diameters(k,:) = [min(diameters), max(diameters), mean(diameters)];
end

% identify the connected components
bins = conncomp(G);
group = mode(bins);
idx = find(bins==group);

% -------------------------------------------------------------------------
% 5 - crate a truss element of each brunch of the graph.
% -------------------------------------------------------------------------
% Initialisation of the matrice for the nodal coordinates [IDn, x,y,z]
Nodes = zeros(length(idx),4);
list_links = zeros(length(link),1);
for n = 1:length(idx)
    Nodes(n,:) = [idx(n), node(idx(n)).comy, node(idx(n)).comx, node(idx(n)).comz];
    list_links(node(idx(n)).links) = 1;
end

% matrix for the nodal occectivity of each element [n1, n2]
idx = find(list_links);
Elements = [link(idx).n1;link(idx).n2]';

Truss_diameters = Truss_diameters(idx,:);
% renumbering "Elements"
for n = 1:size(Nodes,1);
    idx = find(Elements == Nodes(n,1));
    Elements(idx) = n;
    Nodes(n,1)=n;
end
% scaling nodal coordinates accoring to the voxel size
Nodes(:,2:4)=Nodes(:,2:4).*sp;

% removing duplicated truss elements
Elements = unique(Elements,'row');
% for e = 1:size(Elements,1)
%     lug(e) = norm(Nodes(Elements(e,2),2:4)-Nodes(Elements(e,1),2:4));
% end

% ====================== CALCULIX ANALYSIS ================================
% 6.1 - Writing the .inp file for Calculix - STATIC ANALYSIS
% =========================================================================
% setting the file to print
ps = findstr('\',fileraw(1:end-4));
filenameinp = fileraw(ps(end)+1:end-4);
fileinp = sprintf('.\\FEM_models\\%s_Static_model.inp',filenameinp);
fid = fopen(fileinp,'wt');

% Adding the midside point to each beam
Nmax = size(Nodes,1);
Nodes = [Nodes;zeros(size(Elements,1),4)];
Elements = [Elements,zeros(size(Elements,1),1)];
for e = 1:size(Elements,1)
    Nodes(Nmax+e,2:4) = 0.5.*(Nodes(Elements(e,1),2:4)+ Nodes(Elements(e,2),2:4));
    Elements(e,3) = Nmax+e;
end

%%[Diameters,b,idxTruss]=unique(Truss_diameters(:,2));
% - Diameters = list of diameters in ascending order
% - idxTruss = position of the truss diameter in the Diameter array
% -------------------------------------------------------------------------

% List of Nodes
fprintf(fid,'*NODE,NSET=Nall\n');
for n = 1:size(Nodes,1)
    fprintf(fid,'%i,%f,%f,%f\n',n,Nodes(n,2:4));
end
% List of Elements
fprintf(fid,'*ELEMENT,TYPE=B32,ELSET=Eall\n');
for e = 1:size(Elements,1)
    fprintf(fid,'%i,%i,%i,%i\n',e,Elements(e,[1,3,2]));
end
% Element sets (used afterward to specify the truss corss sections)
for e = 1:size(Elements,1)
    fprintf(fid, '*ELSET,ELSET=Elset%i\n',e);
    %%idx = find(idxTruss==e);
    %%fprintf(fid, '%i\n',idx);
    fprintf(fid, '%i\n',e);
end
% Node sets (nodes at Top and Bottom region for the applying BCs)
% nodes at the bottom region
Zmin = min(Nodes(:,4));
Zmax = max(Nodes(:,4));
DeltaZ = Zmax-Zmin;
idx_min = find (Nodes(:,4)<=(Zmin+clp*DeltaZ));
fprintf(fid, '*NSET,NSET=BOTTOM\n');
fprintf(fid, '%i\n',idx_min);
% nodes at the top region
idx_max = find (Nodes(:,4)>=(Zmin+(1-clp)*DeltaZ));
fprintf(fid, '*NSET,NSET=TOP\n');
fprintf(fid, '%i\n',idx_max);

% Boundary Conditions
fprintf(fid,'*BOUNDARY\n');
fprintf(fid, 'BOTTOM,1,3\n');
fprintf(fid, 'TOP,1,2\n');
fprintf(fid, 'TOP,3,3,%f\n',-DeltaZ*alphaDisp);
%fprintf(fid, 'TOP,4,6\n');
% Setting the material
fprintf(fid, '*MATERIAL,NAME=BONE\n');
fprintf(fid, '*ELASTIC\n');
fprintf(fid, '19E3,.3\n');

% Beam sections
%%for e = 1:length(Diameters)
for e = 1:size(Elements,1)   
    fprintf(fid, '*BEAM SECTION, ELSET=Elset%i,MATERIAL=BONE,SECTION=CIRC\n',e);
    %%fprintf(fid, '%5.3f,%5.3f\n',Diameters(e),Diameters(e));
    % min = 1, max = 2, mean = 3
    fprintf(fid, '%5.3f,%5.3f\n',Truss_diameters(e,dia_type),Truss_diameters(e,dia_type));
    % calculating the orientation
    %%truss_dir = [1,0,10];
    V1 = Nodes(Elements(e,2),2:4)-Nodes(Elements(e,1),2:4);
    V2 = V1+(rand(1,3).*norm(V1));
    truss_dir = cross(V1,V2);
    truss_dir = truss_dir./norm(truss_dir);
    fprintf(fid, '%4.3f,%4.3f,%4.3f\n',truss_dir);
end
% setting the analysis
fprintf(fid, '*STEP\n');
fprintf(fid, '*STATIC\n');
fprintf(fid, '1,1\n');
% Boundary Conditions
fprintf(fid,'*BOUNDARY\n');
fprintf(fid, 'BOTTOM,1,3\n');
fprintf(fid, 'TOP,1,2\n');
fprintf(fid, 'TOP,3,3,%f\n',-DeltaZ*alphaDisp);
%fprintf(fid, 'TOP,4,6\n');

% Reporting results on .dat file
fprintf(fid, '*NODE PRINT,NSET=TOP\n');
fprintf(fid, 'RF\n');
%fprintf(fid, 'U\n');
%fprintf(fid, '*EL PRINT,ELSET=Eall\n');
%fprintf(fid, 'S\n');

% reporting resulst on .rfd file (for postprocessing in Calculix)
fprintf(fid, '*NODE FILE, OUTPUT=3D\n');
fprintf(fid, 'U\n');
fprintf(fid, '*EL FILE\n');
fprintf(fid, 'S\n');

fprintf(fid, '*END STEP\n');
fclose(fid);


% ====================== CALCULIX ANALYSIS ================================
% 6.2 - Solving the STATIC ANALYSIS
% =========================================================================
command_string = sprintf('%s\\ccx216.exe -i %s\\%s',FEMsolverPATH,pwd,fileinp(1:end-4));
[status1,result1] = system(command_string);

% ====================== CALCULIX ANALYSIS ================================
% 6.3 - Reading nodal rection forces at TOP nodes
% =========================================================================
filedat = sprintf('%s.dat',fileinp(1:end-4));
fid = fopen(filedat,'rt');
dta = fgets(fid);dta = fgets(fid);dta = fgets(fid);
data = fscanf(fid,'%f',[4,length(idx_max)]);
data = data';
fclose(fid);
RFtot = sum(data,1);

% ====================== CALCULIX ANALYSIS ================================
% 6.4 - Writing the .inp file for Calculix - BUCKLING ANALYSIS
% =========================================================================
fileinp2 = sprintf('.\\FEM_models\\%s_Buckling_model.inp',filenameinp);
fid = fopen(fileinp2,'wt');

% List of Nodes
fprintf(fid,'*NODE,NSET=Nall\n');
for n = 1:size(Nodes,1)
    fprintf(fid,'%i,%f,%f,%f\n',n,Nodes(n,2:4));
end
% List of Elements
fprintf(fid,'*ELEMENT,TYPE=B32,ELSET=Eall\n');
for e = 1:size(Elements,1)
    fprintf(fid,'%i,%i,%i,%i\n',e,Elements(e,[1,3,2]));
end
% Element sets (used afterward to specify the truss corss sections)
%%for e = 1:length(Diameters)
for e = 1:size(Elements,1)
    fprintf(fid, '*ELSET,ELSET=Elset%i\n',e);
    %%idx = find(idxTruss==e);
    %%fprintf(fid, '%i\n',idx);
    fprintf(fid, '%i\n',e);
end
fprintf(fid, '*NSET,NSET=BOTTOM\n');
fprintf(fid, '%i\n',idx_min);
% nodes at the top region
%idx_max = find (Nodes(:,4)>=(Zmin+0.99*DeltaZ));
fprintf(fid, '*NSET,NSET=TOP\n');
fprintf(fid, '%i\n',idx_max);

% Boundary Conditions
fprintf(fid,'*BOUNDARY\n');
fprintf(fid, 'BOTTOM,1,3\n');
%fprintf(fid, 'TOP,4,6\n');

% Setting the material
fprintf(fid, '*MATERIAL,NAME=BONE\n');
fprintf(fid, '*ELASTIC\n');
fprintf(fid, '19E3,.3\n');

% Beam sections
%%for e = 1:length(Diameters)
for e = 1:size(Elements,1)   
    fprintf(fid, '*BEAM SECTION, ELSET=Elset%i,MATERIAL=BONE,SECTION=CIRC\n',e);
    %%fprintf(fid, '%5.3f,%5.3f\n',Diameters(e),Diameters(e));
    % min = 1, max = 2, mean = 3
    fprintf(fid, '%5.3f,%5.3f\n',Truss_diameters(e,dia_type),Truss_diameters(e,dia_type));
    % calculating the orientation
    %%truss_dir = [1,0,10];
    V1 = Nodes(Elements(e,2),2:4)-Nodes(Elements(e,1),2:4);
    V2 = V1+(rand(1,3).*norm(V1));
    truss_dir = cross(V1,V2);
    truss_dir = truss_dir./norm(truss_dir); 
    fprintf(fid, '%4.3f,%4.3f,%4.3f\n',truss_dir);
end
% setting the analysis

fprintf(fid, '*STEP\n');
fprintf(fid, '*BUCKLE\n');
%fprintf(fid, '1,1e-2\n');
fprintf(fid, '1, 0.00005, 4, 100000\n');
fprintf(fid, '*CLOAD\n');
for n = 1:length(idx_max)
    fprintf(fid, '%i,%i,%f\n',data(n,1),1,data(n,2));
    fprintf(fid, '%i,%i,%f\n',data(n,1),2,data(n,3));
    fprintf(fid, '%i,%i,%f\n',data(n,1),3,data(n,4));
end

% Reporting results on .dat file
% fprintf(fid, '*NODE PRINT,NSET=TOP\n');
% fprintf(fid, 'RF\n');
% fprintf(fid, 'U\n');
% fprintf(fid, '*EL PRINT,ELSET=Eall\n');
% fprintf(fid, 'S\n');

% reporting resulst on .rfd file (for postprocessing in Calculix)
fprintf(fid, '*FILE FORMAT, ASCII\n');
fprintf(fid, '*NODE FILE\n');
fprintf(fid, 'U\n');
% fprintf(fid, '*EL FILE\n');
% fprintf(fid, 'S\n');

fprintf(fid, '*END STEP\n');
fclose(fid);

% ====================== CALCULIX ANALYSIS ================================
% 6.5 - Solving the Buckling analysis
% =========================================================================
command_string2 = sprintf('%s\\ccx216.exe -i %s\\%s',FEMsolverPATH,pwd,fileinp2(1:end-4));
[status2,result2] = system(command_string2);

% ====================== CALCULIX ANALYSIS ================================
% 6.6 - Reading The Buckling Load Multiplier
% =========================================================================
filedat = sprintf('%s.dat',fileinp2(1:end-4));
fid = fopen(filedat,'rt');
for r = 1:6
    dta = fgets(fid);
end
BucklingLF = fscanf(fid,'%f',[2,1]);

fclose(fid);

if isnan(RFtot(4))~=1
    fprintf('Si\tSf\tDz [mm]\tRtot [N]\tLPF\n');
    fprintf('%i\t%i\t%5.3f\t%5.3f\t%5.3f\n',...
        Si,Sf,DeltaZ*alphaDisp,RFtot(4),BucklingLF(2));
else
    fprintf('%i\t%i\t%5.3f\tn.a.\tn.a.\n',...
        Si,Sf,DeltaZ*alphaDisp);
end


if 1==10
% =========================================================================
% GRAPHICS
% =========================================================================
% PLOTTING A SLICE OF THE DATASETY
slice = randi([1,size(I,3)],1);
imagesc(I(:,:,slice))

% PLOTTING THE OUTCOME OF THE BINARISATION
h1 = figure;
col=[.7 .7 .8];
hiso = patch(isosurface(BW,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(BW,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(BW,hiso);
alpha(0.1);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
view(45,30);
hold on;

% mapping the distance between the skeleton and the bone surface
scatter3(y,x,z,D(IDskel)*20,D(IDskel));

%plot3(y,x,z,'square','Markersize',1,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(45,30)

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,0);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)
    wl = wl_new;
    skel2 = Graph2Skel3D(node2,link2,w,l,h);
    [A2,node2,link2] = Skel2Graph3D(skel2,0);
    wl_new = sum(cellfun('length',{node2.links}));
end

% display result
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end;
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='r'; % links are red
        end;
        if(node2(i).ep==1)
            col='k';
        end;
        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',5,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(45,30);
grid on; axis on;

% PLOTTING THE GRAPH ------------------------------------------------------
h2 = figure;
plot(G);%,'EdgeLabel',G.Edges.Weight)

% PLOTTING THE TRUSS STRUCTURE---------------------------------------------
h3 = figure;
hold on; axis equal;
view(45,30);
% plotting the mesh of beams
plot3(Nodes(:,2),Nodes(:,3),Nodes(:,4),'.c');
for n = 1:size(Elements,1);
    line([Nodes(Elements(n,1),2),Nodes(Elements(n,2),2)],[Nodes(Elements(n,1),3),Nodes(Elements(n,2),3)],...
    [Nodes(Elements(n,1),4),Nodes(Elements(n,2),4)],'color','m');
end

% h4 = figure;
% volshow(BW);
% h5 = figure;
% volshow(I);
end