clear all;
close all;

load testvol

skel = Skeleton3D(testvol);

figure();
col=[.7 .7 .8];
hiso = patch(isosurface(testvol,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(testvol,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(testvol,hiso);
alpha(0.2);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'square','Markersize',0.1,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(140,80)


% mapping the distance
% https://uk.mathworks.com/help/images/distance-transform-of-a-binary-image.html
% https://uk.mathworks.com/help/images/ref/bwdist.html
C = -1.*(testvol-1);
D = bwdist(C);  % default = Eucledean distance
idx = find(skel);
skelD = zeros(size(skel));
skelD(idx)=D(idx);
scatter3(y,x,z,D(idx),D(idx));
A1 = skel(:,:,13);
A2 = C(:,:,13);
A3 = D(:,:,13);
h2 = figure;
subplot(1,3,1);
imagesc(A1);
subplot(1,3,2);
imagesc(A2);
subplot(1,3,3);
imagesc(A3);