clearvars;
addpath('Utils','-end');
addpath('Utils/quadfit','-end');  
addpath('Utils/geom3d/geom3d/','-end'); 

%% Inputs
param.ptCloudname='';% 'Mangaweka.ply' 'Otira_1cm_grains.ply' 'Test1_clean_registered.ply'
if isempty(param.ptCloudname)==1;
    [param.ptCloudname,param.ptCloudpathname] = uigetfile('*.ply','Select the *.ply point cloud file');
end

%% Loading data
[ptCloud,param]=loadptCloud(param);

%% Algorithm parameters - Compute point cloud size and scaling of the algorithm
% Load parameters
param=defineparameters(ptCloud,param);
% Display algorithm header
displayheader(ptCloud,param);

%% Segment and cluster the point cloud into a point cloud of potential grains
% Find neighbors of point cloud 
[indNeighbors,D]=knnsearch(ptCloud.Location,ptCloud.Location,'K',param.nnptCloud+1);
indNeighbors=indNeighbors(:,2:end);
D=D(:,2:end);
% determine node surface
surface=pi.*min(D,[],2).^2; 
% Compute normals and force them to point towards positive Z
normals = pcnormals(ptCloud,param.nnptCloud);
[normals]=adjustnormals3d(ptCloud.Location(:, 1),ptCloud.Location(:, 2),ptCloud.Location(:, 3),normals,[mean(ptCloud.Location(:,1)), mean(ptCloud.Location(:,2)),10000]);

%% Initial segmentation with Fastscape
[labels,nlabels,labelsnpoint,stack,nstack,ndon,isink]=segment_labels(ptCloud,param,indNeighbors);
cmaplabels=rand(nlabels,3);
    % Plot
    if param.iplot==1;
        pcshow(ptCloud.Location,labels);
        colormap(cmaplabels);
        hold on;
        set(gcf,'color','w');
        set(gca,'color','w');
        axis equal tight;
        hold on;
        plot3(ptCloud.Location(isink,1),ptCloud.Location(isink,2),ptCloud.Location(isink,3),'.r');
        axis off;
    end
    
%% Generate a Pebble structure
for i=1:nlabels;
    ind=find(labels==i);
    Pebble(i).Location=ptCloud.Location(ind,:);
    Pebble(i).ind=ind;
    Pebble(i).surface=surface(ind);
end 

%% Fitting ellipsoids
% fit
[Ellipsoidm]=fitellipsoidtograins(Pebble,param,nlabels);

%% Plot cloud
pcshow(ptCloud.Location,labels);
        colormap(cmaplabels);
        hold on;
        set(gcf,'color','w');
        set(gca,'color','w');
        axis equal tight;
        hold on;
        plot3(ptCloud.Location(isink,1),ptCloud.Location(isink,2),ptCloud.Location(isink,3),'.r');
        axis off;

%% Plot

j = 1
plot_ellipsoid_im(Ellipsoidm(j).p,'EdgeColor',cmaplabels(j,:));

cb = colorbar('north');
set(cb,'position',[.5 .75 .1 .02]);
ylabel(cb,'Labels');

%% Plot python ellipsoid

python_p = [ 41.93209985, 144.26738821, 151.64467324, -94.85891829, ...
    16.70173584, 18.01863136, 530.51530604, -1659.08193363, ...
    -157.51707559, 4774.85246949];
plot_ellipsoid_im(python_p);

%%
python_center = [ 0.19271343,  0.01573514, -0.03469648];
python_radii = [0.46686707, 0.17188535, 0.18014326];
python_rotation = [[-0.9263723 , -0.36679569,  0.0854124 ]
    [-0.30781516,  0.86809616,  0.38943406]
    [ 0.21698891, -0.33446969,  0.91708551]];

cx = python_center(1);
cy = python_center(2);
cz = python_center(3);
ap = python_radii(1);
bp = python_radii(2);
cp = python_radii(3);
R = python_rotation;

plot_ellipsoid(cx, cy, cz, ap, bp, cp, R);





