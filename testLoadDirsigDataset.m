close all
clear all

%some metadata describing the dataset
metadata.data_basename='lidar_cave_jess';
metadata.truth_basename='lidar_pathlength';
metadata.dataLimits=[1500 4000];
sl=load('/home/andreas/Documents/Dropbox/dirsig_scripted/spotlist.mat');
%this is the list of laser spots used in the makepatternsLoop script to
%create the files for the simulations
metadata.spotlist=sl.spotlist;

%the path to the .bn and .img files with the actual data
data_path='/home/andreas/Documents/Dropbox/research/projects/periscope/dirsig_results/20140123FancyCaveWBlocks/';

%the path to the dirsig installation
dirsig_path='/home/andreas/Applications/dirsig/';

ds=loadDirsigDataset(metadata, data_path, dirsig_path);

figure; hold on;
title('Camera Positions')
plot3(ds.cameraPos(:,1), ds.cameraPos(:,2), ds.cameraPos(:,3), 'r.');
for i=1:size(ds.cameraNorm, 1)
    plot3([ds.cameraPos(i,1) ds.cameraPos(i,1)+ds.cameraNorm(i,1)], [ds.cameraPos(i,2) ds.cameraPos(i,2)+ds.cameraNorm(i,2)], [ds.cameraPos(i,3) ds.cameraPos(i,3)+ds.cameraNorm(i,3)], 'g')
end
hold off;

figure; hold on;
title('Laser Positions')
plot3(ds.laserPos(:,1), ds.laserPos(:,2), ds.laserPos(:,3), 'b.');
for i=1:size(ds.laserNorm, 1)
    plot3([ds.laserPos(i,1) ds.laserPos(i,1)+ds.laserNorm(i,1)], [ds.laserPos(i,2) ds.laserPos(i,2)+ds.laserNorm(i,2)], [ds.laserPos(i,3) ds.laserPos(i,3)+ds.laserNorm(i,3)], 'g')
end
hold off;

ds.targetpoints=[0 0 0]; %fwdproject complains if there are no targetspots. maybe that should be fixed.

fwdproject( [0 0 0], [100 100 100], ds);

figure; imagesc(squeeze(ds.data(13,:,:)))