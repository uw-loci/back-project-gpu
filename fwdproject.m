function fwdproject( rightbottomback, lefttopfront, dataset )
%Plots the forward projection
%   requires the right bottom back and left top front coordinates for the
%   projection area, the information in the makeDataset file, and the
%   target location

%make point list for the outline of projection area
x = [rightbottomback(1) lefttopfront(1) lefttopfront(1) rightbottomback(1) rightbottomback(1) rightbottomback(1)...
     rightbottomback(1) rightbottomback(1) rightbottomback(1) rightbottomback(1) lefttopfront(1) lefttopfront(1)...
     lefttopfront(1) lefttopfront(1) lefttopfront(1) lefttopfront(1) lefttopfront(1)...
     lefttopfront(1) rightbottomback(1)]';
y = [lefttopfront(2) lefttopfront(2) rightbottomback(2) rightbottomback(2) lefttopfront(2) lefttopfront(2)...
     rightbottomback(2) rightbottomback(2) lefttopfront(2) lefttopfront(2) lefttopfront(2) lefttopfront(2)...
     rightbottomback(2) rightbottomback(2) lefttopfront(2) lefttopfront(2) lefttopfront(2)...
     rightbottomback(2) rightbottomback(2)]';
z = [rightbottomback(3) rightbottomback(3) rightbottomback(3) rightbottomback(3) rightbottomback(3) lefttopfront(3)...
     lefttopfront(3) rightbottomback(3) rightbottomback(3) lefttopfront(3) lefttopfront(3) rightbottomback(3)...
     rightbottomback(3) lefttopfront(3) lefttopfront(3) rightbottomback(3) lefttopfront(3)...
     lefttopfront(3) lefttopfront(3)]';

%plot the scene
figure; hold on;
axis equal; 
title ('Forward Projection')
plot3(x, y, z, 'b-'); %Plot outline of projection area
plot3(dataset.laserOrigin(1), dataset.laserOrigin(2), dataset.laserOrigin(3), 'rx');%laser origin
plot3(dataset.cameraOrigin(1), dataset.cameraOrigin(2), dataset.cameraOrigin(3), 'kx');%Camera origin

plot3(dataset.cameraPos(:,1), dataset.cameraPos(:,2), dataset.cameraPos(:,3), 'k.'); %camera positions
plot3(dataset.laserPos(:,1), dataset.laserPos(:,2), dataset.laserPos(:,3), 'r.') %laser positions


plot3(dataset.targetpoints(:,1), dataset.targetpoints(:,2), dataset.targetpoints(:,3), 'mo') %Target position

end

