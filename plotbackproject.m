function plotbackproject( rightbottomback, lefttopfront, dataset, Wraw, pts, threshold )
%Plots the forward projection
%   requires the right bottom back and left top front coordinates for the
%   projection area, the information in the makeDataset file, data from the
%   backprojection, and the coordinates

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
title ('Back Projection')
plot3(x, y, z, 'b-'); %Plot outline of projection area

if dataset.laserOrigin(3) < 100
    plot3(dataset.laserOrigin(1), dataset.laserOrigin(2), dataset.laserOrigin(3), 'rx');%laser origin
end

if dataset.cameraOrigin(3) < 100
    plot3(dataset.cameraOrigin(1), dataset.cameraOrigin(2), dataset.cameraOrigin(3), 'kx');%Camera origin
end

plot3(dataset.cameraPos(:,1), dataset.cameraPos(:,2), dataset.cameraPos(:,3), 'k.'); %camera positions
plot3(dataset.laserPos(:,1), dataset.laserPos(:,2), dataset.laserPos(:,3), 'r.') %laser positions

%Plotting normalized results from the backprojection

%change how plotting works depending on if W or Wraw are passed
if size(size(Wraw),2)==2
    for i=1:size(Wraw,1)
        if Wraw(i) > threshold
           plot3(pts(i,1), pts(i,2), pts(i,3), 'm.')
        end
    end
    
else
   for i=1:size(Wraw,1)
       for j=1:size(Wraw,2)
           for k=1:size(Wraw,3)
               if Wraw(i,j,k) > threshold 
                  plot3(pts(i,j,k,1), pts(i,j,k,2), pts(i,j,k,3), 'm.')
               end
           end
       end
   end
 title 'Backproject Ray Filter' 
end


a = isfield(dataset, 'targetpoints');

if a == 1
    plot3(dataset.targetpoints(:,1), dataset.targetpoints(:,2), dataset.targetpoints(:,3), 'bo') %Target position
end
end

