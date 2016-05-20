% Create the target and Set Parameters for the case where no data is
% available

c = 2997924598; %m/s

% Dataset
laserPosFBB = [0,0,-2]; 
laserPosRTB = [1, 1, -2];
lsgridsize = 0.25; %distance between laser points

dataset.cameraPos = [0.5,.5,-2; 1,.5,-2];
dataset.laserOrigin = [0,0,0];
dataset.cameraOrigin = [0.75,0,0];
dataset.t = 1000; %Number of time points in the data
dataset.t0 = 4; %m Distance pulse traveled between the laser origin and the detector
dataset.deltat = 5e-12 * c; %m
laserNorm = [0,0,1];
cameraNorm = [0,0,1];

%Make grid defining the projection area
rightbottomback=[0,0.25,-1.5];
lefttopfront=[1,0.75,-0.5];
gridsize=.01; %distance between points

% Target
target = [0.3, 0.5, -1]; %location of the object 
normal = [0, 0, 1]; %ray defining the direction the object is facing
area = .01; %size of target
camarea = .01; %area of camera looking at the wall
camapp = .5^2*pi; %camera aperature area
tr = 30e-12; %time resolution (s)

% target = [0.5, 0.5, -1.5]; %location of the object 
% normal = [0, 0, 1]; %ray defining the direction the object is facing
% area = .01; %size of target
% camarea = .01; %area of camera looking at the wall
% camapp = 0.0125^2*pi; %camera aperature area
% tr = 30e-12; %time resolution (s)

%Define number of points (x and y) in target
ts_x = 2; 
ts_y = 1;

%Define distance between points (m)
d = 0.3;

%unit vectors defining orientation of the target
uv_x = [1 0 0];
% uv_x = [cos(pi/6) 0 -sin(pi/6)];
uv_y = [0 1 0];
% uv_y = [0 cos(pi/6) sin(pi/6)] %angled unit vector

%% Making the dataset
lpts = (abs(laserPosRTB-laserPosFBB)./lsgridsize)+[1 1 1];
nlpts = lpts(1).*lpts(2).*lpts(3);

dataset.cameraNorm = zeros(size(dataset.cameraPos));
for i = 1:size(dataset.cameraPos,1)
    dataset.cameraNorm(i,:)=cameraNorm;
end

dataset.laserNorm = zeros(nlpts, 3);
dataset.laserPos = zeros(nlpts,3);
  vx=[1 0 0].*lsgridsize;
  vy=[0 1 0].*lsgridsize;
 
 j=1;
 for x=1:lpts(1)
     for y=1:lpts(2)
             dataset.laserPos(j,:)=laserPosFBB+(x-1).*vx+(y-1)*vy;
             dataset.laserNorm(j,:)=laserNorm;
             j=j+1;
     end
 end
 
dataset.data = zeros(size(dataset.laserPos,1), size(dataset.cameraPos,1), dataset.t);
