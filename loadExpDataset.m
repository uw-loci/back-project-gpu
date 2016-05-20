function [ dataset ] = loadExpDataset( data_path, pic_path, cam_Cal, parameters, dataset )
%read and load the data from the experimental setup
%   Detailed explanation goes here

%Load data from files in data_path
display('loading data...')

dataUnpadded = dlmread([data_path 'data.txt']);
[firstbounce, ~, ~] = TCSPCread([data_path 'first_bounce.dat']);

peaks = round(size(firstbounce, 1)/parameters.period); % # of peaks displayed on hydraharp at a time
nlaserPos = size(dataUnpadded, 1); % # of laser positions
% nlaserPos = 1; % # of laser positions
ncameraPos = size(dataset.cameraPos, 1); %# fo camera positions

figure; plot(firstbounce); title 'First Bounce';
hold on; plot(parameters.fbduration,0:70, 'r');
legend('first bounce','fb duration')

%return

%Populate variables
data = zeros(size(dataUnpadded, 1), size(dataUnpadded,2)+10000);
data(:,1:size(dataUnpadded, 2)) = dataUnpadded;

fb=zeros(peaks, parameters.fbduration); %first bounce
for i = 1:peaks
    fb(i,:) = firstbounce(((i-1)*parameters.period+1):((i-1)*parameters.period + parameters.fbduration));
end

figure; hold on; title 'fb'
for i = 1:peaks
        if i == 1
            plot(fb(1,:), 'r')
        end
        if i == 2
            plot(fb(2,:), 'b') 
        end
        if i == 3
            plot(fb(3,:), 'm')
        end
        if i == 4
            plot(fb(4, :), 'g')
        end
end
legend('1', '2', '3', '4')
hold off; 

[~, pos]=max(sum(fb,1)); %adds the peaks displayed at the same time to minimize noise, finds the peak
fbIndex = pos;
offset = fbIndex + ceil(parameters.offset/dataset.deltat); %offset in time bins
% offset = 1300;
dataset.t0 = sqrt(sum((dataset.laserOrigin-dataset.cameraPos).^2))+sqrt(sum((dataset.cameraOrigin-dataset.cameraPos).^2)) + parameters.offset; 

dataset.data = zeros(nlaserPos, ncameraPos, dataset.t);

for i = 1:ncameraPos
for j = 1:nlaserPos
    for k = 1:size(fb,1)
        dataset.data(j,i,:) = dataset.data(j,:) + data(j,((offset + (k-1)*parameters.period)...
            :(offset +(k-1)*parameters.period+dataset.t-1)));
    end
end
end
% dataset.data = zeros(nlaserPos, dataset.t);
% for j = 1:nlaserPos
%     for k = 1:size(fb,1)
%         dataset.data(j,:) = dataset.data(j,:) + data(j,((offset + (k-1)*parameters.period)...
%             :(offset +(k-1)*parameters.period+dataset.t-1)));
%     end
% end


%--------------------not edited yet--- from loadapdCalCam file

laserPos = zeros(2, nlaserPos);
for i = 1:nlaserPos
    camImg = imread([pic_path num2str(i-1) '.jpg']);
    camImg(:,610:end,:)=0;
    camImg(460:end,320:end,:)=0;
    camImg=camImg(end:-1:1,:,:);
    
    imgThresh = squeeze(camImg(:,:,2));
    imgThresh(imgThresh<0.95*max(imgThresh(:))) = 0;
    projx = squeeze(sum(imgThresh, 1));
    projy = squeeze(sum(imgThresh, 2));
    centroidx = sum(projx.*(1:size(camImg, 2)))/sum(projx);
    centroidy=sum(projy'.*(1:size(camImg,1)))/sum(projy);
    
    laserPos(1,i) = centroidx;
    laserPos(2,i) = centroidy; 
end    

camcal=generateCameraCalibration(cam_Cal);

worldpoints=zeros(2, nlaserPos);
for i=1:nlaserPos
    worldpoints([2 1],i)=wcfromccv2(camcal.world, camcal.cam, laserPos([2 1],i));
end

origincam=[504.1042 480-436.8251];
originworld([2 1])=wcfromccv2(camcal.world, camcal.cam, origincam([2 1]));

worldpoints=worldpoints([2 1],:);
worldpoints(2,:)=worldpoints(2,:);
originworld=originworld([2 1]);
originworld(2)=originworld(2);

figure; plot(squeeze(dataset.data)')
figure; hold on; plot(laserPos(1,:), laserPos(2,:), 'g.')
plot(origincam(1), origincam(2), 'b.'); hold off;
figure; hold on; plot(worldpoints(2,:), worldpoints(1,:), 'r.')
plot(originworld(2), originworld(1), 'b.'); hold off; 

worldpoints(1,:)=worldpoints(1,:)-originworld(1);
worldpoints(2,:)=worldpoints(2,:)-originworld(2);

figure; plot(worldpoints(2,:), worldpoints(1,:), 'r.')

dataset.laserPos = zeros(nlaserPos, 3);
dataset.laserNorm = dataset.laserPos;
dataset.laserNorm(:,3) = 1;

for i=1:nlaserPos
    dataset.laserPos(i, :)=[worldpoints(2,i) worldpoints(1,i) 0].*100;
end

dataset.cameraNorm = zeros(ncameraPos, 3);
dataset.cameraNorm(:,3) = 1;
end

