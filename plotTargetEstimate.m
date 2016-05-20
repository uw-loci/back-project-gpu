function [ texpected, d ] = plotTargetEstimate( dataset, target, lp )
%time plot of the experimental data
%   Detailed explanation goes here

texpected = zeros(size(dataset.data,2),1); %time bin we expect to see the target
for i = 1: size(dataset.data,2)
    d = sqrt(sum((dataset.laserOrigin-dataset.laserPos(lp,:)).^2));
    d = d+sqrt(sum((dataset.laserPos(lp,:)-target).^2));
    d = d+sqrt(sum((dataset.cameraPos(i,:)-target).^2));
    d = d+sqrt(sum((dataset.cameraPos(i,:)-dataset.cameraOrigin).^2));
    texpected(i,1)=(d-dataset.t0)*(1/dataset.deltat);
end


figure; hold on;
colormap(gray);
imagesc(squeeze(dataset.data(lp,:,:)));
for i=1:size(texpected,1)
    plot(texpected(i), i, 'ro')
end
end

