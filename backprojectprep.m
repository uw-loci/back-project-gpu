function [ bpinput] = backprojectprep(dataset, pts)
%Prepares backprojection inputs to form needed by the cpp file
%   Detailed explanation goes here

bpinput.data = zeros(dataset.t*size(dataset.cameraPos,1), size(dataset.laserPos,1));


for i = 1:size(dataset.data, 1)
    for k = 1:size(dataset.data,2) 
        bpinput.data(((k-1)*size(dataset.data,3)+1):((k)*size(dataset.data,3)), i)= squeeze(dataset.data(i, k, :));
    end
end

bpinput.pts = pts'; %form needed for c++, [x; y; z]
bpinput.laserPos = dataset.laserPos';
bpinput.laserNorm = dataset.laserNorm';
bpinput.cameraPos = dataset.cameraPos';
bpinput.cameraNorm = dataset.cameraNorm';
bpinput.laserOrigin = dataset.laserOrigin';
bpinput.cameraOrigin = dataset.cameraOrigin';

end

