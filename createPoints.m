function [ dataset ] = createPoints(dataset, target, normal, area, camarea, camapp, tr)
%creates data for points in the target object
%t0 is the distance that corresponds with index 0.  That is not the first
%element of the data array (that is index = 1)
    
c = 2997924598; %m/s 
numpixel = tr/dataset.deltat*c; %number of pixels at full width half max
sigma = numpixel/2*sqrt(2*log(2)); %for the gaussian

for lp = 1:size(dataset.laserPos,1)
    dlolp = sqrt(sum((dataset.laserPos(lp,:,:)-dataset.laserOrigin).^2));
    dlpt = sqrt(sum((target-dataset.laserPos(lp,:,:)).^2));
    uvlpt = (target-dataset.laserPos(lp,:,:))./dlpt; %unit vector of laserposition-target
    
    for cp = 1:size(dataset.cameraPos, 1)
        dtcp = sqrt(sum((dataset.cameraPos(cp,:,:)-target).^2));
        dcpco = sqrt(sum((dataset.cameraOrigin-dataset.cameraPos(cp,:,:)).^2));
        
        distfactor = area/(2*pi*dlpt^2)*camarea/(2*pi*dtcp^2)*camapp/(2*pi*dcpco^2);
        cosfactor = dot(-normal, uvlpt)*dot(normal, -uvlpt);
        
        dist = dlolp+dlpt+dtcp+dcpco; 
        dist = dist-dataset.t0;
        index = ceil(double(dist)/double(dataset.deltat));
         if (index > 0) && (index < dataset.t)
             dataset.data(lp, cp, :) = squeeze(dataset.data(lp, cp, :)) + ...
             squeeze(distfactor*cosfactor* exp(-1/2.*(((1:dataset.t)-index)/sigma).^2))';
         end
         
    end
    
end

end

