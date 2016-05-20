function [centroidx, centroidy]=extractPointfromRegion(camim)
projx=squeeze(sum(camim, 1));
projy=squeeze(sum(camim, 2));

centroidx=sum(projx.*(1:size(camim,2)))/sum(projx);
centroidy=sum(projy'.*(1:size(camim,1)))/sum(projy);
