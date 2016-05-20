function [pts, pixels] = makegrid( rightbottomback, lefttopfront, gridsize )
%Creates the grid defining the backprojection area
%   distances are in meters


pixels=ceil(abs(rightbottomback-lefttopfront)./gridsize);
vx=[1 0 0].*gridsize;
vy=[0 1 0].*gridsize;
vz=[0 0 1].*gridsize;
npixels=pixels(1)*pixels(2)*pixels(3);
pts=zeros(npixels,3);
j=1;
for x=1:pixels(1)
    for y=1:pixels(2)
        for z=1:pixels(3)
            pts(j,:)=rightbottomback+(x-1).*vx+(y-1)*vy+(z-1)*vz;
            j=j+1;
        end
    end
end
end

