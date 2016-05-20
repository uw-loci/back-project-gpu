function [ W ] = rayFilter( W, threshold )
%filter along the rays going away from the wall
%   finds the maximum value in each row.  If it is less than the threshold,
%   it sets equal to zero, if greater than the threshold, it keeps it.

%***make sure that the rays are drawn perpendicular to the surface light is
%bouncing off (ie wall)


%assume there can only be one point in the direction perpendicular to the
%plane of the wall

for j = 1:size(W,2)
    for i = 1:size(W,3) 
        b = W(:,j, i); %written as (z, y, x)
		[val, index] = max(b);
        if val<threshold
            b(:,1)=0;
        else
            b(:,1)=0;
            b(index,1)=val;
        end
		W(:, j, i)=b;
    end
end


%doing a ray filter in the other directions, set the value to zero if it is
%less than 90% of the max
for k = 1:size(W,1)
    for i = 1:size(W,3) 
        b = W(k,:, i); %written as (z, y, x)
		val = max(b);
        for l = 1:size(b,2)
            if b(1,l)<=(0.9*val)
                b(1,l)=0;
            end
            W(k, :, i)=b;
   
        end
    end
end

for k = 1:size(W,1)
    for j = 1:size(W,2) 
        b = squeeze(W(k,j, :)); %written as (z, y, x)
		val = max(b);
        for l = 1:size(b,1)
            if b(l,1)<=(0.85*val)
                b(l,1)=0;
            end
            W(k, j, :)=b;
   
        end
    end
end


end

