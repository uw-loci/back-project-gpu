function [ dataset ] = createTarget( dataset, target, normal, area, camarea, camapp, tr, ts_x, ts_y, d, uv_x, uv_y )
%Populates the dataset.data for the target object
%   All distances are in meteres
%   Cooridnate Origin = at the galvos (time and laser origin)


%Create a list of target points and simulate the data
dataset.targetpoints = zeros(ts_x.*ts_y,3);
i=1;
for tgx=1:ts_x
    for tgy=1:ts_y
        dataset.targetpoints(i,:)=target+uv_x.*d.*((tgx-1)/1)+uv_y.*d.*((tgy-1)/1);
        dataset=createPoints(dataset, dataset.targetpoints(i,:), normal, area, camarea, camapp, tr);
        i=i+1;
    end
end


end

