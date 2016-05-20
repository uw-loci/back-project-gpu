function dataset = loadDirsigDataset(metadata, data_path, dirsig_path)
%packages data collected in dirsig simulations for backprojection
%   Detailed explanation goes here

    addpath([dirsig_path 'extra' filesep 'matlab-lidar' filesep]);
    
    sl=metadata.spotlist;
    
    % load laser positions
    f=fopen([data_path metadata.truth_basename '_' num2str(sl(1,1)) '_' num2str(sl(1,2)) '.img'], 'r');
    [data_path metadata.truth_basename '_' num2str(sl(1,1)) '_' num2str(sl(1,2)) '.img']
    truth_cube=fread(f, 'double');
    fclose(f);
    
    lidat=readProtoLidar([data_path metadata.data_basename '_' num2str(sl(1,1)) '_'  num2str(sl(1,2)) '.bin']);
    tTotal = size(lidat.task.pulse.signal,1);
    truth=reshape(truth_cube, [], size(lidat.task.pulse.signal,2), size(lidat.task.pulse.signal,3));
    
    dataset.laserOrigin=lidat.task.pulse.header.platformLocation';
    dataset.cameraOrigin=lidat.task.pulse.header.platformLocation';
    
    disp('Computing normals ...')
    dataset.cameraPos=zeros(size(lidat.task.pulse.signal,2)*size(lidat.task.pulse.signal,3), 3);
    dataset.cameraNorm=zeros(size(lidat.task.pulse.signal,2)*size(lidat.task.pulse.signal,3), 3);
    gradientx=squeeze(truth(12,1:(end-1),:)-truth(12,2:end,:));
    gradienty=squeeze(truth(12,:,1:(end-1))-truth(12,:,2:end));
    normals=zeros(3, size(gradientx,1)+1, size(gradientx,2));
    normals(1,1:(end-1),:)=gradientx;
    normals(2,:,1:(end-1))=gradienty;
    normals(3,:,:)=1;
    normals(:,end,:)=normals(:,end-1,:);
    normals(:,:,end)=normals(:,:,end-1);
    for j=1:(size(lidat.task.pulse.signal,2)*size(lidat.task.pulse.signal,3))
        dataset.cameraPos(j,:)=[truth(10, mod(j, 64)+1, floor((j-1)/64)+1) truth(11, mod(j, 64)+1, floor((j-1)/64)+1) truth(12, mod(j, 64)+1, floor((j-1)/64)+1)];
        dataset.cameraNorm(j,:)=[normals(1, mod(j, 64)+1, floor((j-1)/64)+1) normals(2, mod(j, 64)+1, floor((j-1)/64)+1) normals(3, mod(j, 64)+1, floor((j-1)/64)+1)];
        dataset.cameraNorm(j,:)=dataset.cameraNorm(j,:)./sqrt(sum(dataset.cameraNorm(j,:).^2));
    end
    dataset.laserPos=zeros(size(sl,1),3);
    dataset.laserNorm=zeros(size(sl,1),3);
    for i=1:size(sl,1)
        dataset.laserNorm(i,1)=normals(1, ceil(sl(i,1)./4), ceil(sl(i,2)./4));
        dataset.laserNorm(i,2)=normals(2, ceil(sl(i,1)./4), ceil(sl(i,2)./4));
        dataset.laserNorm(i,3)=normals(3, ceil(sl(i,1)./4), ceil(sl(i,2)./4));
        dataset.laserNorm(i,:)=dataset.laserNorm(i,:)./sqrt(sum(dataset.laserNorm(i,:).^2));
    end
    if isfield(metadata, 'dataLimits') && (numel(metadata.dataLimits)==2)
        dataset.t=metadata.dataLimits(2)-metadata.dataLimits(1)+1;
    else
        dataset.t=size(lidat.task.pulse.signal,1);
        metadata.dataLimits(1)=1;
        metadata.dataLimits(2)=dataset.t;
    end
    
    disp(['Setting data limits to [' num2str(metadata.dataLimits(1)) ' ' num2str(metadata.dataLimits(2)) '] ...'])
    dataset.data=zeros(size(sl,1), size(lidat.task.pulse.signal,2)*size(lidat.task.pulse.signal,3), dataset.t);
   
    if exist([data_path 'noise.mat'], 'file')==0
       dataset.data = ones(size(dataset.data));
           
    else        
        for i=1:size(sl,1)
            lidat=readProtoLidar([data_path metadata.data_basename '_' num2str(sl(i,1)) '_'  num2str(sl(i,2)) '.bin']);
            disp(['Loading ' metadata.data_basename '_' num2str(sl(i,1)) '_'  num2str(sl(i,2)) '.bin'])
            dataset.laserPos(i,:)=[squeeze(truth(10,1,1)) squeeze(truth(11,1,1)) 0]+[squeeze(truth(10,end,end)-truth(10,1,1))/256*sl(i,1) squeeze(truth(11,end,end)-truth(11,1,1))/256*sl(i,2) 0]+truth(12,ceil(sl(i,1)/4),ceil(sl(i,2)/4))*[0 0 1];%128 128
    %         dataset.laserPos(i,:)=[squeeze(truth(10,1,1)) squeeze(truth(11,1,1)) 0]+[squeeze(truth(10,end,end)-truth(10,1,1))/256*sl(i,2) squeeze(truth(11,end,end)-truth(11,1,1))/256*sl(i,1) 0]+truth(12,ceil(sl(i,2)/4),ceil(sl(i,1)/4))*[0 0 1];%128 128
            for j=1:(size(lidat.task.pulse.signal,2)*size(lidat.task.pulse.signal,3))
                dataset.data(i,j,:)=lidat.task.pulse.signal(metadata.dataLimits(1):metadata.dataLimits(2), mod(j, 64)+1, floor((j-1)/64)+1);
            end
        end
    end
%     dataset.targetpoints=[];
%     dataset.t=size(lidat.task.pulse.signal,1);
    dataset.deltat=(lidat.task.pulse.header.timeGateStop-lidat.task.pulse.header.timeGateStart)/tTotal*2.99e8;
    dataset.t0=lidat.task.pulse.header.timeGateStart*2.99e8+(metadata.dataLimits(1)-1)*dataset.deltat;
    
    %tshiftd=gateStart*2.99e8*100+450*tpp+2301*tpp;%+470*tpp;
    
end

