clear all
close all
clc

c = 2.99e8; %m/s

dataSource = 'none';
psf = false; 
%{ 
   for dataSource
      'dirsig' = data from dirsig
      'experimental' = data from experimental setup
      'none' = no data, simulate the target
    
    psf: get the pointspread function
    
%}


%% Define Variables and File Locations

%For Filtered Backprojection
threshold = .4;
iterations = 1;


%Loading Data
switch dataSource;
    case 'dirsig'
        setParametersDirsig;
        
    case 'experimental'
        setParametersExperimental;
        
    case 'none'
        setParametersNone;
end
%% load PSF Parameters
if  psf==true 
    setParametersPsf; %Write this file    
end


%%  Run main code

switch dataSource
    case 'dirsig'
        % makes the dataset for dirsig resutls
        dataset = loadDirsigDataset(metadata, data_path, dirsig_path);
        dataset.t0 = dataset.t0+48; %fudge factor

        % threshold the data before backprojection
%         dataset.data(dataset.data < 0.5*max(dataset.data(:))) = 0;
        
        % Create projection area 
        [pts, pixels] = makegrid(rightbottomback, lefttopfront, gridsize);
        
        % Plot the scene
        plotScene( rightbottomback, lefttopfront, dataset);

        figure; imagesc(squeeze(dataset.data(13,:,:)))
%         for i = 1:size(dataset.laserPos,1)
        [texpected, dist] = plotTargetEstimate(dataset, target_est, 13);
%         end      

    case 'experimental'
        dataset = loadExpDataset(data_path, pic_path, cam_Cal, parameters, dataset);
                 
        % Create projection area 
        [pts, pixels] = makegrid(rightbottomback, lefttopfront, gridsize);
        
        % Plot the scene
        plotScene( rightbottomback, lefttopfront, dataset);

        
    case 'none'
       % Create the forward projection of the target
         dataset = createTarget(dataset, target, normal, area, camarea,...
                   camapp, tr, ts_x, ts_y, d, uv_x, uv_y );
        
       % Create projection area 
         [pts, pixels] = makegrid(rightbottomback, lefttopfront, gridsize);
        
       % Plot fwd projection
         fwdproject( rightbottomback, lefttopfront, dataset );
end

%% PSF
 if psf == true
     
      targetPsf = (rightbottomback+lefttopfront)./2;

    % Create forward projection of the pointspread function
      dataset_psf = createTarget(dataset, targetPsf, normal, area,...
                    camarea, camapp, tr, 1, 1, d, uv_x, uv_y );
    % Transpose dataset for the backprojection            
      bpinput_psf = backprojectprep(dataset_psf, pts);  
        
    % Running backprojection for the pointspread function
      Wraw_psf = backproject_fast(double(bpinput_psf.data), double(bpinput_psf.laserPos), ...
                 double(bpinput_psf.pts), double(bpinput_psf.cameraPos), ...
                 double(bpinput_psf.laserOrigin), double(bpinput_psf.cameraOrigin),...
                 double(dataset_psf.deltat), double(dataset_psf.t0), ...
                 double(bpinput.laserNorm), double(bpinput.cameraNorm));
                    
    if iterations > 1                   
        for i = 2:iterations
            Wraw_psf = backproject_fast_iterate(double(bpinput_psf.data),...
                       double(bpinput_psf.laserPos), double(bpinput_psf.pts),... 
                       double(bpinput_psf.cameraPos), double(bpinput_psf.laserOrigin),...
                       double(bpinput_psf.cameraOrigin), double(dataset_psf.deltat), ...
                       double(dataset_psf.t0), Wraw_psf, ...
                       double(bpinput.laserNorm), double(bpinput.cameraNorm));
         end
    end
        
    % normalize the raw data for psf
    Wraw_psf=Wraw_psf-min(Wraw_psf(:));
    Wraw_psf=Wraw_psf./max(Wraw_psf(:));

    px=pixels(1);
    py=pixels(2);
    pz=pixels(3);

    ptsReshaped=reshape(pts, [pz, py, px,3]);
    W_psf = reshape(Wraw_psf,[pz, py, px]);

 end

dataset.data=dataset.data+rand(size(dataset.data))*0.1*max(dataset.data(:));
%% Backproject
%Transpose dataset files for backprojection
bpinput = backprojectprep(dataset, pts); 
 
%call backprojection for Target
    %time per pixel and timeshift are in meters   
Wraw = backproject_fast(double(bpinput.data), double(bpinput.laserPos), ...
                        double(bpinput.pts), double(bpinput.cameraPos), ...
                        double(bpinput.laserOrigin), double(bpinput.cameraOrigin),...
                        double(dataset.deltat), double(dataset.t0), ...
                        double(bpinput.laserNorm), double(bpinput.cameraNorm));

% Wraw = backproject_fastOld(double(bpinput.data), double(bpinput.laserPos), ...
%                         double(bpinput.pts), double(bpinput.cameraPos), ...
%                         double(bpinput.laserOrigin), double(bpinput.cameraOrigin),...
%                         double(dataset.deltat), double(dataset.t0));

display('reshaping first itteration...')
px=pixels(1);
py=pixels(2);
pz=pixels(3);
W = reshape(Wraw,[pz, py, px]);
slice=squeeze(W(:,ceil(size(W,2)/2),:));
figure; imagesc(slice)
slice=slice./max(slice(:));
imwrite(slice.*64, jet, ['iteration' num2str(i) '.tif'], 'tif')
% figure; 
% for i=1:size(W,2)
%     subplot(4,5,i);
%     imagesc(squeeze(W(:,i,:)));
%     if i == 3
%       title ('Number of iterations: 1')
%     end
% end
                    
if iterations > 1                   
   for i = 2:iterations
     Wraw = backproject_fast_iterate(double(bpinput.data), double(bpinput.laserPos), ...
                        double(bpinput.pts), double(bpinput.cameraPos), ...
                        double(bpinput.laserOrigin), double(bpinput.cameraOrigin),...
                        double(dataset.deltat), double(dataset.t0), Wraw, ...
                        double(bpinput.laserNorm), double(bpinput.cameraNorm));
                    

    W = reshape(Wraw,[pz, py, px]);
    slice=squeeze(W(:,ceil(size(W,2)/2),:));
    figure; imagesc(slice)
    slice=slice./max(slice(:));
    imwrite(slice.*64, jet, ['iteration' num2str(i) '.tif'], 'tif')
%     figure; 
%     for j=1:size(W,2)
%         subplot(4,5,j);
%         imagesc(squeeze(W(:,j,:)));
%         if j == 3
%           title (['Number of iterations: ' num2str(i)])
%         end
%     end
   end
end   
return
%normalize the raw data
display('normalize raw data')

% if exist([data_path 'noise.mat'], 'file')==0 && strcmp(dataSource, 'dirsig')==1
%     save([data_path 'noise.mat'], 'Wraw');
%     return
%     
% elseif strcmp(dataSource, 'dirsig')==1
%     Elipses = load([data_path 'noise.mat']);
%     Elipses.Wraw(Elipses.Wraw <= 2) = 1;
% %     Wraw = Wraw./Elipses.Wraw;
% end

Wraw=Wraw-min(Wraw(:));
Wraw=Wraw./max(Wraw(:));

plotbackproject( rightbottomback, lefttopfront, dataset, Wraw, pts, threshold )

return 
%images of raw data   

px=pixels(1);
py=pixels(2);
pz=pixels(3);

ptsReshaped=reshape(pts, [pz, py, px,3]);
W = reshape(Wraw,[pz, py, px]);


% W(:,:,1:(end-9)) = W(:,:,10:end);
% W(:,:,(end-10):end)=0;
% figure; imagesc(squeeze(W(:,:, 2)));
% figure; 
% for i=1:size(W,1)
%     subplot(2,5,i);
%     imagesc(squeeze(W(i,:,:)));
%     if i == 3
%       title (['Number of iterations: ' num2str(iterations)])
%     end
% end

% figure; imagesc(squeeze(W(2,:, :)));
% figure; imagesc(squeeze(W(:, 1, :)));


%% Derivative Filter
display('Calculating first derivatives...')
%first derivatives
dWz = zeros(size(W));
dWz(1:(end-1),:,:) = W(1:(end-1),:,:)-W(2:end,:,:);
dWz(end,:,:)=dWz((end-1),:,:);

dWy = zeros(size(W));
dWy(:,1:(end-1),:) = W(:,1:(end-1),:)-W(:,2:end,:);
dWy(:,end,:)=dWy(:,(end-1),:);

dWx = zeros(size(W));
dWx(:,:,1:(end-1)) = W(:,:,1:(end-1))-W(:,:,2:end);
dWx(:,:,end)=dWx(:,:,(end-1));

%Second Derivatives
display('Calculating second derivatives...')
ddWz = zeros(size(W));
ddWz(1:(end-1),:,:) = dWz(1:(end-1),:,:)-dWz(2:end,:,:);
ddWz(end,:,:)=ddWz((end-1),:,:);

ddWy = zeros(size(W));
ddWy(:,1:(end-1),:) = dWy(:,1:(end-1),:)-dWy(:,2:end,:);
ddWy(:,end,:)=ddWy(:,(end-1),:);

ddWx = zeros(size(W));
ddWx(:,:,1:(end-1)) = dWx(:,:,1:(end-1))-dWx(:,:,2:end);
ddWx(:,:,end)=ddWx(:,:,(end-1));

dW = dWz + dWy + dWx;
ddW = ddWz + ddWy + ddWx;
ddW(ddW>0)=0;

%Not Edited
display ('nearist neighbor filter')
density = 3;

Wnnth=-ddW;
Wnnth = Wnnth./max(Wnnth(:));

% threshold = 0.41;
Wnnth(Wnnth>threshold)=1;
Wnnth(Wnnth<=threshold)=0;
Wnnt=zeros(size(Wnnth));
for z=2:size(Wnnth,1)-1
    for y=2:size(Wnnth,2)-1
        for x=2:size(Wnnth,3)-1
            c=squeeze(sum(sum(sum(Wnnth(z-1:z+1, y-1:y+1, x-1:x+1)))));
            Wnnt(z,y,x)=c>density;
        end
    end
end

ptsr = zeros(size(pts));
Wnntr = Wnnt(:);
Wnntra = size(Wnntr);

figure; hold on
for i=find((Wnntr)>0.5)
    plot3(pts(i,1), pts(i,2), pts(i,3), 'r.')    
end
hold off


for i=find((Wnntr)>0.5)
    ptsr(i,:) = [pts(i,1) pts(i,2) pts(i,3)];
    Wnntra(i) = Wnntr(i);
end

plotbackproject( rightbottomback, lefttopfront, dataset, Wnntr, ptsr, threshold );
figure; imagesc(squeeze(W(:,(size(W,2)/2),:)))
return

%% Ray filter
W = rayFilter(W, threshold);
plotbackproject( rightbottomback, lefttopfront, dataset, W, ptsReshaped, threshold )

%% Deconvolution
Wd = fftn(W)./fftn(W_psf);
Wd = ifftn(Wd);
Wd = ifftshift(Wd);

Wd(Wd<0)=0;
Wd=Wd-min(Wd(:));
Wd=Wd./max(Wd(:));
figure; imagesc(squeeze(abs(Wd(:,:, 15))));
figure; imagesc(squeeze(abs(Wd(20,:, :))));
figure; imagesc(squeeze(abs(Wd(:, 10, :))));
plotbackproject( rightbottomback, lefttopfront, dataset, Wd, ptsReshaped, threshold);
