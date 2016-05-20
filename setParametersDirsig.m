
platform = 'Andreas';

% Set Parameters for loading data from Dirsig simulation
    datafile = '20150311MoundWSphere-shifted';
    dataLimit = [2500 4000];
    
    target_est = [50, 0, 50];
    rightbottomback = [-100, -100, 0];
    lefttopfront = [100, 100, 100];
    gridsize = 6; %m

   
%% Origanize data  
% Datapaths for different Platforms
    switch platform
        
        case 'Jess'
            data_path = ['c:\Users\Jess\Dropbox\dirsig_results\' datafile '\']; %the path to the .bn and .img files with the actual data   
            sl=load('c:\Users\Jess\Dropbox\dirsig_scripted\spotlist.mat'); %laser spots from makepatternsLoop
            dirsig_path='c:\Program Files\DIRSIG 4\'; %the path to the dirsig installation   

        case 'Andreas'
            data_path = ['/home/andreas/Documents/Dropbox/research/projects/periscope/dirsig_results/' datafile '/']; %the path to the .bn and .img files with the actual data   
            sl=load('/home/andreas/Documents/Dropbox/dirsig_scripted/spotlist.mat'); %laser spots from makepatternsLoop
            dirsig_path='/home/andreas/Applications/dirsig/'; %the path to the dirsig installation   

    end
       
% setUp metadata class describing the dataset
    metadata.data_basename='lidar_cave_jess';
    metadata.truth_basename='lidar_pathlength';
    metadata.dataLimits=dataLimit;
    metadata.spotlist = sl.spotlist;