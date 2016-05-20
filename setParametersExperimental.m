platform = 'Jess';

% Set Parameters for loading Experimental Data
    datafile = '20140829-Tcloseandtwopatches_nocalscreen';
    updateParameters = true; %choose to load or update parameters
    pf = 54.93e6; % pulse frequency hz
    
    target_est = [0, 44, 55]; %remove later since this is the same as parameters.target
    rightbottomback = [-150, -10, 10];
    lefttopfront = [150, 80, 180];
    gridsize = 2; %cm    
        
% If parameters file already exists in the data folder, it will load that instead 
    parameters.fbduration=8000; % Length of first bounce in pixels
    parameters.duration=14000; % Length of the pulse, pixels
%     parameters.offset = 51; %cm from the wall 
    parameters.offset = 0; %cm from the wall
    
    parameters.target=[0 44 55]; % cm
    parameters.period=ceil(1/pf*1e12); % when next pulse should begin, pixels (time bins)
    parameters.tfudgecm=0;
%     parameters.deltalaser=[-9 3 -7]; % Correction Factor , cm
    parameters.deltalaser=[0 0 0]; % Correction Factor , cm
    parameters.deltacameraPos=[0 0 0]; % Correction Factor, cm
    
%define equipment positions and the number of time bins --- distances in cm
    dataset.laserOrigin=[-51.0 13.2 149.9] + parameters.deltalaser; 
    dataset.cameraOrigin=[-23.3 18.4 117.4];
    dataset.cameraPos=[-30.1 36.1 0] + parameters.deltacameraPos; 
    dataset.deltat = 1e-12*c*100; % size of the time bin, cm
    dataset.t = parameters.duration; %number or time bins
    

%% Origanize data  
% Datapaths for different Platforms
    switch platform
        
        case 'Jess'
            data_path = ['c:\Users\Jess\My Documents\LOCI\MIRIS\Results\' datafile '\']; 
            pic_path = ['c:\Users\Jess\My Documents\LOCI\MIRIS\Results\' datafile '\Pics\'];
            cam_Cal = 'C:\Users\Jess\My Documents\LOCI\MIRIS\camCal\calibration.bmp';
    end
    
    if exist([data_path 'parameters.mat'], 'file') && updateParameters == true
        b = input('Overwritting saved parameters.  Press y to continue overwritting or n to load current parameters\n', 's');
        if b == 'n'
            updateParameters = false;
        end
    end
            
    switch updateParameters
        
        case true           
            save([data_path 'parameters.mat'], 'parameters');
            display('Updating parameters...');
        
        case false
            load ([data_path 'parameters.mat']);
            display('Loading parameters...');
    end
    