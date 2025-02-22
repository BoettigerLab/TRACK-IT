%% Trace Analysis
% this version runs without supervision as a *script* rather a GUI/function
% The script is organized into blocks to facilite readability and
% troubleshooting. 


%% Block1 path definitions
% these are obviously system specific, you will need to update them to
% reflect the locations of your current data. 

% Where you would like to save the output data
saveFolder = SetFigureSavePath('U:\Manuscripts\Jude Live Imaging\Data\traces_out\','makeDir',false);
npp = 108; % nanometers per pixel

% Where your hdf5 files from spot-fitting exists
JudeData01 = '\\BlabServer1\JudeData01\';
JudeData02 = '\\BlabServer1\JudeData02\';
JudeData03 = '\\BlabServer1\JudeData03\';

% server2 filepaths
NAS02_Vol1 = '\\169.254.32.249\NAS02_Vol1\';
NAS02_Vol2 = '\\169.254.32.249\NAS02_Vol2\';
NAS02_Vol3 = '\\169.254.32.249\NAS02_Vol3\';
NAS02_Vol4 = '\\169.254.32.249\NAS02_Vol4\';

% distance in 5k of the labeled pair in each experiment
dis_kb = [5,5,20,20,55,55,70,70,134,134,260,260,407,407,799,799,2030,2030,12795,12795,73567,73567];
max_seps =  floor((10*dis_kb).^(1/3)); % anticipated separation based on genomic distance (following the 1/3 scaling) 
dataFolders ={ ... 
                [NAS02_Vol3,'Jude\2024-01-03_1kb\']; 
                [NAS02_Vol3,'Jude\2024-01-05_1kb\']; 
                [NAS02_Vol4,'Jude\2023-11-28_16kb\'];  
                [NAS02_Vol4,'Jude\2023-11-26_16kb\'];  
                [NAS02_Vol4,'Jude\2023-11-30_50kb\']; 
                [NAS02_Vol4,'Jude\2023-12-01_50kb\']; 
                [NAS02_Vol4,'Jude\2023-11-22_76kb\'];  
                [NAS02_Vol4,'Jude\2023-11-24_76kb\']; 
                [NAS02_Vol2,'Jude\2023-12-04_130kb\']; 
                [NAS02_Vol2,'Jude\2023-12-06_130kb\'];      
                [NAS02_Vol2,'Jude\2023-12-18_255kb\'];    
                [NAS02_Vol2,'Jude\2023-12-20_255kb\'];  
                [NAS02_Vol2,'Jude\2023-12-26_400kb\'];   
                [NAS02_Vol2,'Jude\2023-12-28_400kb\'];    
                [NAS02_Vol3,'Jude\2023-12-30_800kb\']; 
                [NAS02_Vol3,'Jude\2024-01-01_800kb\']; 
                [NAS02_Vol2,'Jude\2023-12-22_2mb\'];    
                [NAS02_Vol2,'Jude\2023-12-24_2mb\']; 
                [NAS02_Vol3,'Jude\2024-01-07_12mb\'];  
                [NAS02_Vol3,'Jude\2024-01-13_12mb\'];  
                [NAS02_Vol3,'Jude\2024-01-15_73mb\']; 
                [NAS02_Vol3,'Jude\2024-01-17_73mb\']; 
            };
                  

dataFolders = {
 [NAS02_Vol2,'Jude\2024-04-10_1157\']; % intraTAD 100 kb (G9) 
 [NAS02_Vol2,'Jude\2024-04-08_1157\'];  % intraTAD 100 kb (G9)   
 [JudeData02,'2024-03-21_1137\'];  % interTAD 100 kb (E7)
 [JudeData02,'2024-03-23_1137\'];  % interTAD 100 kb (E7)
 [JudeData02,'2024-03-29_1107\'];  %  interTAD 400 kb (A11)
 [JudeData02,'2024-03-31_1107\'];  %  interTAD 400 kb (A11)
 [JudeData03,'2024-04-18_1187\'];  % intraTAD 400 kb (E9)
 [JudeData03,'2024-04-16_1187\']}; % intraTAD 400 kb (E9)
dis_kb = [100,100,100,100,400,400,400,400];
max_seps =  floor((10*dis_kb).^(1/3)); % anticipated separation based on genomic distance (following the 1/3 scaling) 


%% main processing loop
numDataFolders = length(dataFolders);
for dd = 1:numDataFolders  
    dataFolder = dataFolders{dd};
    maxSep = max_seps(dd);
    hfiles = FindFiles([dataFolder,'*pars*.hdf5'],'fullPath',false);
    bin_tag = ['_pars_',hfiles{1}(18:23)]; 
    disp(bin_tag)
    saveOverlayTraces = false;
    showPlots = false;    
    t_tot = tic;
    sampleDax = FindFiles([dataFolder,'sample*C1.dax']);
    for movie_index = 1:length(sampleDax)      
       trajFiles = FindFiles([saveFolder,'dis3D_10Hz_','exp',num2str(dd,'%02d'),'_well',num2str(movie_index,'%02d'),'*.csv']); % 
       if ~isempty(trajFiles)
           continue
       end    
       disp(['analyzing ',sampleDax{movie_index}])
    
    %% Parameters (most of these do not need to be changed)

    % File names: 
    varIn = {'daxFile1',sampleDax{movie_index},...   % this changes for every FOV  
            'bin_tag',bin_tag,...
            'saveFolder', saveFolder,...
            'alignment_file', [dataFolder,'alignmentData.txt'],...
            'chrom_correct_file',[dataFolder,'tform3D_chromatic.mat']};
    
    
    for morePars = 1  % shorthand, allows this to be collapsed for easier reading     
        defaults = cell(0,3);
        % extra figures to monitor
        defaults(end+1,:) = {'Fig_zLink','integer',0}; %#ok<*SAGROW> % 0 for off
        defaults(end+1,:) = {'Fig_zFit','integer',0}; % 0 for off
        defaults(end+1,:) = {'fig_mergeT','integer',0}; % 0 for off
        defaults(end+1,:) = {'saveMovies','boolean',false}; %  save folder
        % figures for xy-link
        defaults(end+1,:) = {'movAveFig','integer',0}; % per trace
        defaults(end+1,:) = {'spotMapFig','integer',0}; % per trace
        defaults(end+1,:) = {'traceFig','integer',0}; %  per trace
    
        % ID pairs / seed points   
        defaults(end+1,:) = {'maxSep','nonnegative',maxSep}; % 20  % in pixels 
        defaults(end+1,:) = {'seedBinResolution','integer',7}; % 
        defaults(end+1,:) = {'minSpotsPerBin','integer',50}; % 
        defaults(end+1,:) = {'minSpotsPerTrace','integer',400}; % 
        defaults(end+1,:) = {'imSize','array',[]}; %    
        defaults(end+1,:) = {'removeDots','positive',[]}; % manually ID some spots to remove by index
    
        % xy-Linking
        defaults(end+1,:) = {'movAveSteps','positive',10}; % Number of windows in which to split the trace in when computing the moving average. Note, multiplied by z-depth    
        defaults(end+1,:) = {'movAveStepMaxFold','positive',3}; % max fold change in step size relative to local average, used by Remove Jumps
        defaults(end+1,:) = {'moveAveMaxPixStep','positive',7}; % max absolute step size (in pixels) for moving avearage trace, used by Remove Jumps
        defaults(end+1,:) = {'movAveLocalStep','positive',6}; % window size to determine local jumps  (both rounds use this filter)
        % parameters for coarse linking of whole field
        defaults(end+1,:) = {'movAveMaxGap','positive',4}; % Downsampling scale to determine spot clusters (spots more than this number of pixels apart will be considered as separate clusters)    
        defaults(end+1,:) = {'minObsPerAve','positive',10}; % minimum number of observations in the moving average window to be counted as a valid obs in the coarse downsampling
        defaults(end+1,:) = {'maxDistToSeed','positive',35}; % maximum distance to seed point to be counted in the coarse linking
        % parameters for refined measurement of moving average
        defaults(end+1,:) = {'windowR','positive',35}; % 'radius' of window around the original seed point in which valid localizations may occur
        defaults(end+1,:) = {'maxDistFromRoughAve','positive',5}; % Distance in pixels from the rough moving average trace in constructing the refined trace (no binning this time, just a restricted xy window)
        defaults(end+1,:) = {'minObsPerAve2','positive',6}; % minimum number of observations in the moving average window to be counted as a valid obs
        defaults(end+1,:) = {'coarseInterpWindow','positive',20}; % window size of points used for interpolation filter (allows recovery in fine trace of time blocks not detected in coarse trace)
        % parameters for trace assembly
        defaults(end+1,:) = {'maxDistFromMovingAve','positive',4}; % Max distance an acceptable spot can be from its local-time-moving average centroid reported as fold-change relative to moving-average stepsize
        defaults(end+1,:) = {'maxStep','positive',1.5}; % distance in pixels, the max distance from the last localization to be added to the trace (permisive)   
        defaults(end+1,:) = {'maxStepFrame1','positive',30}; % distance in pixels, the max distance from the starting point to start counting (starting point is first non-nan in the time average)
        
    
       % merge sister traces xyz/C 
       pars.figSisLabel = 0; % plot demultiplexed sisters  5
    
        % Merging traces over time
        defaults(end+1,:) = {'maxTraceSep','positive',30};  %       maxTraceSep = 30;     
        defaults(end+1,:) = {'maxOverlapFrac','nonnegative',.01}; %  maxOverlapFrac  = 0.01;
        defaults(end+1,:) = {'maxStepVarIncrease','positive',2};  % fold change increase in step-size variation after trace merge maxStepVarIncrease = 2; 
        defaults(end+1,:) = {'minPointsPerTrace','positive',100};  %  minPointsPerTrace = 100;
       
        % save data-tables
        defaults(end+1,:) = {'saveTables','boolean',true}; %  
        defaults(end+1,:) = {'npp','positive',108}; %  nm per pixel xy
        
        % save Movies
        defaults(end+1,:) = {'cropRadius','integer',15}; % Radius of image around located spot to crop for movies.
    
        % filepath
        defaults(end+1,:) = {'saveFolder','string',''}; 
        defaults(end+1,:) = {'daxFile1','string',''}; 
        defaults(end+1,:) = {'alignment_file','string',''}; % 
        defaults(end+1,:) = {'chrom_correct_file','string',''}; % 
        defaults(end+1,:) = {'bin_tag','string','_2d_iters'}; % 
        defaults(end+1,:) = {'framesToLoad','integer',0}; % 0 = load all frames
        defaults(end+1,:) = {'zDepth','integer',0};  % 0 = autodetect   
        defaults(end+1,:) = {'loadOverlayMovie','boolean',false};
        defaults(end+1,:) = {'verbose','boolean',true};
        defaults(end+1,:) = {'veryverbose','boolean',false};

        pars = ParseVariableArguments(varIn,defaults,mfilename);
        pars.saveOverlayTraces = saveOverlayTraces;   
    end
    
    %% Main function
    
    %% stp 1, load data
    for stp1 = 1   % shorthand, allows this to be collapsed for easier reading
            disp('running...')
            % putting the uigetfiles in the Execute is better
            if isempty(pars.daxFile1)
                [daxName,dataFolder] = uigetfile('*.dax','select a dax movie to load');
                daxFile1 = [dataFolder,daxName];
                
            else
                daxFile1 = pars.daxFile1;
                [dataFolder,~] = fileparts(daxFile1);
            end
            if ~strcmp(dataFolder(end),filesep)
                dataFolder = [dataFolder,filesep]; %#ok<AGROW>
            end
    
            % === Get source folder and calibration files
            % there are some fall-back files for testing. Obviously these
            % hard-coded paths will need updating, though in general you
            % shouldn't be using these anyway, you should provide a new
            % alignment file with every experiement. 
            if isempty(pars.alignment_file)
                folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';  % update me
                pars.alignment_file = [folderCal2,'alignmentData.txt'];
            else
                if pars.verbose
                    disp(['using camera alignment file: ',pars.alignment_file]);
                end
            end
            if isempty(pars.chrom_correct_file)
                folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';
                pars.chrom_correct_file = [folderCal2,'tform3D.mat'];
            else
                if pars.verbose
                    disp(['using camera alignment file: ',pars.chrom_correct_file]);
                end
            end
            alignT = readtable(pars.alignment_file);
            alignS = table2struct(alignT);
            
            % error checking -- make sure files exist
            if ~exist(pars.alignment_file,'file')
                error(['unable to find ',pars.alignment_file, ', check filepath']);
            end
            if ~exist(pars.chrom_correct_file,'file')
                error(['unable to find ',pars.chrom_correct_file, ', check filepath']);
            end
    
            % create a save folder if none was passed
            if isempty(pars.saveFolder)
                pars.saveFolder = [dataFolder,'Analysis\'];
                if ~exist(pars.saveFolder,'dir')
                    mkdir(pars.saveFolder)
                end
            end
    
            % Auto-compute the z-scan cycle from the off-set file
           
            offFile =  regexprep(daxFile1,'_C1.dax','.off');
            offTable = ReadTableFile(offFile,'delimiter',' ');
            stageZ = offTable.offset;
                % auto determine cycle
             if pars.zDepth == 0
                    currZ = 0; oldZ = -inf; 
                    [~,z0] =min(stageZ(1:100));
                    z = z0;
                    nRounds = 100;
                    zDep = nan(nRounds,1);
                    for n=1:nRounds
                        while currZ > oldZ
                            oldZ = stageZ(z);
                            z=z+1;
                            currZ = stageZ(z);
                        end
                        zDep(n)= z-z0;
                        z0 = z; currZ = 0; oldZ = -inf;  % reset
                    end
                    pars.zDepth = median(zDep);
               disp(['system determines the stack height = ',num2str(pars.zDepth)]);
               % figure(10); clf;  % for troubleshooting, remove comment
               % subplot(2,1,1); plot(stageZ(1:100),'.-');
               %  title(['stackHeight = ',num2str(pars.zDepth)'])
               %  subplot(2,1,2); plot(stageZ,'.-');
            end
    
    
            daxFile2 = regexprep(daxFile1,'C1','C2'); % 
            binName1 = regexprep(daxFile1,'.dax',[pars.bin_tag,'.hdf5']); % 
            binName2 = regexprep(daxFile2,'.dax',[pars.bin_tag,'.hdf5']);
           
    
            % === load the images 
            f = floor(height(offTable)/2)+1; % 1
            [im1f,info1] = ReadDax(daxFile1,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);  % note the +3 since its a 4 step series 
            [im2f,info2] = ReadDax(daxFile2,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);
            im1 = max(im1f,[],3);
            im2 = max(im2f,[],3);
            im2 = ApplyReg(fliplr(im2),alignS);
            overlayImage0 = cat(3,im1,im2); % higher contrast 
            overlayImage0 = IncreaseContrast(overlayImage0,'high',.9999,'low',.1);
            figure(1); clf; Ncolor(overlayImage0);  axis image;
    
            % === Load the fit data
            if pars.framesToLoad <= 0 
                nFrames = info1.number_of_frames;  %
            else
                nFrames = pars.framesToLoad;
            end
            pars.nFrames = nFrames; % save for letter
    
            % nFrames = 200; % shorten for testing
            if pars.verbose
                disp('loading channel 1 fits:')
            end    
            fits1 = LoadDaoFits(binName1,'verbose',pars.verbose)
    
            if pars.verbose
                disp('loading channel 2 fits:')
            end
    
            fits2 = LoadDaoFits(binName2,'verbose',pars.verbose) 
    
    
            % === Apply camera alignment and chromatic correction
            fits2c = Register2CamFitTable(fits2,'alignment_file',pars.alignment_file,'chrom_correct_file',pars.chrom_correct_file);
    
            % ==== Load movie for overlays
            %   This loads T frames which the next step will crop to and
            %   assemble as a movie with a slider bar.  This can be a useful
            %   way to browse the data, but is unrelated to the analysis. 
            for p=1:double(pars.loadOverlayMovie) %  (made an -if statement into a loop as matlab script editor has a neat function to hide loops, but not if statements).
                % Load an overlay movie 
                T = 50; % number of frames to show in the movie (slows down processing but useful for validation).   
                selFrames = floor(linspace(1,nFrames-pars.zDepth,T));
                [h,w,~]=size(overlayImage0);
                overlayMovie = zeros(h,w,2,T,'uint16');
                if pars.verbose
                    disp('loading movie frames...')
                end
                for i=1:T
                    if pars.verbose
                        disp([num2str(100*(i)/T),'% complete']);
                    end
                    f=selFrames(i);
                    im1 = ReadDax(daxFile1,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);  % note the +3 since its a 4 step series 
                    im2 = ReadDax(daxFile2,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);
                    im1 = max(im1,[],3);
                    im2 = max(im2,[],3);
                    im2 = ApplyReg(fliplr(im2),alignS);
                    overlayMovie(:,:,1,i) = IncreaseContrast(im1,'high',.99999,'low',.1);
                    overlayMovie(:,:,2,i) = IncreaseContrast(im2,'high',.99999,'low',.1);
                end
                
                if showPlots
                    for f=1:T  % play movie backwards (to end on brightest frame); 
                        figure(1); clf;
                        Ncolor(overlayMovie(:,:,:,T-f+1));  axis image;
                        pause(.01);
                    end
                end
            end
    end 
    
    %% Explore results
    % fig1= MovieSlider(overlayMovie); % optional function to explore
    % results with a simple figure and slider bar to change frame. 
    
    %% step 2: ID pairs
    % 
    fits1(fits1.xsigma>2.5,:) = []; % throw out bad fits
    fits2c(fits2c.xsigma>2.5,:) = [];
    
    for stp2 = 1 
    % (making this a loop of size 1 allows the step to be collapsed using
    % the "-" sign on the left in the matlab text editor, which is
    % sometimes convienent). 
            [h_im,w_im,~]=size(overlayImage0);
            % === auto-select centers for dancing spots 
            [xy1,ptMap1] = TraceCenters({fits1.x},{fits1.y},...
                'minSpotsPerBin',pars.minSpotsPerBin,...
                'minSpotsPerTrace',pars.minSpotsPerTrace,...
                'binResolution',pars.seedBinResolution,...
                 'imSize',[h_im,w_im],...
                'showPlots',1);  % ,'showPlots',true
            pars.minSpotsPerBin2 = pars.minSpotsPerBin;
            [xy2,ptMap2] = TraceCenters({fits2c.x},{fits2c.y},...
                'minSpotsPerBin',pars.minSpotsPerBin2,...
                'minSpotsPerTrace',pars.minSpotsPerTrace,...
                'binResolution',pars.seedBinResolution,...
                'imSize',[h_im,w_im],...
                'showPlots',1);
    
            if showPlots
                figure(4); clf;
                Ncolor(IncreaseContrast(cat(3,ptMap1,ptMap2),'high',.999));
                hold on; plot(xy1(:,1)/pars.seedBinResolution,xy1(:,2)/pars.seedBinResolution,'ro','MarkerSize',20)
                hold on; plot(xy2(:,1)/pars.seedBinResolution,xy2(:,2)/pars.seedBinResolution,'co','MarkerSize',20)
            end
    
             
            % === pair selected spots
            % keep only dancing spots that have partners in the other chromatic channel
            [matched,cost] = MatchPoints(xy1,xy2,'maxDist',pars.maxSep);
            m =matched(cost<pars.maxSep,:);      
            nMatched = size(m,1);
            xy1p = xy1(m(:,1),1:2);
            xy2p = xy2(m(:,2),1:2); % indexed in order of starting point
                
            % ---- remove pairs if requested (based on previous run of step)
            if ~isempty(pars.removeDots)
                xy1p(pars.removeDots,:) = [];
                xy2p(pars.removeDots,:) = [];
                nMatched = size(xy1p,1);
            end
    
            % ---- view
                if exist('overlayMovie','var')
                    tLapseMov = max(overlayMovie,[],4);
                else
                    tLapseMov = overlayImage0;
                end
                f1 = figure(1); clf; Ncolor(tLapseMov);  axis image;
                hold on; plot(xy1(:,1),xy1(:,2),'ro','MarkerSize',20);
                hold on; plot(xy2(:,1),xy2(:,2),'bs','MarkerSize',20);
        
                hold on; plot(xy1p(:,1),xy1p(:,2),'yo','MarkerSize',20);
                hold on; plot(xy2p(:,1),xy2p(:,2),'gs','MarkerSize',20);
        
                
                sNum = cellstr(num2str((1:nMatched)'));
                text(xy1p(:,1)+25,xy1p(:,2),sNum,'color','w'); hold on;
                f1.Position = [0,0,1200,1200]; pause(.1);
        
                SetFigureSavePath(pars.saveFolder,'makeDir',true); 
                [~,movieName] = fileparts(sampleDax{movie_index});
                SaveFigure(f1,'name',['exp',num2str(dd,'%02d'),'_initPairs_',movieName],'formats',{'png'},'overwrite',true);
            
            % --------
    
    end
    
    
    %% Step 3 - LinkXY
    
    for stp3 = 1    
        % the same analysis is applied to both color channels
        % for transparency and troubleshooting, I've written this out
        t_link = tic;
        SetFigureSavePath(pars.saveFolder);
          [spotTrace1,spotTrace2,spotsPerStep] = LinkXYdoublets(fits1,fits2c,xy1p,xy2p,pars.zDepth,...
           'nFrames',nFrames,'parameters',pars);
        t_link = toc(t_link);
        disp(['xy link complete in ',num2str(t_link),' seconds'])
    end
    [nSpots,tObs,dDims] = size(spotTrace1);  %  dDims = 8 data dims:  x,y,z,frame,height,bkd,sigma,significance
    zDepth = pars.zDepth;
    
    
    %% Step 4: Assign sisters in XYZ and match sisters across color channels
     spotTrace1_out = spotTrace1;
     spotTrace2_out = spotTrace2;
      % view the traces
       dims = {'x','y'};
       sisRatios = zeros(nSpots,2);
       for s=1:nSpots  % 
        % detect sisters
            sisRatio_c1 = DetectSisters(squeeze(spotTrace1(s,:,:,:)),'zDepth',pars.zDepth);
            sisRatio_c2 = DetectSisters(squeeze(spotTrace2(s,:,:,:)),'zDepth',pars.zDepth);
            sisRatios(s,:) = [sisRatio_c1,sisRatio_c2];
    
            if showPlots
                figure(1); clf; % s=101 % - clear doublet
                for d=1:2
                    subplot(2,1,d);
                    plot(squeeze(spotTrace1(s,:,d,1)),'.'); hold on;
                    plot(squeeze(spotTrace1(s,:,d,2)),'.')
                    plot(squeeze(spotTrace2(s,:,d,1)),'.'); hold on;
                    plot(squeeze(spotTrace2(s,:,d,2)),'.'); ylabel(dims{d})
                    legend('1a','1b','2a','2b'); title(['doublet trace spot ',num2str(s), ' sis rs=',num2str(sisRatio_c1,3),' ',num2str(sisRatio_c2,3)]);
                    pause(.01);
                end
                if sisRatio_c1 > 0.2 || sisRatio_c2 > 0.2
                    pause(1);
                    [spotTrace1_out(s,:,:,:),spotTrace2_out(s,:,:,:)] = AssignSistersXYZC(...
                        squeeze(spotTrace1(s,:,:,:)),squeeze(spotTrace2(s,:,:,:)),'traceNum',s,...
                        'figSisLabel',pars.figSisLabel) ;
                end
            end
       end    
       g2score = max([sisRatio_c1,sisRatio_c2],2);
       isG2 = g2score > 0.5;
    
    %% Step 5: Fit Z
    % this is a bit slow 
    spotStacks = {spotTrace1_out,spotTrace2_out};
    for stp5 = 1
        t_zfit = tic;
        pntArray = cell(2,2);
        for c=1:2
            pntArray_cs1 = cell(nSpots,1); % a cell array makes this easier / more mem-efficient to parfor loop  
            pntArray_cs2 = cell(nSpots,1); % a cell array makes this easier / more mem-efficient to parfor loop  
            spotStacks_s1 = cell(nSpots,1);
            spotStacks_s2 = cell(nSpots,1);
            for s=1:nSpots % prep the parfor loop
                spotStacks_s1{s} =  squeeze(spotStacks{c}(s,:,:,1)); % only take sister 1   ===================
                spotStacks_s2{s} =  squeeze(spotStacks{c}(s,:,:,2)); % only take sister 2   ===================
            end
            parfor s=1:nSpots %  parfor this  % s=33
                % sort array into z-stack
                currStack1 = nan(nFrames/zDepth,8,zDepth);  %  a 3D matrix, t_Obs x d_Dims x z_Depth.    %  dDims = 8 data dims:  x,y,z,frame,height,bkd,sigma,significance   
                currStack2 = nan(nFrames/zDepth,8,zDepth);  %  a 3D matrix, t_Obs x d_Dims x z_Depth. 
                for z =1:zDepth
                    currStack1(:,:,z) =spotStacks_s1{s}(z:zDepth:end,:) ; 
                    currStack2(:,:,z) =spotStacks_s2{s}(z:zDepth:end,:) ; 
                end
                % fitStack = currStack;
                fitStack1 = FitZ(currStack1,'spotNum',s,'showFig',0); % pars.Fig_zFit
                
                fitStack2 = FitZ(currStack2,'spotNum',s,'showFig',0); % pars.Fig_zFit
                % sort z-stack back into an array
                pntArray_cs1{s} = nan(nFrames,8);
                pntArray_cs2{s} = nan(nFrames,8);
                for z=1:zDepth
                    pntArray_cs1{s}(z:zDepth:end,:) = fitStack1(:,:,z);
                    pntArray_cs2{s}(z:zDepth:end,:) = fitStack2(:,:,z);
                end
            end
            % reorganize into array
            pntArray{1,c} = permute(cat(3,pntArray_cs1{:}),[3,1,2]);
            pntArray{2,c} = permute(cat(3,pntArray_cs2{:}),[3,1,2]);
        end
        t_zfit = toc(t_zfit);
        disp(['z link complete in ',num2str(t_zfit),' seconds']);
    end
    
    
    %% Step 6: ID traces to merge in time
    %  pntArray are the original traces
    %  pntArrayA contains the the merged traces
    % cell array 2-sisters x 2-colors, each containing a 3D matrix: 
    %  c_Cells x t_Obs x d_Dims.
     %  dDims = 8 data dims:  x,y,z,frame,height,bkd,sigma,significance

     maxSisRatios = max(sisRatios,[],2) ;
    for stp6 = 1
        % figure(3); clf; Ncolor(tLapseMov);  axis image; %  figures just for troubleshooting 
        % sNum = cellstr(num2str( (1:size(xy1p,1))' ));
        % text(xy1p(:,1),xy1p(:,2)+5,sNum,'color','y','FontSize',12); hold on;
      
        sis =1 ;
        [nCells,~,~] = size(pntArray{sis,1});
        % --- Merge traces
        % should look at pntArrays for proximity and temporal complimentarity vs.
        % overlap. 
        %   we do this with pdist to allow for arbitary goup size (unlike knn)
        %   and to avoid potential binning splits of adjacent groups. 
        % find seed points that are nearby
        dMap = squareform(pdist(xy1p));
        dMap = triu(dMap);
        dMap(dMap==0) = nan; % ignore self
        dMap(dMap>pars.maxTraceSep) = nan;
        % figure(7); clf; imagesc(dMap); colorbar;  clim([0,10]);
        dMap(dMap<=pars.maxTraceSep) = 1;
        dMap = ~isnan(dMap);
        % figure(7); clf; imagesc(dMap); colorbar; % clim([0,1]);
        hasData = find(sum(dMap,2)>0);
        nGrps = length(hasData);
        grpList = cell(nGrps,1);
        for n=1:nGrps
            currSpot = hasData(n);
            sameGroup = find(dMap(hasData(n),:));
            grpList{n} = [currSpot,sameGroup];
            % remove from map to avoid duplicates
            dMap(sameGroup,:) = false;
            dMap(:,sameGroup) = false;
        end
        grpSize = cellfun(@length,grpList);
        grpList(grpSize<=1) = [];
        
    
        % Merge traces
        toRemove = false(nCells,2);  % updated from nSpots
        pntArrayA = pntArray;
        for c=1:2 % loop over colors
            if pars.fig_mergeT
                figure(7+c); clf;
            end
            % figure(8); clf;
            nD = length(grpList);
            for d=1:nD  % d=19
                if pars.fig_mergeT
                    figure(7+c); clf; subplot(2,1,1); % subplot(nD,1,d);
                    for g=1:length(grpList{d})
                        plot( squeeze(pntArray{sis,c}(grpList{d}(g),:,1)),'.-' ); hold on; 
                        ylabel(grpList{d}(:))
                    end
                    legend();
                end
                % look for overlap
                traceOlap = pntArray{sis,c}(grpList{d},:,1) > 0;
                % figure(8); subplot(nD,1,d); imagesc(traceOlap);
                % --- merge traces 
                % criteria: 
                %      to merge, the enteries should not overlap more than 1%
                %      both traces should have more than 1% of total
                %      Also, post merge, the magnitude of the step-size should not
                %      increase substantially
                obsPerTrace = sum(traceOlap,2);
                totPoints = sum(obsPerTrace);
                uniqueFrames = sum(traceOlap,1)==1;
                totUnique = sum(uniqueFrames);
                if totUnique/totPoints >(1-pars.maxOverlapFrac) &&  min(obsPerTrace/totPoints) > pars.maxOverlapFrac
                    for g=2:length(grpList{d})
                        mergeArray = cat(4,pntArrayA{sis,c}(grpList{d}(1),uniqueFrames,:),  pntArrayA{sis,c}(grpList{d}(g),uniqueFrames,:));
                        newEntry = nansum(mergeArray,4);  
                        % look at step size variation to decide if we keep this merge.    
                        stp1 = nanmean(abs(diff(pntArrayA{sis,c}(grpList{d}(1),uniqueFrames,1))));
                        stpG = nanmean(abs(diff(pntArrayA{sis,c}(grpList{d}(g),uniqueFrames,1))));
                        stpNew = nanmean(abs(diff(newEntry(1,:,1))));  % stepsize mean of merged distribution 
                        stpOrig = nanmean([stp1,stpG]);    % stepsize mean of original distributions 
                        if stpNew/stpOrig < pars.maxStepVarIncrease
                            pntArrayA{sis,c}(grpList{d}(1),uniqueFrames,:) = newEntry;
                            toRemove(grpList{d}(g),c) = true;
                            if pars.fig_mergeT
                                figure(7+c);  subplot(2,1,2); 
                                plot( squeeze(pntArrayA{sis,c}(grpList{d}(1),:,1)),'.-' ); hold on;   % +2
                                title(['spot channel ',num2str(c)]); pause(1);
                            end
                        end
                    end
                end
            end
        end
        disp(['number of merged traces found for chn1 & chn2 = ',num2str(sum(toRemove))])
                
        
       %  %--- show merged traces
       %  f2= figure(2); clf; Ncolor(tLapseMov);  axis image;
       %  sNum = cellstr(num2str( (1:size(xy1p,1))' ));
       %  keep = ~toRemove(:,1);
       %  text(xy1p(keep,1),xy1p(keep,2),sNum(keep),'color','w','FontSize',15); hold on;
       %  text(xy1p(~keep,1),xy1p(~keep,2),sNum(~keep),'color','r','FontSize',15); hold on;
       %  SaveFigure(f2,'name',['removedPairs_',movieName],'formats',{'png'},'overwrite',true);
       % 
       % % for trouble-shooting
       % % add overlay of colored spots; 
       %  xx = pntArray{sis,1}(:,:,1); xx=xx(:);
       %  yy = pntArray{sis,1}(:,:,2); yy=yy(:);
       %  tt = (1:length(xx))';
       %  scatter(xx,yy,[],tt,'o','SizeData',6); colormap('jet'); colorbar;
       % 
       %  xx = pntArray{sis,2}(:,:,1); xx=xx(:);
       %  yy = pntArray{sis,2}(:,:,2); yy=yy(:);
       %  scatter(xx+15,yy+15,[],tt,'s','SizeData',6);
        
    
    
        % also remove traces with too few points:       
        tooFew = false(nCells,2);
        for c=1:2
            tooFew(:,c) = sum(~isnan(pntArrayA{sis,c}(:,:,1)),2)<pars.minPointsPerTrace; 
        end
        disp(['removing from chn1 & chn2 ' , num2str(sum(tooFew)), ' traces with fewer than ',num2str(pars.minPointsPerTrace),' observations']);
        toRemove = toRemove | tooFew;
       
    
        % actually do the removal
        for c=1:2
            pntArrayA{sis,c}(toRemove(:,c),:,:) = [];
            
        end
        
        maxSisRatios(any(toRemove,2)) = [];
        % now we need to re-pair points
        %    some of the traces may have only merged in one channel
        %    some of the traces may have only been too sparse in one channel
    
        % pntArrayB is re-paired version of pntArrayA
        % === Pair / match traces 
        %  previous step matched based on all spots. This step re-matches the data
        %  using only the linked spots. It is also necessary because the Linking
        %  alogrithm is not guarenteed to preserve order [though if it doesn't drop
        %  unlinked starting points, I think I could force it to prserve order].
        xy1 = squeeze(nanmean(pntArrayA{sis,1}(:,:,1:2),2)); % average over frames
        xy2 = squeeze(nanmean(pntArrayA{sis,2}(:,:,1:2),2)); % average over frames
        [matched,cost] = MatchPoints(xy1,xy2,'maxDist',pars.maxSep);
        m = matched(cost<pars.maxSep,:);      
        nMatched = size(m,1);
        pntArrayB{1,1} = pntArrayA{sis,1}(m(:,1),:,:);
        pntArrayB{1,2} = pntArrayA{sis,2}(m(:,2),:,:);
        xy1b = xy1(m(:,1),:);
        xy2b = xy2(m(:,2),:);  
        nPts = size(xy2b,1);
        sNum = cellstr(num2str( (1:nPts)' ));
    end
    
%%

    %% Step 7: compute distances
    

    nCells = size(pntArrayB{sis,1},1);
    dis3D = nan(nCells,tObs); % 10 Hz data summary
    dis2D = nan(nCells,tObs);% 10 Hz data summary
    
    T = tObs;
    disM2 = cell(nCells,1); % 2 Hz data summary
    disM3 = cell(nCells,1);% 2 Hz data summary
    xyzT = cell(nCells,1); % 2 Hz x1 y1 z1 x2 y2 z2
    w = movie_index; % well number  
    for c=1:nCells  % should parfor this   c = 215
        % save absolute x,y position so we can return to the original movies.
        x1 =  pntArrayB{sis,1}(c,:,1)'*npp;  % xc = nanmean(x1); x1 = x1-xc;
        y1 =  pntArrayB{sis,1}(c,:,2)'*npp; % yc = nanmean(y1);  y1 = y1 - yc;
        z1 =  pntArrayB{sis,1}(c,:,3)';    zc = nanmean(z1);  z1 = z1 - zc;
        h1 =  pntArrayB{sis,1}(c,:,5)';
        bkd1 = pntArrayB{sis,1}(c,:,6)';
        t1 = (1:length(x1))';
        x2 =  pntArrayB{sis,2}(c,:,1)'*npp;  % x2 = x2 - xc;
        y2 =  pntArrayB{sis,2}(c,:,2)'*npp;  % y2 = y2 - yc;
        z2 =  pntArrayB{sis,2}(c,:,3)';    z2 = z2 - nanmean(z2);
        h2 =  pntArrayB{sis,2}(c,:,5)';
        bkd2 = pntArrayB{sis,2}(c,:,6)';
        t2 = (1:length(x2))';
        
        
         % Remove noise based on typical jump
        % (this step is slow)
        xy1 = RemoveJumps([x1,y1],'localRegion',10,'maxDiff',3,'maxAbsStep',200,'removeLoners',true);
        xy2 = RemoveJumps([x2,y2],'localRegion',10,'maxDiff',3,'maxAbsStep',200,'removeLoners',true);
        
         % compute distance 2D
        xy1s = fillmissing(xy1,'linear','EndValues','nearest','maxGap',10);
        xy2s = fillmissing(xy2,'linear','EndValues','nearest','maxGap',10);
        dxy = sqrt(sum( (xy1s - xy2s).^2,2));
        nd1 = isnan(xy1(:,1));  % take the original nans
        nd2 = isnan(xy2(:,1));
        if sum(nd1) >= sum(nd2)  % reversed
            nd = nd1;
        else
            nd = nd2;
        end
        dxy(nd,:) = nan; % add them into the saved distances
        xy1s(nd,:) = nan;
        xy2s(nd,:) = nan;
        hasData = ~isnan(z1) & ~isnan(z2);
        z1 = z1 - nanmean(z1(hasData));
        z2 = z2 - nanmean(z2(hasData));
        z1i=z1; 
        z2i=z2;
        z1(isnan(z1i))=z2(isnan(z1i));
        z2(isnan(z2i))=z1(isnan(z2i));
        xyz1 = [xy1s(:,1),xy1s(:,2),z1];
        xyz2 = [xy2s(:,1),xy2s(:,2),z2];
        % compute distance 3D
        %   note, we don't actually add any new data points, but if the
        %   two spots are seen in different z-frames we use the average
        %   z of where the second was last seen to call the distance.
        xyz1s = fillmissing(xyz1,'linear','EndValues','nearest' ,'maxGap',200); % ,'maxGap',10
        xyz2s = fillmissing(xyz2,'linear','EndValues','nearest' ,'maxGap',200); % ,'maxGap',10
        dxyz = sqrt(sum( (xyz1s - xyz2s).^2,2));
        nd1 = isnan(xyz1(:,1)); % take the original nans
        nd2 = isnan(xyz2(:,1));
        if sum(nd1) >= sum(nd2)  % reversed
            nd = nd1;
        else
            nd = nd2;
        end
        dxyz(nd,:) = nan; % add them into the saved distances  
        
        xyz1o = xyz1s;
        xyz1o(nd,:) = nan;
        xyz2o = xyz2s;
        xyz2o(nd,:) = nan;
        time = (1:T)'*.1; % convert to seconds         
        saveTable = table(time);
        saveTable.x1 = xyz1o(:,1);
        saveTable.y1 = xyz1o(:,2);
        saveTable.z1 = z1;
        saveTable.x2 = xyz2o(:,1);
        saveTable.y2 = xyz2o(:,2);
        saveTable.z2 = z2;
        saveTable.dxyz = dxyz;
        
        % compute 2 hz position (once per z-stack)
        dxyz2 = nan(1,T/pars.zDepth);
        ts = 1:5:T;
        ds = zeros(T/pars.zDepth,1);
        for t=1:T/pars.zDepth
            tt = ts(t);
            [dxyz2(t),idx] = min(dxyz(tt:tt+pars.zDepth-1)); % keep only 1 point per z-stack 
            ds(t) = tt-1+idx;
        end
         % correct stage drift with offset file
        z1(isnan(z1) & ~isnan(saveTable.x1)) = 0;
        z2(isnan(z2) & ~isnan(saveTable.x2)) = 0;
        z_off = stageZ*5e3; % convert to nm 
        z1 = z1 + z_off;
        z2 = z2 + z_off;
        saveTable.z1 = z1;
        saveTable.z2 = z2; 
        tempDis3D(d,1:T/pars.zDepth) = dxyz2;
        tempDis3D(d,T/pars.zDepth+1) = sisPerc; % be sure to chop this off
        saveTableOut = saveTable(ds,:);
        skip = isnan(saveTableOut.dxyz);
        saveTableOut(skip,:) = []; % no need to save rows with missing data 
        saveTableOut.cellID = cellID*ones(height(saveTableOut),1);
        saveTableOut.well = w*ones(height(saveTableOut),1);
        saveTableOut.expID = dd*ones(height(saveTableOut),1);
        
        % save table
        e = ceil(dd/2);
        if w <= 4 % first for wells are +dTag
            wo = w;
            name = 'plus-dTag';
        elseif w >4
            wo = w-4;
            name = 'untreated';
        end
        if isG2(c)
            cell_cycle = 'G2';
        else
            cell_cycle = 'G1';
        end        
        saveTableOut.name = repmat(name,height(saveTableOut),1);
        saveTableOut.dis_kb = repmat(dis_kb(dd),height(saveTableOut),1);
        saveTableOut.cell_cycle = repmat(cell_cycle,height(saveTableOut),1);

        saveTableName = ['exp',num2str(dd,'%02d'),'_well',num2str(wo,'%02d'),...
                '_cell',num2str(c),'_',num2str(dis_kb(dd)),'kb_',name,'.csv'];
       writetable(saveTableOut,[saveFolder,saveTableName]);
    end
           

    %%  Step9: save some image stats
    if saveOverlayTraces
        for stp9 = 1  
            % ==== create & save a tiled 2-color image array of the spot-data 
            %   using the the middle frame 
            T = size(overlayMovie,4);
            tLapseMov = overlayMovie(:,:,:,floor(T/2)); % middle frame
            for p=0:double(pars.loadOverlayMovie) %   (made an -if statement into a loop as matlab script editor has a neat function to hide loops, but not if statements).  
                nP = size(xy1c,1);
                xf = round(xy1c(:,1));
                yf = round(xy1c(:,2));
                r= pars.spotDisplayRadius;
                [h,w,~] = size(tLapseMov); % (this one is big, we keep it in 'data' to avoid duplication)   
                imTile1 = zeros(r,r,nP,'uint16');
                imTile2 = zeros(r,r,nP,'uint16');
                for p=1:nP % loop over points
                    x1 = max([1,xf(p)-r]);
                    x2 = min([w,xf(p)+r]);
                    y1 = max([1,yf(p)-r]);
                    y2 = min([h,yf(p)+r]);
                    xL = x2-x1+1;
                    yL = y2-y1+1;
                    if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                        disp(([y1,y2,x1,x2]))
                        disp(([h,w]));
                        disp('debug here')
                    end
                    im1 = squeeze(tLapseMov(y1:y2,x1:x2,1));
                    im2 = squeeze(tLapseMov(y1:y2,x1:x2,2));
                    im1 = IncreaseContrast(im1,'high',.99999,'low',.1);
                    im2 = IncreaseContrast(im2,'high',.99999,'low',.1);
                    imTile1(1:yL,1:xL,p) = im1;
                    imTile2(1:yL,1:xL,p) = im2;
                end
                % the display is the slow
                t=1;
                [~,labelOffsets,imT1] = TileImage(imTile1(:,:,:),'colormap',gray(256),'numRows',10);
                [~,~,imT2] = TileImage(imTile2(:,:,:),'colormap',gray(256),'numRows',10);            
                imO = Ncolor(cat(3,imT1,imT2));
                f1 = figure(1); clf; imagesc(imO);
                nums = cellstr(num2str( (1:nP)' ));
                text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
            end
            f1.Position = [0 0 1200 1200];  pause(.1); % keep a constant size for export;
            SaveFigure(f1,'name',['singlePointRC_',movieName],'formats',{'png'},'overwrite',true);
        
        
            % === create a tiled 2-color superposition of all time as an array
            tLapseMov = max(overlayMovie,[],4);
            for p=0:double(pars.loadOverlayMovie) %   (made an -if statement into a loop as matlab script editor has a neat function to hide loops, but not if statements).  
                nP = size(xy1c,1);
                xf = round(xy1c(:,1));
                yf = round(xy1c(:,2));
                r= pars.spotDisplayRadius;
                [h,w,~] = size(tLapseMov); % (this one is big, we keep it in 'data' to avoid duplication)   
                imTile1 = zeros(r,r,nP,'uint16');
                imTile2 = zeros(r,r,nP,'uint16');
                for p=1:nP
                    x1 = max([1,xf(p)-r]);
                    x2 = min([w,xf(p)+r]);
                    y1 = max([1,yf(p)-r]);
                    y2 = min([h,yf(p)+r]);
                    xL = x2-x1+1;
                    yL = y2-y1+1;
                    if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                        disp(([y1,y2,x1,x2]))
                        disp(([h,w]));
                        disp('debug here')
                    end
                    im1 = squeeze(tLapseMov(y1:y2,x1:x2,1));
                    im2 = squeeze(tLapseMov(y1:y2,x1:x2,2));
                    im1 = IncreaseContrast(im1,'high',.99999,'low',.1);
                    im2 = IncreaseContrast(im2,'high',.99999,'low',.1);
                    imTile1(1:yL,1:xL,p) = im1;
                    imTile2(1:yL,1:xL,p) = im2;
                end
                % the display is the slow
                t=1;
                [~,labelOffsets,imT1] = TileImage(imTile1(:,:,:),'colormap',gray(256),'numRows',10);
                [~,~,imT2] = TileImage(imTile2(:,:,:),'colormap',gray(256),'numRows',10);            
                imO = Ncolor(cat(3,imT1,imT2));
                f1 = figure(1); clf; imagesc(imO);
                nums = cellstr(num2str( (1:nP)' ));
                text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
            end
            f1.Position = [0 0 1200 1200];  pause(.1); % keep a constant size for export;
            SaveFigure(f1,'name',['timeLapseRC_',movieName],'formats',{'png'},'overwrite',true);
            
            % === create a tiled array of the 2D trajectories of each point 
            %   in 2-colors as a simple plot
            [nPts,tObs,dDat] = size(pntArrayB{1});
            gap = 50;
            f2 = figure(2); clf;
            for n=1:nPts
                [rowNum,colNum] = ind2sub([10,ceil(nPts/10)],n);
                x1 = pntArrayB{1}(n,:,1) - nanmean(pntArrayB{1}(n,:,1) ) + labelOffsets(n,1); %#ok<*NANMEAN>
                y1 = pntArrayB{1}(n,:,2) - nanmean(pntArrayB{1}(n,:,2) ) + labelOffsets(n,2);% rowNum*gap;
                x2 = pntArrayB{2}(n,:,1) - nanmean(pntArrayB{2}(n,:,1) ) + labelOffsets(n,1);% + colNum*gap;
                y2 = pntArrayB{2}(n,:,2) - nanmean(pntArrayB{2}(n,:,2) ) + labelOffsets(n,2);% + rowNum*gap;
                plot(x1(~isnan(x1)),y1(~isnan(x1)),'.-','color',[0 .5 0],'linewidth',.1); hold on;
                plot(x2(~isnan(x2)),y2(~isnan(x2)),'.-','color',[.5 0 .5],'linewidth',.1); hold on;
                set(gca,'YDir','reverse')
                % pause(.3);
            end
            f2.Position = [0 0 1200 1400];  pause(.1);  % keep a constant size for export;
            SaveFigure(f2,'name',['traceDataMG_',movieName],'formats',{'png'},'overwrite',true);
        
        
        % === Create 2-pannel tiled array of color-as-time of all the spots
            for p=0:double(pars.loadOverlayMovie)
                nP = size(xy1c,1);
                xf = round(xy1c(:,1));
                yf = round(xy1c(:,2));
                r= pars.spotDisplayRadius;
                [h,w,~] = size(tLapseMov); % (this one is big, we keep it in 'data' to avoid duplication)   
                imTile1 = cell(nP,1); % zeros(r,r,nP,T,'uint16');
                imTile2 = cell(nP,2); % zeros(r,r,nP,T,'uint16');
                for p=1:nP
                    x1 = max([1,xf(p)-r]);
                    x2 = min([w,xf(p)+r]);
                    y1 = max([1,yf(p)-r]);
                    y2 = min([h,yf(p)+r]);
                    xL = x2-x1+1;
                    yL = y2-y1+1;
                    if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                        disp(([y1,y2,x1,x2]))
                        disp(([h,w]));
                        disp('debug here')
                    end
                    im3D = zeros(2*r+1,2*r+1,T,'uint16');
                    im1 = squeeze(overlayMovie(y1:y2,x1:x2,1,:));
                    [hi,wi,~] = size(im1);
                    im3D(1:hi,1:wi,:) = IncreaseContrast(im1,'high',.9999,'low',.1);
                    imTile1{p} = Ncolor(im3D,'colorMax',true,'colormap','cool');
            
                    im3D = zeros(2*r+1,2*r+1,T,'uint16');
                    im2 = squeeze(overlayMovie(y1:y2,x1:x2,2,:));
                    [hi,wi,~] = size(im2);
                    im3D(1:hi,1:wi,:) = IncreaseContrast(im2,'high',.9999,'low',.1);
                    imTile2{p} = Ncolor(im3D,'colorMax',true,'colormap','cool');
            
                end
                % figure(3); clf; TileImageStack(cat(4,imTile1{:}));
                imStk1 = cat(4,imTile1{:}); % x y RGB cell
                imStk2 = cat(4,imTile2{:}); % x y RGB cell
            
                imRGB1 = cell(3,1); 
                imRGB2 = cell(3,1); 
                for c=1:3
                [~,labelOffsets,imRGB1{c}] = TileImage(squeeze(imStk1(:,:,c,:)),'colormap',gray(256),'numRows',10);
                [~,labelOffsets,imRGB2{c}] = TileImage(squeeze(imStk2(:,:,c,:)),'colormap',gray(256),'numRows',10);
                end
                rgb1 = cat(3,imRGB1{:});
                rgb2 = cat(3,imRGB2{:});
                nums = cellstr(num2str( (1:nP)' ));
            end
            f3 = figure(3); clf; 
            f3.Position = [0 0 3000 1000];  pause(.1);
            subplot(1,2,1); imagesc(rgb1); text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
            subplot(1,2,2); imagesc(rgb2); text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
            set(gca,'color','k');  pause(.1);
            SaveFigure(f3,'name',['timeLapseColor_',movieName],'formats',{'png'},'overwrite',true);
        
        
        % === For all points, create a time-overlay trace
        if pars.saveOverlayTraces
            for stp9b=1
                % this is a bit slow
                [nPts,tObs,dDat] = size(pntArrayB{1});
                frameT = (1:tObs);
                f3= figure(3); clf;
                 r = 20;
                xf = round(xy1c(:,1));
                yf = round(xy1c(:,2));
                for n=1:nPts   
                        x1 = max([1,xf(n)-r]);
                        x2 = min([w,xf(n)+r]);
                        y1 = max([1,yf(n)-r]);
                        y2 = min([h,yf(n)+r]);
                        xL = x2-x1+1;
                        yL = y2-y1+1;
                        if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                            disp(([y1,y2,x1,x2]))
                            disp(([h,w]));
                            disp('debug here')
                        end
                        im3D = zeros(2*r+1,2*r+1,T,'uint16');
                        im1 = squeeze(overlayMovie(y1:y2,x1:x2,1,:));
                        [hi,wi,~] = size(im1);
                        im3D(1:hi,1:wi,:) = IncreaseContrast(im1,'high',.9999,'low',.1);
                        rgb1 = Ncolor(im3D,'colorMax',true,'colormap','cool');
                
                        im3D = zeros(2*r+1,2*r+1,T,'uint16');
                        im2 = squeeze(overlayMovie(y1:y2,x1:x2,2,:));
                        [hi,wi,~] = size(im2);
                        im3D(1:hi,1:wi,:) = IncreaseContrast(im2,'high',.9999,'low',.1);
                        rgb2 = Ncolor(im3D,'colorMax',true,'colormap','cool');
                    x1 = pntArrayB{1}(n,:,1) - nanmean(pntArrayB{1}(n,:,1) ) + r; % labelOffsets(n,1); %#ok<*NANMEAN>
                    y1 = pntArrayB{1}(n,:,2) - nanmean(pntArrayB{1}(n,:,2) ) + r; %  labelOffsets(n,2);% rowNum*gap;
                    x2 = pntArrayB{2}(n,:,1) - nanmean(pntArrayB{2}(n,:,1) ) + r; % labelOffsets(n,1);% + colNum*gap;
                    y2 = pntArrayB{2}(n,:,2) - nanmean(pntArrayB{2}(n,:,2) ) + r; %labelOffsets(n,2);% + rowNum*gap;
                    f3 = figure(3); clf;
                    subplot(1,2,1); imagesc(rgb1); hold on;
                    plot(x1(~isnan(x1)),y1(~isnan(x1)),'-','color',[0 .5 0],'linewidth',.1); hold on;
                    plot(x2(~isnan(x2)),y2(~isnan(x2)),'-','color',[.5 0 .5],'linewidth',.1); hold on;
                    scatter(x1(:),y1(:),[],frameT(:),'o','SizeData',6); colormap('jet'); hold on; set(gca,'YDir','reverse')
                    subplot(1,2,2); imagesc(rgb2); hold on;
                    plot(x1(~isnan(x1)),y1(~isnan(x1)),'-','color',[0 .5 0],'linewidth',.1); hold on;
                    plot(x2(~isnan(x2)),y2(~isnan(x2)),'-','color',[.5 0 .5],'linewidth',.1); hold on;
                    scatter(x2(:),y2(:),[],frameT(:),'o','SizeData',6); colormap('jet');  hold on; set(gca,'YDir','reverse')
                   f3.Position = [0 0 2000 800];   pause(.1); % keep a constant size for export;
                   SaveFigure(f3,'name',['traceDataColorCode_',movieName,'_spot',num2str(n)],'formats',{'png'},'overwrite',true,'verbose',false);
                   pause(.1);
                end
            end
        
        end
        end
    
    end
    %% Step 10, save mini-movies (optional)
    for stp10 =1
            if pars.saveMovies
                disp(['loading raw movie data...']);
                t_sav = tic;
                daxName1 = daxFile1;
                daxName2 = regexprep(daxName1,'C1.dax','C2.dax');
                [nCells,nFrames,~] = size(pntArrayB{1});
                xf = nan(nCells,1);
                yf = nan(nCells,1);
                r = pars.cropRadius;
                nT = ceil(nFrames/pars.zDepth);
                im5D = zeros(2*r+1,2*r+1,pars.zDepth,2,nT,nCells,'uint16');
                t=0;
                figure(1); clf; 
                k = 1; % counter of each frame in a batch
                b = 1; % counter for the batch number
                batchSize = 100;
               for f = 1:1:nFrames
                        % read movie in chunks.
                        %   when running in parpool this should reduce the number of 
                        %   simultaneous read/write commands trying to pull from the 
                        %   same disk.  
                        if k==1
                            ff = batchSize*(b-1)+k; 
                            fe = min(ff+batchSize-1,nFrames);
                            % disp(ff)
                            im1_dax = ReadDax(daxName1,'startFrame',ff,'endFrame',fe,'verbose',false); 
                            im2_dax = ReadDax(daxName2,'startFrame',ff,'endFrame',fe,'verbose',false);
                            if pars.verbose
                                disp([num2str(100*ff/nFrames),'% movie loaded'])
                            end
                        else
                            if k==batchSize
                                k=0; % reset
                                b=b+1;
                            end
                        end
                         % disp([f,k,b,ff])
                        im1 = im1_dax(:,:,k+1);
                        im2 = im2_dax(:,:,k+1);
            %             im1 = ReadDax(daxName1,'startFrame',f,'endFrame',f,'verbose',false); 
            %             im2 = ReadDax(daxName2,'startFrame',f,'endFrame',f,'verbose',false);
                        im2 = ApplyReg(fliplr(im2),alignS);
                        im3 = cat(3,im1,im2);
                        [h,w] = size(im1);
                        % im3 = IncreaseContrast(im3,'high',.99995,'low',.001); % constant contrast
                        
                        % figure(1); clf; Ncolor(im3);
                    o = rem(f,pars.zDepth); 
                    if o==0
                        o=pars.zDepth; % keeping indexes straight
                    elseif o==1
                        t=t+1; % starting new series
                    end
                    for s=1:nCells % s =4
                        try
                       % update centering
                            x = round(pntArrayB{1}(s,f,1)); % in pixels :)
                            y = round(pntArrayB{1}(s,f,2));
                            if ~isnan(x)
                                xf(s) = x;
                                yf(s) = y;
                            end
                        % plot
                        if ~isnan(xf(s))
                            x1 = max([1,xf(s)-r]);
                            x2 = min([w,xf(s)+r]);
                            y1 = max([1,yf(s)-r]);
                            y2 = min([h,yf(s)+r]);
                            imS = im3(y1:y2,x1:x2,:);
                            xL = x2-x1+1;
                            yL = y2-y1+1;
                            if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                                disp(([y1,y2,x1,x2]))
                                disp(([h,w]))
                                disp('debug here')
                            end
                            im5D(1:yL,1:xL,o,:,t,s) = imS; 
                        end
                        catch er
                            warning(['error processing spot ',num2str(s), ' frame ',num2str(f)])
                            disp(([y1,y2,x1,x2]))
                            warning(er.message);
                            warning(er.getReport);
                            disp('debug here')
                        end
                    end
                    % set(gcf,'color','k');
                    % pause(.01);
                    k=k+1;
                end
                
            
                %% save files
                [~,daxName] = fileparts(daxFile1);
                daxName = regexprep(daxName,'_C1.dax','');
                saveMovieFolder = [pars.saveFolder,daxName,'_movies\'];
                disp(['saving cropped movies as i5d...']);
                if ~exist(saveMovieFolder,'dir')
                    mkdir(saveMovieFolder);
                end          
                for s=1:nCells
                    im5Dout = im5D(:,:,:,:,:,s);
                    % [dH,dW,dZ,dC,dT] = size(im5Dout);
                     WriteImage5D(im5Dout,[saveMovieFolder,'spot_',num2str(s,'%03d'),'.i5d'],'lociX',xf(s),'lociY',yf(s));
                end    
                t_sav = toc(t_sav)
                disp(['spent ',num2str(t_sav/60,3),' minutes saving movies']);
            end
    
    end
    
    
    end
    t_tot = toc(t_tot);
    disp(['completed in ',num2str(t_tot/60,3),' minutes']);
end
