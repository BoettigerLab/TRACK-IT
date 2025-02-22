
    %% stats
    % Hansen
    Hansen_scope = 'Zeiss Airy scan'; % point-scanning confocal
    Hansen_probe1 = 10e3; %  '224x Tet'; % 10e3; % 224x Tet 20-50 bp spacing and PGK promotor hygromycin resistence = 10kb;
    Hansen_probe2 = 1e3; % + spreading to unknown extent.  
    Hansen_chromatic = NaN; % not reported
    Hansen_zstep = 250; % nm
    Hansen_z_frames = 30;
    Hansen_frame_rate = 20; 
    Hansen_color_mode = 'sequential';
    Hansen_time_per_z = 0.54182; % 2 colors-per-z, 30-z, 250 um steps.  (see supp)
    %   1.3e-6*584*584*30; % =.44/z  13 s motion blur (or 26 s, which doesn't jibe), 
    %    This is per-pixel dwell time, so it must be repeated per z-step 
    Hansen_tot_frames = 365; %  movie we recorded 365 frames [total frames]
    Hansen_loc_error =200; % 186-308 nm, mostly ~200 nm.  ~40-50nm x/y and 117-211 z  (see Table 3 in supp)
    % G_WT = 0.0024 um/s^(1/2);  RMS-dist= 357 nm; 
    
    Gregor_scope = 'Leica SP5'; % point-scanning confocal
    Gregor_frame_rate = 28; % s or 5
    Gregor_frame_rate2 = 5; % s or 5
    Gregor_frame_blur = NaN; % not reported. 8 um stack is acquired during the 28s interval. 
    Gregor_z_frames = 25; % 334 nm steps. 
    Gregor_zstep = 334; % 
    Gregor_loc_error = 180; % pm6 'based on localization control construct'
    Gregor_color_mode = 'sequential'; % 405 (blue-eve+En) and 591 (red-PP7) together, then green (promoter)
    Gregor_chromatic_eror = 46; % norm([21 ,21,35]); % std = 21 dX 35 dZ
    Gregor_chromatic_method = 'tetraspek beads 3D'; % std = 21 dX 35 dZ
    Gregor_probe1 = 'eve-mRNA 11 kb'; % FigS1  expanded endogenous with selection markers
    Gregor_probe2 = 'ParB/ParS'; % 1.3 kb spacer+ 1 kb + spreading
    Gregor_tot_frames = 99;
    Gregor_tot_frames_fast = 143; 
    
    % Mach...Giorgetti
    Luca_scope = 'Nikon Ti-E in HiLo with Delta EMCCD'; % wide field
    Luca_zstep = 300; % nm
    Luca_z_frames = 21; % 6 um stack
    Luca_frame_rate = 30; % s
    Luca_frame_blur = 0.050*21;
    Luca_tot_frames = 400; % 
    Luca_loc_error = 130 ; % pm (average distance between registered tetraspek spots)
    Luca_chromatic_method  = 'tetraspek beads 2D'; % plus manual z of 40 nm
    Luca_color_mode = 'simultaneous 2-cam';
    Luca_probe1 = '120x LacO';
    Luca_probe2 = '140x TetO';
    
    % Alexander
    Alex_scope = 'Nikon + spinning disk'; % CSU-22  EMCCD (Photometric Evolve Delta or Andor iXon)
    Alex_z_frames = 21;  % 21-28
    Alex_z_stp = 300; % 
    Alex_frame_blur = 1.6; % 0.030 * 21
    Alex_chromatic_error =  67; %   norm([norm([12 10 36]), norm([16 16 50])])
    Alex_loc_error = 67; % limited only by chromatic size and SNR matched beads 
    Alex_frame_rate_s = 20 ; % s 
    Alex_tot_frames = 80; % 
    Alex_probe1 = '224x TetO 10kb'; %  10e3; % 
    Alex_probe2 =  '144x cuO 10 kb'; 10e3;%
    Alex_color_mode = 'simultaneous 2-cam';
    
    % Murre
    murre_scope = 'AP OMX 4 EMCCDs'; % (widefield)
    murre_z_frames = 20 ; % 20-30
    murre_z_stp = 500; % nm
    murre_exposure = 0.020 ; % 10 - 20 ms
    murre_chromatic_error = NaN; % not acknowledged (possibly explains the data)
    murre_loc_err = 100 ; % nm % norm([22 22 95])
    murre_loc_err_method = 'mixed GFP/SNAP array';
    murre_probe1 = '360x TetO 15kb'; %  ?360 TetO (original not clearly stated in this work, check the ealier)
    murre_probe2 = '360x TetO-ortho 15kb'; %  360 TetO-orthog
    murre_color_mode = 'simultaneous 2-cam';
    murre_frame_rate = 40;
    murre_frame_rate2 = 2;
    murre_tot_frames = 200;
    murre_tot_frames_fast = 200;
    
    jude_loc_error = 25; % chrom only - step error separate
    jude_tot_frames = 18e3/5;
    jude_frame_rate = 0.5;
    
    %% Gabriele...Hansen Science 2022
    hFolder = 'U:\GenomeData\ByPublication\GabrielleHansen2022\';
    datasets = FindFiles([hFolder,'*.tsv']);
    drop = contains(datasets,'unfiltered');
    filtData = datasets(~drop);
    h_sel_e = 2;
    % 27 = no TAD (10 kb) 
    % 36 = WT
    % 65 = no CTCF sites 
    % aux CTCF 4 hr
    % aux RAD21 4hr
    % aux RAD21 2hr
    
    % select a subset of the experiments to display, using the key above
    % filtData = filtData([1,2,3,6,8,9]); 
     % filtData = filtData([2,4,7,10])  % dataset 3 is a bit off, probably a  real effect of sensitivity to the leaky auxin in the RAD21 
    %
    nExp = length(filtData);
    hansen_traj = cell(1,nExp);
    hansen_traj2D= cell(1,nExp);
    hansen_names = {'10 kb','505 kb','no CTCF-sites','CTCF 0h','CTCF 2h','CTCF 4h','RAD21 0h','RAD21 2h','RAD21 4h','WAPL 0h','WAPL 4h'}
    for e=1:nExp 
        disp(['loading ',filtData{e}]);
        HansenWT_table =ReadTableFile(filtData{e}); 
        T = max(HansenWT_table.t);
        N = length(unique(HansenWT_table.id));
        traj = nan(N,T); % N x T
        X1 = nan(N,T);
        X2 = nan(N,T);
        Y1 = nan(N,T);
        Y2 = nan(N,T);
        Z1 = nan(N,T);
        Z2 = nan(N,T);
        for n=1:N
            i = HansenWT_table.id==n-1;
            dist = HansenWT_table.dist(i);
            time = HansenWT_table.t(i)+1;
            traj(n,time) = dist;
            X1(n,time) = HansenWT_table.x(i);
            Y1(n,time) = HansenWT_table.y(i);
            Z1(n,time) = HansenWT_table.z(i);
            X2(n,time) = HansenWT_table.x2(i);
            Y2(n,time) = HansenWT_table.y2(i);
            Z2(n,time) = HansenWT_table.z2(i);
        end
        dis3D = 1e3*sqrt( (X1-X2).^2 + (Y1-Y2).^2 + (Z1-Z2).^2 );
        disXY = 1e3*sqrt( (X1-X2).^2 + (Y1-Y2).^2  );
        disXZ = 1e3*sqrt( (X1-X2).^2 + (Z1-Z2).^2 );
    
        missingData = isnan(dis3D(:));
        % sum(missingData)/length(missingData);
        nObs = sum(~isnan(dis3D),2);
        [~,idx] = sort(nObs,'descend');
        hansen_traj{e} = dis3D(idx,:);
        hansen_traj2D{e} = disXY(idx,:);
    end
    
    
    %% Mach...Giorgetti Nat. Gen. 2022 data
    % two spot traces    % THESE NAMES ARE NOT CLEAR TO ME.  
    % condition                 deg             ncells      ntracks     na       clone name
    % LacO+TetO		        rad21 0min	        214	        238	        NA	        1A2
    % LacO+TetO	            rad21 120min	    277	        300         NA	
    % LacO+TetO +CTCF	    rad21  0min 	    152	        164	        NA	        1B1
    % LacO+TetO +CTCF	    rad21 120min	    248     	269	        NA	
    % LacO+TetO -prom        rad21 0min         155	        166	        NA  	    1F4
    % LacO+TetO -prom   	rad21 120min	    170	        183 	    NA	
    lucaFolder = 'U:\GenomeData\ByPublication\MachGerogetti2022\Mach_et_al_NG_2022\two_colors\';
    FindFiles([lucaFolder,'*1A2*_0min*']) ; % 
    
    lucaExps{1} = FindFiles([lucaFolder,'*1A2_0min*']) ; % 
    lucaExps{2} = FindFiles([lucaFolder,'*1A2_120min*']);  % 
    lucaExps{3} = FindFiles([lucaFolder,'*1B1_0min*']) ; % 
    lucaExps{4} = FindFiles([lucaFolder,'*1B1_120min*']) ; % with CTCF, no  Cohesin
    lucaExps{5} = FindFiles([lucaFolder,'*1F4_0min*']) ; % 
    lucaExps{6} = FindFiles([lucaFolder,'*1F4_120min*']) ; % with CTCF, no  Cohesin
    
    nExp = length(lucaExps);
    luca_names = {'WT','no RAD21','+CTCF-sites','+CTCF-sites, no RAD21','no prom','no prom, no RAD21'};
    luca_traj = cell(nExp,1);
    luca_traj2D = cell(nExp,1);
    for e=1:nExp % e=1
        lucaCSVs = lucaExps{e};
          T =400;
            N = 150;
            traj = nan(N,T); % N x T
            traj2 = nan(N,T); % N x T
            c=0;
        for f=1:length(lucaCSVs) % by FOV
            lucaTable = readtable(lucaCSVs{f});
            cellIDs = unique(lucaTable.cell);
            nCells = length(cellIDs);
            for n=1:nCells
                i = lucaTable.cell==cellIDs(n);
                dist = lucaTable.distance(i);
                dist2 = sqrt( lucaTable.x(i).^2 + lucaTable.y(i).^2 );
                time = lucaTable.frame(i)+1;
                c=c+1;
                traj(c,time) = dist;
                traj2(c,time) = dist2;
            end
        end
        size(traj)
        traj(traj==0) = nan;
        traj2(traj2==0) = nan;
        % sort data
        % missingData = isnan(traj(:));
        % sum(missingData)/length(missingData);
        nObs = sum(~isnan(traj),2);
        [~,idx] = sort(nObs,'descend');
        luca_traj{e} = 1e3*traj(idx,:);
        luca_traj2D{e} = 1e3*traj2(idx,:);
    end
    
    %% Bruckner...Gregor 2023
    gFolder = 'U:\GenomeData\ByPublication\Gregor_Science2023\';
    datasets = FindFiles([gFolder,'*.csv']);
    filtData = datasets;
    nExp = length(filtData);
    
    % for e=1:nExp  % e=5
    % e = 5;  % 150 kb norm homie
    % e = 6;  % 150kb, 5s res instead of 30
    % e = 7;  % 150 kb no-homie (pretty much off, maybe some stray burst)  
    % e = 1;% ~30 kb? 
    %  e = 2;  % 30 kb no-homie (some expression, proximity correlation much weaker (100 nm closer to like 3 nm closer, from 4x more expression to like 1.02x more expression)  
    
    tg_dists = [58,82,88,149,190,595,3327];
    gregor_traj = cell(nExp,1);
    gregor_traj2D = cell(nExp,1);
    gregor_names = {'58','58 no homie','82','88','149','149 5s','149 no homie','190','595','3327'};
    for e=1:nExp % e=5
        disp(['loading ',filtData{e}]);
        gTable =readtable(filtData{e}); 
        cellIDs = unique(gTable.cell_id);
        nCells = length(cellIDs);
        T = max(gTable.time_point)+1;
        traj = nan(nCells,T);
        traj2D = nan(nCells,T);
        expr =  nan(nCells,T);
        N = nCells;
        X1 = nan(N,T);
        X2 = nan(N,T);
        Y1 = nan(N,T);
        Y2 = nan(N,T);
        Z1 = nan(N,T);
        Z2 = nan(N,T);
        for c=1:nCells
            xy = gTable{gTable.cell_id==cellIDs(c),[9,10]};  % 2D
            xyz = gTable{gTable.cell_id==cellIDs(c),[9,10,11]};  % 3D 2pt
            xyz1 = gTable{gTable.cell_id==cellIDs(c),[3,4,5]};  % 3D 1pt
            xyz2 = gTable{gTable.cell_id==cellIDs(c),[6,7,8]};  % 3D 1pt
            dist = sqrt(sum(xyz.^2,2));
            time = gTable{gTable.cell_id==cellIDs(c),2}+1;
            traj(c,time) = dist;
            traj2D(c,time) = sqrt(sum(xy.^2,2));
            % time = gTable.time_point+1;
            X1(c,time) = xyz1(:,1);
            Y1(c,time) = xyz1(:,2);
            Z1(c,time) = xyz1(:,3);
            X2(c,time) = xyz2(:,1);
            Y2(c,time) = xyz2(:,2);
            Z2(c,time) = xyz2(:,3);
            rna = gTable{gTable.cell_id==cellIDs(c),12};
            expr(c,time) = rna;
        end  
        dis3D = sqrt( (X1-X2).^2 + (Y1-Y2).^2 + (Z1-Z2).^2 );
        disXY = sqrt( (X1-X2).^2 + (Y1-Y2).^2  );
        disXZ = sqrt( (X1-X2).^2 + (Z1-Z2).^2 );
        
        % sort data
        % missingData = isnan(dis3D(:));
        % sum(missingData)/length(missingData);
        nObs = sum(~isnan(dis3D),2);
        [~,idx] = sort(nObs,'descend');
        gregor_traj{e} = dis3D(idx,:);   
        gregor_traj2D{e} = disXY(idx,:);   
    end
    
    
    %% Murre Data - pro-B cells (ABL-transformed)
    % short = imaged every 2s for 400s (200)
    % long = imaged every 40s for 4800s 
    % folder = 'U:\GenomeData\ByPublication\MurreDudko2019\Source Data\RawData\ContractedCells\ShortTermImaging\';
    folder = 'U:\GenomeData\ByPublication\MurreDudko2019\Source Data\RawData\NonContractedCells\LongTermImaging\';
    datFiles = FindFiles([folder,'G1*.dat']);
    numCells = length(datFiles);
    murre_names = {'long','fast'};
    Tmax = 200;
    dis_VDJ = nan(numCells,Tmax);
    time_VDJ = nan(numCells,Tmax);
    dis3D = nan(numCells,Tmax);
    dis2D = nan(numCells,Tmax);
    for c=1:numCells % d=1
        try
            g1t = readtable([folder,'\G1_',num2str(c),'.dat']);
            r1t = readtable([folder,'\R1_',num2str(c),'.dat']);
            [timeM,t1,t2] = intersect(g1t{:,4},r1t{:,4});
            g1 = g1t{t1,1:3};
            r1 = r1t{t2,1:3};
            dis = sqrt( (g1(:,1)-r1(:,1)).^2 + (g1(:,2)-r1(:,2)).^2+ (g1(:,3)-r1(:,3)).^2);
            dis2 = sqrt( (g1(:,1)-r1(:,1)).^2 + (g1(:,2)-r1(:,2)).^2  );
            dis_VDJ(c,1:length(dis)) = dis;
            time_VDJ(c,1:length(dis)) = timeM;
            dis3D(c,timeM) = dis;
            dis2D(c,timeM) = dis2;
        catch
        end
    end
    % missingData = isnan(dis3D(:));
    % sum(missingData)/length(missingData);
    nObs = sum(~isnan(dis3D),2);
    [~,idx] = sort(nObs,'descend');
    murre_traj{1} = dis3D(idx,:)*1e3;    
    murre_traj2D{1} = dis2D(idx,:)*1e3;    
    
    % short term
    folder = 'U:\GenomeData\ByPublication\MurreDudko2019\Source Data\RawData\ContractedCells\ShortTermImaging\';
    datFiles = FindFiles([folder,'G1*.dat']);
    numCells = length(datFiles);
    
    dis_VDJ = nan(numCells,Tmax);
    time_VDJ = nan(numCells,Tmax);
    dis3D = nan(numCells,Tmax);
    dis2D = nan(numCells,Tmax);
    for c=1:numCells % d=1
        try
            g1t = readtable([folder,'\G1_',num2str(c),'.dat']);
            r1t = readtable([folder,'\R1_',num2str(c),'.dat']);
            [timeM,t1,t2] = intersect(g1t{:,4},r1t{:,4});
            g1 = g1t{t1,1:3};
            r1 = r1t{t2,1:3};
            dis = sqrt( (g1(:,1)-r1(:,1)).^2 + (g1(:,2)-r1(:,2)).^2+ (g1(:,3)-r1(:,3)).^2);
            dis2 = sqrt( (g1(:,1)-r1(:,1)).^2 + (g1(:,2)-r1(:,2)).^2  );
            dis_VDJ(c,1:length(dis)) = dis;
            time_VDJ(c,1:length(dis)) = timeM;
            dis3D(c,timeM) = dis;
            dis2D(c,timeM) = dis2;
        catch
        end
    end
    % missingData = isnan(dis3D(:));
    % sum(missingData)/length(missingData);
    nObs = sum(~isnan(dis3D),2);
    [~,idx] = sort(nObs,'descend');
    murre_traj{2} = dis3D(idx,:)*1e3;    
    murre_traj2D{2} = dis2D(idx,:)*1e3;    
    
    %% Alexander Sox2 data
    trackTable = readtable('U:\GenomeData\ByPublication\Alexander2019\elife-41769-supp4-v1.csv');
    dataTypes = unique(trackTable{:,1});
    nExp_A = length(dataTypes);
    alex_traj = cell(nExp_A,1);
    alex_traj2D = cell(nExp_A,1);
    alex_sel = 7; % Sox2-SCR -ESC
    alex_names = {'cntrl-cntrl ESC','cntrl-cntrl MES','cntrl-cntrl NPC','SCR-cntrl ESC','SCR-cntrl MES','SCR-cntrl NPC',...
        'Sox2-SCR esc','Sox2-SCR MES','Sox2-SCR NPC','Sox2-SCRdel ESC','SOX2-del-SCR ESC'};
    for e=1:nExp_A   
        isCurr =  contains(trackTable{:,1},dataTypes{e});
        tab = trackTable(isCurr,:);
        [name,idx] = unique(tab.Locus_ID,'stable');  
        tMax = 80;
        % tstep = 20; % s
        nSpots = length(name);
        dis3D = nan(nSpots,tMax);
        dis2D = nan(nSpots,tMax);
        for n = 1:nSpots-1
            spotTable = tab(idx(n):idx(n+1)-1,:);
            t = spotTable.C1_T_Value_frame ;
            dis3D(n,t) = spotTable.Corrected_XYZ_Distance_um*1e3; % convert to nm
            dis2D(n,t) = spotTable.Corrected_XY_Distance_um*1e3; % convert to nm
        end
        % missingData = isnan(dis3D(:));
        % sum(missingData)/length(missingData);
        nObs = sum(~isnan(dis3D),2);
        [~,idx] = sort(nObs,'descend');
        alex_traj{e} = dis3D(idx,:);    
        alex_traj2D{e} = dis2D(idx,:);    
    end
    
    
    
    
    
    %%================================ combine==============================%%
    
    %% COMBINE DATA
    papers = {hansen_traj, luca_traj, alex_traj, murre_traj, gregor_traj};
    expNames = {hansen_names, luca_names,alex_names,murre_names,gregor_names};
    papersRecent = {hansen_traj, luca_traj, gregor_traj};
    papers2D = {hansen_traj2D, luca_traj2D, alex_traj2D, murre_traj2D, gregor_traj2D};
    papersFast = {murre_traj{2},gregor_traj{2}};
    fav_exp = [2,1,7,1,2]; % 3
    fav_name = {'mESC-505kb','mESC-150kb','mESC-120kb','proB-2.7Mb','fly-54kb'};
    fast_exp = [0,0,0,2,6]; % 5s and 2s resolution
    tstep = [20,30,20,40,28];
    
    
   
    
    %% step-size analysis, all experiments
    totExp = 45;
    expLabel = cat(2,hansen_names,' ', luca_names,' ',alex_names,' ',murre_names,' ',gregor_names,' ');
    
    
    traceCnt = zeros(totExp,1);
    obsCnt = zeros(totExp,1);
    medDist = zeros(totExp,1);
    
    
    tgDist_all = [58    58    82    88     149    149  149     190         595        3327];
    
    genDist = [20,505*ones(1,length(hansen_traj)-1),0,...
        150*ones(1,length(luca_traj)),0,...
        166-111,166-13,166*ones(1,9),0,...
        2.7e3,2.7e3,0,...
        tgDist_all,0];
    
    timeRes = [Hansen_frame_rate*ones(1,length(hansen_traj)),0,...
        Luca_frame_rate*ones(1,length(luca_traj)),0,...
        Alex_frame_rate_s*ones(1,length(alex_traj)),0,...
        murre_frame_rate,murre_frame_rate2,0,Gregor_frame_rate*ones(1,5),...
        Gregor_frame_rate2,Gregor_frame_rate*ones(1,4),0];
    