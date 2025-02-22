%% load data
% dataset 1
tic
dataFolder = 'U:\Manuscripts\Jude Live Imaging\Data\traces_v13\';
dis_kb = [5,20,55,70,134,260,407,799,2030,12795,73567];
nE = length(dis_kb);
data3D = cell(2,11);
for e=1:nE    
    datafile = [dataFolder,'Dis3D_2Hz_noTreatment_',num2str(dis_kb(e)),'kb.csv'];
    data3D{1,e} = readmatrix(datafile);
    datafile = [dataFolder,'Dis3D_2Hz_dTag_',num2str(dis_kb(e)),'kb.csv'];
    data3D{2,e} = readmatrix(datafile);
end

% load exp 2
dis_kb2 = [100,100,400,400];
tadPos ={'intra','inter','inter','intra'};
data3D_2 = cell(2,4); 
for e=1:4    
    datafile = [dataFolder,'Dis3D_2Hz_Reg2_noTreatment_',tadPos{e},'_',num2str(dis_kb2(e)),'kb.csv'];
    data3D_2{1,e} = readmatrix(datafile);
    datafile = [dataFolder,'Dis3D_2Hz_Reg2_dTag_',tadPos{e},'_',num2str(dis_kb2(e)),'kb.csv'];
    data3D_2{2,e} = readmatrix(datafile);
end
toc

%% load published data for comparison data
% C:\Users\Alistair.BOETTIGERSERVER\Desktop\code\boettiger-notebook\LiveORCA\Figs_v13\
CompAllPub_SuppFig1

%% all cells, next obs  - Inter-quartile range
exp_F2D = zeros(2,11);
exp_FnD = zeros(2,11,2);
exp_nCells = zeros(2,11);
exp_nObs = zeros(2,11);
T = 3600;
st = 1;
for d=1:2
    for e=1:11
        dis3D = data3D{d,e}(:,1:T);
        nObs = sum(~isnan(dis3D),2);
        dis3D(nObs<30,:) = [];
        nC = size(dis3D,1);
        exp_nCells(d,e) = nC;
        exp_nObs(d,e) = sum(nObs);
        stp1_h = nan(nC,1);
        iqr_h = nan(nC,1);
        for c=1:size(dis3D,1)
            trac = dis3D(c,:);
            trac(isnan(trac)) = [];
            stp1_h(c) = nanmedian(abs(trac(:,1+st:st:end)-trac(:,st:st:end-st)),2);
            mx = .75; % (length(trac)-4)/length(trac);
            mn =  .25;% 4/length(trac);
            iqr_h(c) = quantile(trac,mx,2) - min(trac);
        end
        exp_F2D(d,e) = nanmedian(stp1_h./iqr_h);
        exp_FnD(d,e,1) = nanmedian(stp1_h);
        exp_FnD(d,e,2) = nanmedian(iqr_h);
    end
end
stp1 =  exp_FnD(:,:,1);
nanmedian(stp1(:))
nanmedian(exp_F2D(:))
 

exp2_F2D = zeros(2,4);
exp2_FnD = zeros(2,4,2);
exp2_nCells = zeros(2,4);
exp2_nObs = zeros(2,4);
T = 3600;
st = 1;
for d=1:2
    for e=1:4
        dis3D = data3D_2{d,e}(:,1:T);
        nObs = sum(~isnan(dis3D),2);
        dis3D(nObs<30,:) = [];
        nC = size(dis3D,1);
        exp2_nCells(d,e) = nC;
        exp2_nObs(d,e) = sum(nObs);
        stp1_h = nan(nC,1);
        iqr_h = nan(nC,1);
        for c=1:size(dis3D,1)
            trac = dis3D(c,:);
            trac(isnan(trac)) = [];
            stp1_h(c) = nanmedian(abs(trac(:,1+st:st:end)-trac(:,st:st:end-st)),2);
            mx = .75;  % interquartile range - 75th percentile
            mn =  .25; %  25th percentile
            iqr_h(c) = quantile(trac,mx,2) - min(trac);  % IQR
        end
        exp2_F2D(d,e) = nanmedian(stp1_h./iqr_h);
        exp2_FnD(d,e,1) = nanmedian(stp1_h);
        exp2_FnD(d,e,2) = nanmedian(iqr_h);
    end
end


pub_F2D = nan(5,12);
pub_FnD = nan(5,12);
pub_nCells = nan(5,13);
pub_nObs  = nan(5,13);
for p=1:5
    nE = length(papers{p});
    for e=1:nE
        dis3D = papers{p}{e};
        dis3D(dis3D==0) = nan;
        nObs = sum(~isnan(dis3D),2);
        dis3D(nObs<10,:) = [];
        nC = size(dis3D,1);
        pub_nCells(p,e) = nC;
        pub_nObs(p,e) = sum(nObs);
        stp1_h = nan(nC,1);
        iqr_h = nan(nC,1);
        for c=1:size(dis3D,1)
            trac = dis3D(c,:); trac(trac==0) = nan;
            trac(isnan(trac)) = [];
            stp1_h(c) = nanmedian(abs(trac(:,1+st:st:end)-trac(:,st:st:end-st)),2);
            mx = .75; % interquartile range
            mn = .25; %  
            iqr_h(c) = quantile(trac,mx,2) - quantile(trac,mn,2);
        end
        pub_F2D(p,e) = nanmedian(stp1_h./iqr_h); % nanmedian  nanmean
        pub_FnD(p,e,1) = nanmedian(stp1_h);
        pub_FnD(p,e,2) = nanmedian(iqr_h);
    end
end

pub_FnD(pub_FnD==0) = nan;
pub_F2D(pub_F2D==0) = nan;

pub_stp = pub_FnD(:,:,1)';
nanmean(pub_stp,1)

figure(100); clf;
exp_stp = exp_FnD(:,:,1)';
exp_rng = exp_FnD(:,:,2)';
subplot(1,2,1); bar([nanmean(pub_FnD([1:5],:,2),2);nanmean(exp_rng(:))]); hold on; 
ylabel('3D dist. (nm)'); set(gca,'XTick',1:6,'XTickLabel',paperNames([1:5,8])); 
subplot(1,2,1); bar([nanmean(pub_FnD([1:5],:,1),2);nanmean(exp_stp(:))]); box off;
subplot(1,2,2); bar([nanmean(pub_F2D([1:5],:),2);nanmean(exp_F2D(:))]);
box off;
set(gca,'XTick',1:6,'XTickLabel',paperNames([1:5,8]));
ylabel('F to D ratio'); 

figure(101); clf;
exp_stp = exp_FnD(:,:,1)';
subplot(1,2,1); bar([nanmean(pub_FnD([1:5],:,1),2);nanmean(exp_stp(:))]); box off;
subplot(1,2,2); bar([nanmean(pub_F2D([1:5],:),2);nanmean(exp_F2D(:))]);  box off;

figure(102); clf; 
pub_stp = pub_FnD(:,:,2)';
exp_stp = exp_FnD(:,:,2)';
bar([pub_stp(pub_stp~=0 & ~isnan(pub_stp)); exp_stp(:)  ]); hold on;
pub_stp = pub_FnD(:,:,1)';
exp_stp = exp_FnD(:,:,1)';
bar([pub_stp(pub_stp~=0 & ~isnan(pub_stp)); exp_stp(:)  ]);
box off;
ylabel('3D dist. (nm)');



all_GenDist = [20,505*ones(1,length(hansen_traj)-1),...
    150*ones(1,length(luca_traj)),...
    166-111,166-13,166*ones(1,9),...
    2.7e3,2.7e3,...
    tgDist_all,...
    dis_kb,dis_kb,...
    dis_kb2,dis_kb2];
nE = length(all_GenDist);


timeRes = [Hansen_frame_rate*ones(1,length(hansen_traj)),...
    Luca_frame_rate*ones(1,length(luca_traj)),...
    Alex_frame_rate_s*ones(1,length(alex_traj)),...
    murre_frame_rate,murre_frame_rate2,...
    Gregor_frame_rate*ones(1,5), Gregor_frame_rate2,Gregor_frame_rate*ones(1,4),...
    .5*ones(1,22),  .5*ones(1,8)];

all_nCells = [pub_nCells(pub_nCells~=0 & ~isnan(pub_nCells)); exp_nCells(:);  exp2_nCells(:)];
all_nObs = [pub_nObs(pub_nObs~=0 & ~isnan(pub_nObs)); exp_nObs(:); exp2_nObs(:)];

pub_stp = pub_FnD(:,:,2)';
exp_stp = exp_FnD(:,:,2)';
exp2_stp = exp2_FnD(:,:,2)';
all_D = [pub_stp(pub_stp~=0 & ~isnan(pub_stp));exp_stp(:);exp2_stp(:)];
pub_stp = pub_FnD(:,:,1)';
exp_stp = exp_FnD(:,:,1)';
exp2_stp = exp_FnD(:,:,1)';
all_F = [pub_stp(pub_stp~=0 & ~isnan(pub_stp)); exp_stp(:); exp2_stp(:)];

figure(103); clf;
cmapDist = hsv(ceil(log2(100e3))-2 + 2);
pub_ColorID = ceil(log2(all_GenDist) - 2);
for e=1:nE
    figure(103);
    subplot(6,1,1); 
    bar(e,all_D(e),'FaceColor',cmapDist(pub_ColorID(e),:)); hold on;
    fc =  .5*cmapDist(pub_ColorID(e),:);
    bar(e,all_F(e),'FaceColor',fc); hold on;
    ylabel('distance (nm)'); title('median frame step-size and dynamic range (median IQR)'); box off;

    figure(103);
    subplot(6,1,2);
    bar(e,all_F(e)./all_D(e),'FaceColor',cmapDist(pub_ColorID(e),:)); hold on;
    ylabel('F-to-D ratio'); title('F-to-D ratio'); box off;


    subplot(6,1,3); 
    bar(e,1./timeRes(e),'FaceColor',cmapDist(pub_ColorID(e),:)); hold on;
    ylabel('Stacks per second'); title('imaging speed'); box off;

      subplot(6,1,4); 
    bar(e,all_nCells(e),'FaceColor',cmapDist(pub_ColorID(e),:)); hold on;
    ylabel('number of traces'); title('number of cells traced'); box off;

    
      subplot(6,1,5); 
    bar(e,all_nObs(e),'FaceColor',cmapDist(pub_ColorID(e),:)); hold on;
    ylabel('number of Obs'); title('total number of observations'); box off;

        subplot(6,1,6);
    semilogy(1,1,'w.'); hold on;
    bar(e,all_GenDist(e),'FaceColor',cmapDist(pub_ColorID(e),:)); hold on;
    ylabel('genomic distance (kb)');  box off; title('probe separation')


end

 subplot(6,1,6);
studyBar = [1,sum(~isnan(pub_stp))];
studyClr = lines(5); 
for s=1:5
plot(  [sum(studyBar(1:s)),sum(studyBar(1:s+1))],[1e5,1e5],'-','color',studyClr(s,:),'linewidth',3);
end



%% Tot obs
allDat1 = cat(1,data3D{:,:});
domain1_N = sum(~isnan(allDat1(:)))
allDat = cat(1,data3D{1,:});
n_WT11lines = sum(~isnan(allDat(:)))
allDat = cat(1,data3D{2,:});
n_dTag11lines = sum(~isnan(allDat(:)))

allDat2 = cat(1,data3D_2{:,:});
domain2_N = sum(~isnan(allDat2(:)))

domain1_N + domain2_N
domain1_C = size(allDat1,1)
domain2_C = size(allDat2,1)

tObs = zeros(1,5);
tCells = zeros(1,5);
for p=1:5
    pData = cellfun(@(x) sum(~isnan(x(:))),papers{p}) ; %  cat(1,papers{p}{:});
    cData = cellfun(@(x) size(x,1),papers{p}) ;
    tObs(p) = sum(pData);
    tCells(p) = sum(cData);
end

figure(1); clf;
subplot(1,2,1);
bar([tCells,domain1_C+domain2_C]); box off;
text((1:6)-.5,[tCells,domain1_C+domain2_C]+.1e4,cellstr(num2str([tCells,domain1_C+domain2_C]',2)))
set(gca,'XTick',1:6,'XTickLabel',paperNames([1:5,8])); ylabel('Tot. Cells'); 

subplot(1,2,2);
bar([tObs,domain1_N+domain2_N]); box off;
text((1:6)-.5,[tObs,domain1_N+domain2_N]+.1e7,cellstr(num2str([tObs,domain1_N+domain2_N]',2)))
set(gca,'XTick',1:6,'XTickLabel',paperNames([1:5,8]));
ylabel('Tot. Obs.'); 


%% original 11 lines
T = 3600;
cmap = hsv(12);
ObsFracE = zeros(11,T);
figure(1); clf;
for e=1:11
    ee = 12-e;
    dis3D = cat(1,data3D{1,ee}(:,1:T),data3D{2,ee}(:,1:T)); % combine same cell line 
    nObs = sum(~isnan(dis3D),2);
    [nObs,idx] = sort(nObs,'descend');
    dis3D = dis3D(idx,:);
    nCells = length(nObs);
    nMax = min([nCells,100]); % first 100 lines or max, to avoid bias by  
    dis3D = dis3D(1:nMax,:);
    nC = size(dis3D,1);
    ObsFrac = sum(~isnan(dis3D),1)/nC;
    ObsFrac = ObsFrac./max(ObsFrac);
    ObsFracE(e,:) = ObsFrac;
    plot(ObsFrac,'color',cmap(ee,:)); hold on;
    xlim([300,3600]); ylim([0,1]); box off;
end
legend('location','southwest')
ylabel('fraction detected')
xlabel('frame number')


figure(3); clf; 
for p=1:5
    nPE = length(papers{p});
        cmap2 = hsv(nPE+2);
    for e=1:nPE
        subplot(1,5,p);
        dis3D =  papers{p}{e};
        dis3D(dis3D==0) = nan;
        nObs = sum(~isnan(dis3D),2);
        [nObs,idx] = sort(nObs,'descend');
        dis3D = dis3D(idx,:);
        nCells = length(nObs);
        nMax = min([nCells,100]);
        dis3D = dis3D(1:nMax,:);
        nC = size(dis3D,1);
        ObsFrac = sum(~isnan(dis3D),1)/nC;
        ObsFrac = ObsFrac./max(ObsFrac);
        ObsFrac(ObsFrac==0)= NaN;
        figure(3); % subplot(11,1,e); 
        plot(ObsFrac,'color',cmap2(e,:)); hold on;
        % xlim([300,3600]); 
        ylim([0,1]); box off;
    end
     legend(expNames{p},'location','eastoutside','FontSize',5)
    ylabel('fraction detected')
    xlabel('frame number')
end



figure(2); clf; cmap2 = hsv(7);
for e=1:4
    dis3D = cat(1,data3D_2{1,e}(:,1:T),data3D_2{2,e}(:,1:T)); % combine same cell line 
    nObs = sum(~isnan(dis3D),2);
    [nObs,idx] = sort(nObs,'descend');
    dis3D = dis3D(idx,:);
    nCells = length(nObs);
    nMax = min([nCells,100]);
    dis3D = dis3D(1:nMax,:);
    nC = size(dis3D,1);
    ObsFrac = sum(~isnan(dis3D),1)/nC;
    ObsFrac = ObsFrac./max(ObsFrac);
    ObsFrac(ObsFrac==0)= NaN;
    figure(2); % subplot(11,1,e); 
    plot(ObsFrac,'color',cmap2(2+e,:)); hold on;
    xlim([300,3600]); ylim([0,1]); box off;
end
legend('location','southwest')
ylabel('fraction detected')
xlabel('frame number')


