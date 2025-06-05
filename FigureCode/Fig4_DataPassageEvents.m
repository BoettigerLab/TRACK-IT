%% Fig 4, plot experimental data

%% load data
% dataset 1
dataFolder = 'U:\Manuscripts\Jude Live Imaging\Data\traces_v13\';
tic
dis_kb = [5,20,55,70,134,260,407,799,2030,12795,73567];
data3D = cell(2,11);
for e=1:11 % e=11    
    datafile = [dataFolder,'Dis3D_2Hz_noTreatment_',num2str(dis_kb(e)),'kb.csv'];
    data3D{1,e} = readmatrix(datafile);
    datafile = [dataFolder,'Dis3D_2Hz_dTag_',num2str(dis_kb(e)),'kb.csv'];
    data3D{2,e} = readmatrix(datafile);
end


%%

dTagStatus = {'-dTag','+dTag'};

es = 4:8;
N = length(es);
procIdStartEnd = cell(2,N);
procTraceStats = cell(2,N);
extSpeeds = cell(2,N);
procTraces = cell(2,N);
figure(4); clf;
dd = 0;
for i=1:N
    for d=1:2
        dd=dd+1;
        e=es(i); 
        disM = data3D{d,e}; 
        [traces,traceIdStartEnd,traceStats] = FindProcessive(disM,'minD',2*nanmedian(disM(:)),...
             'minLen',16,'maxBackStep',.26,'maxBackDist',.26,'maxMiss',.2,'maxLen',35,'showPlots',0,'maxD',5*nanmedian(disM(:)));
        if isempty(traceStats)
            continue
        end
        procIdStartEnd{d,i} = traceIdStartEnd;
        procTraceStats{d,i} =traceStats;
        procTraces{d,i} = traces;
        extSpeeds{d,i} = [traceStats.passageSize]*40/1000./[traceStats.totSteps]*(2/1);  % (nm/fr)*(fr/s)
    end
end


%% plot
figure(1); clf;
lineWidth = 1;
markerSize =6;
showTrace = {[],[3,4,1,7,2],[1,2,3,4,5],[1,2,3,4,5],[]; 
             [], 1,          1,     [],   [] };

% showTrace = {[],[3,1],[1,5],[5,2,4],[]; 
%              [], 1,          1,     [],   [] };

s = 0;
for i=2:N-1
    s=s+1;
    for d=1:2
         e=es(i);
        disM = data3D{d,e};
        clr_pre_nLE = [1 .8 .8];
        clr_pre_yLE = [.8 .8 1];
        clr_pas_yLE = [.1 .6 1];
        clr_pas_nLE = [1 .6 .1];
        
        % Plot traces
        figure(1);  % subplot(N,2,(i-1)*2+1); % cla;
        nT = size(disM,2);
        traceIdStartEnd = procIdStartEnd{d,i}; 
        nC = size(traceIdStartEnd,1);
        cmap = hsv(nC+1);
        k=0;
        pre = 10;
        win = 60;
        post = win -pre; 
        lg = [];
        names = {};
        for c = showTrace{d,i} %  1:nC    % c =5
            figure(4); subplot(1,N-2,s);  % subplot(N,2,(i-1)*2+1);
            t0 = max([1,traceIdStartEnd(c,2)-pre+1]);
            t1 = min([nT,traceIdStartEnd(c,2)+post]);
            yy = disM(traceIdStartEnd(c,1),t0:t1);
            yy = yy - min(yy); 
            t = 1:win;
            if t0 < 2 % pad to keep the alignment
                pad = pre - traceIdStartEnd(c,2);
                yy = [nan(1,pad),yy];
            end

            if isnan(yy(1)); y00 = yy(~isnan(yy)); yy(1) = y00(1); end;

            tout = t(~isnan(yy));
            if tout(1)~=pre
                plot(tout,yy(~isnan(yy)),'.-','color',clr_pre_yLE,'linewidth',lineWidth,'markerSize',markerSize);  hold on;   
                t_ext = traceIdStartEnd(c,3)-traceIdStartEnd(c,2);
                t = t(pre:pre+t_ext);  yy = yy(pre:pre+t_ext);
                lg(c) = plot(t(~isnan(yy)),yy(~isnan(yy)),'.-','linewidth',lineWidth,'markerSize',markerSize,'color',clr_pas_yLE);  hold on; 
                names{c} = num2str(c);
                ylim([0,1200])
            pause(.1);
            end
        end
        box off; 
        xlabel('frame (0.5s)'); 
        ylabel('3D distance change (nm)')
        title([ num2str(dis_kb(e)) ])
    end
end
%% extrusion speed

figure(1); clf; 
yyaxis left;
semilogy(1,1,'w.'); hold on;
violin(extSpeeds(1,2:4),'bandwidth',.1,'plotMean',false,'lineColor','none','faceColor',.7*ones(1,3));
hold on; BoxPlotCell(extSpeeds(1,2:4));
for i=2:4
    sp = extSpeeds{1,i};
    x = i-1.5+ .5*ones(length(sp),1) +  .3*(0.5-rand(length(sp),1));
    plot(x,sp,'k.','markerSize',10);
end
ylim([.1,10]);
box off;
aveExt = nanmean(cat(2,extSpeeds{1,2:4}));
medExt = nanmedian(cat(2,extSpeeds{1,2:4}))
title(['Ave ext speed =',num2str(aveExt,2)]);
ylabel('extrusion speed (kb/s)');
set(gca,'XtickLabel',cellstr(num2str(dis_kb(es(2:4))')))

yyaxis right;
set(gca,'YScale','log');
ylim(1000/40*[.1,10]);
ylabel('extrusion speed (nm/s)')
numProc = cellfun(@length,extSpeeds);
numProc(1,5)=numProc(1,5)-1;  numProc(2,2)=numProc(2,2)-1; % no pre-frames to validate trace

% chi-square 
n1 = sum(numProc(1,2:4));
n2 = sum(numProc(2,2:4));
N1 = sum(cellfun(@(x) size(x,1), data3D(1,es(2:4))));
N2 = sum(cellfun(@(x) size(x,1), data3D(2,es(2:4))));
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

figure(5); clf;
imagesc(numProc); colorbar; GetColorMap('CyanToBlack');
set(gca,'YTick',1:2,'YTickLabel',{'untreated','+dTag'},...
    'XTick',1:5,'XTickLabel',dis_kb(es),'FontSize',10)
box off;

title(['number of processive traces. \chi^2 p-value=',num2str(pval,2)],'FontSize',10)

figure(5); clf;
bar(1:N,numProc); ylabel('processive traces');
set(gca,'XTick',1:5,'XTickLabel',dis_kb(es),'FontSize',10); box off;
xlabel('separation'); legend({'untreated','+dTag'})
title(['number of processive traces. \chi^2 p-value=',num2str(pval,2)],'FontSize',10)
    
%%
mean_speed = cellfun(@mean,extSpeeds(1,2:4))
std_speed = cellfun(@std,extSpeeds(1,2:4))
n_speeds = cellfun(@length,extSpeeds(1,2:4))
SEM_speed = std_speed./sqrt(n_speeds)

all_speed = cat(2,extSpeeds{1,2:4});
mean_all = mean(all_speed )
std_all = std( all_speed)/sqrt(length(all_speed))
