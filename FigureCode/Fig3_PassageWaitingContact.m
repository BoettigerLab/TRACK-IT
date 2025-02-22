%% includes domain 1 and domain 2
% includes supp Figs on threshold effects and G1 vs G2

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

% load exp 2
dis_kb2 = [100,100,400,400];
dis_kb2_sign =  [95,-108,-409,407];
tadPos ={'intra','inter','inter','intra'};
data3D_2 = cell(2,4); 
for e=1:4    
    datafile = [dataFolder,'Dis3D_2Hz_Reg2_noTreatment_',tadPos{e},'_',num2str(dis_kb2(e)),'kb.csv'];
    data3D_2{1,e} = readmatrix(datafile);
    datafile = [dataFolder,'Dis3D_2Hz_Reg2_dTag_',tadPos{e},'_',num2str(dis_kb2(e)),'kb.csv'];
    data3D_2{2,e} = readmatrix(datafile);
end

toc

%%  Compute
theta = 50;

[nt_wt_50,   nt_ci_wt50, nt_cdfs_50] = CompContactWaitingKM(data3D(1,:),'d_contact',theta,'cI',.01);
[dTag_wt_50, dTag_ci_wt50, dTag_cdfs_50] = CompContactWaitingKM(data3D(2,:),'d_contact',theta,'cI',.01);

[p3_freq_nt,p3_freq_nt_ci] = CompPassageFreq(data3D(1,:),'d_contact',theta);
[p3_freq_dTag,p3_freq_dTag_ci] = CompPassageFreq(data3D(2,:),'d_contact',theta);

[cd_nt,cd_nt_ci,cd_all_nt] = ContactDurationKM(data3D(1,:),'theta',theta);
[cd_dTag,cd_dTag_ci,cd_all_dTag] = ContactDurationKM(data3D(2,:),'theta',theta);

[nt_wt_500,   nt_ci_wt500] = CompContactWaitingKM(data3D(1,:),'d_contact',500);
[dTag_wt_500, dTag_ci_wt500] = CompContactWaitingKM(data3D(2,:),'d_contact',500);

%% Plot

cmap = hsv(12);
cmap_dtag = .3*cmap + .7;
gry1 = [.4 .4 .4];
gry2 = [.8 .8 .8];
spf = 0.5; % seconds per frame 
tb = 600; % TAD border
tc = [.75 0 0]; % TAD border line color
figure(1); clf;

subplot(1,4,1);
sel = 1:8;
x = dis_kb(sel); 
y1=spf*nt_wt_50(sel,3);
xe = [dis_kb(sel); dis_kb(sel)];
ye1 = spf*nt_ci_wt50(sel,:,3)';
y2 = spf*dTag_wt_50(sel,3);
ye2 = spf*dTag_ci_wt50(sel,:,3)';
clr1 = cmap(sel,:);
clr2 = cmap_dtag(sel,:);
for i=1:length(sel)
    loglog(x(i),y1(i),'o','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1(i),'o','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye1(:,i),'k-');
    loglog(x(i),y2(i),'^','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
     loglog(x(i),y2(i),'^','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2(:,i),'k-');
end
sc = max([nt_wt_50(8,3)./nt_wt_500(8,1),dTag_wt_50(8,3)./dTag_wt_500(8,1)]);
sel = 8:11; 
x = dis_kb(sel); 
y1= sc*spf*nt_wt_500(sel,1);
xe = [dis_kb(sel); dis_kb(sel)];
ye1 = sc*spf*nt_ci_wt500(sel,:,1)';
y2 = sc*spf*dTag_wt_500(sel,1);
ye2 = sc*spf*dTag_ci_wt500(sel,:,1)';
clr1 = cmap(sel,:);
clr2 = cmap_dtag(sel,:);
for i=1:length(sel)
    loglog(x(i),y1(i),'o','MarkerSize',2,'color',clr1(i,:),'linewidth',7); hold on;
     loglog(x(i),y1(i),'o','MarkerSize',7,'color','k','linewidth',1); hold on;
    plot(xe(:,i),ye1(:,i),'k-');
    loglog(x(i),y2(i),'^','MarkerSize',2,'color',clr2(i,:),'linewidth',5); hold on;
     loglog(x(i),y2(i),'^','MarkerSize',6,'color','k','linewidth',1); hold on;
    plot(xe(:,i),ye2(:,i),'k-');
end
plot([tb,tb],[1,5e6],'--','color',tc);
rectangle('Position',[14,6,1.2e3,.2e4],'edgecolor',.7*ones(1,3));
box off; 
xlim([4,1e5]);  ylim([4,.8e6]);
xlabel('genome separation (kb)');
ylabel('search time (s)');
title('search time')

%
subplot(1,4,2);
sel = 1:8;
x = dis_kb(sel); 
y1=spf*nt_wt_50(sel,3);
xe = [dis_kb(sel); dis_kb(sel)];
ye1 = spf*nt_ci_wt50(sel,:,3)';
y2 = spf*dTag_wt_50(sel,3);
ye2 = spf*dTag_ci_wt50(sel,:,3)';
clr1 = cmap(sel,:);
clr2 = cmap_dtag(sel,:);
for i=1:length(sel)
    loglog(x(i),y1(i),'o','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1(i),'o','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye1(:,i),'k-');
    loglog(x(i),y2(i),'^','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
    loglog(x(i),y2(i),'^','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2(:,i),'k-');
end
box off;
xlim([4,2e3]);  ylim([4,.3e4]);
plot([tb,tb],[1,5e6],'--','color',tc);
xlabel('genome separation (kb)');
ylabel('search time (s)');
title('search time')

%
subplot(1,4,3);
fph = 2*60*60;  % frames to hours
x = dis_kb'; 
y1 = fph*p3_freq_nt;
xe = [dis_kb; dis_kb];
ye1 = fph*p3_freq_nt_ci';
y2 = fph*p3_freq_dTag;
ye2 = fph*p3_freq_dTag_ci';
clr1 = cmap;
clr2 = cmap_dtag;
for i=1:length(x)
    loglog(x(i),y1(i),'o','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1(i),'o','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye1(:,i),'k-');
    loglog(x(i),y2(i),'^','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
    loglog(x(i),y2(i),'^','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2(:,i),'k-');
end
 plot([tb,tb],[1e-2,1e2],'--','color',tc);
 xlim([4,1e5]);
box off;
xlabel('genomic separation (kb)');
ylabel('search frequency (events/hr)');
set(gcf,'color','w');
title('search frequency')

% 
subplot(1,4,4);
x = dis_kb'; 
y1 = spf*cd_nt(:,:,2);
xe = [dis_kb; dis_kb];
ye1 = spf*cd_nt_ci(:,:,:,2);
y2 = spf*cd_dTag(:,:,2);
ye2 = spf*cd_dTag_ci(:,:,:,2);
clr1 = cmap;
clr2 = cmap_dtag;
sz = 10;
for i=1:length(x)
    loglog(x(i),y1(i),'o','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1(i),'o','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye1(:,i),'k-');
    loglog(x(i),y2(i),'^','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
    loglog(x(i),y2(i),'^','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2(:,i),'k-');
end
plot([tb,tb],[.1,10],'--','color',tc);
box off;
xlabel('genomic separation (kb)');
ylabel('contact duration (s)');
ylim([.4 0.8]);
xlim([4,1e5]); 
title('contact duration')



figure(1); subplot(1,4,1);
  fit_dTAG = QuickLogFit(dis_kb(2:10),[spf*dTag_wt_50(2:8,3); sc*spf*nt_wt_500(9:10,1)]);
 fit_NTnear = QuickLogFit(dis_kb(2:7),[spf*nt_wt_50(2:7,3)],'plotOpts',{'r-'});
y0 = .06; % ref lines
y1 = .03;
plot(dis_kb(2:11),y1*dis_kb(2:11).^(2),'--','color',gry2);
plot(dis_kb(2:11),y0*dis_kb(2:11).^(5/3),'--','color',gry1);

subplot(1,4,2);
plot(dis_kb(2:8),y1*dis_kb(2:8).^(2),'--','color',gry2);
plot(dis_kb(2:8),y0*dis_kb(2:8).^(5/3),'--','color',gry1);
plot(dis_kb(2:8),10^fit_dTAG.p2*dis_kb(2:8).^fit_dTAG.p1,'b-')
plot(dis_kb(2:8),10^fit_NTnear.p2*dis_kb(2:8).^fit_NTnear.p1,'r-')


%% legend
figure(3); clf;
for i=1:11
    plot(i,0,'o','color',cmap(i,:),'MarkerSize',2,'linewidth',10); hold on;
    plot(i,0,'o','color','r','MarkerSize',10); hold on;
    
    plot(i,1,'o','color',cmap_dtag(i,:),'MarkerSize',2,'linewidth',10); hold on;
    plot(i,1,'o','color','b','MarkerSize',10); hold on;
end
ylim([-1 2]); xlim([0,12]); box off;

%% corrected vs. uncorrected

figure(2); clf;
sel = 1:8;
x = dis_kb(sel); 
y1=spf*nt_wt_50(sel,3);
y1r=spf*nt_wt_50(sel,1);
xe = [dis_kb(sel); dis_kb(sel)];
ye1 = spf*nt_ci_wt50(sel,:,3)';
ye1r = spf*nt_ci_wt50(sel,:,1)';
y2 = spf*dTag_wt_50(sel,3);
y2r = spf*dTag_wt_50(sel,1);
ye2 = spf*dTag_ci_wt50(sel,:,3)';
ye2r = spf*dTag_ci_wt50(sel,:,1)';
clr1 = cmap(sel,:);
clr2 = cmap_dtag(sel,:);
for i=1:length(sel)
    subplot(1,3,2);
    loglog(x(i),y1(i),'o','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1(i),'o','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye1(:,i),'k-');

     loglog(x(i),y1r(i),'s','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1r(i),'s','MarkerSize',13,'color','k'); hold on;
    plot(xe(:,i),ye1r(:,i),'k-');

    subplot(1,3,3);
    loglog(x(i),y2(i),'^','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
    loglog(x(i),y2(i),'^','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2(:,i),'k-');

    loglog(x(i),y2r(i),'v','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
    loglog(x(i),y2r(i),'v','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2r(:,i),'k-');
end


for i=1:length(sel)
    subplot(1,3,1);

    loglog(x(i),y1r(i),'s','MarkerSize',2,'color',clr1(i,:),'linewidth',10); hold on;
    loglog(x(i),y1r(i),'s','MarkerSize',13,'color','k'); hold on;
    plot(xe(:,i),ye1r(:,i),'k-');

    loglog(x(i),y2r(i),'v','MarkerSize',2,'color',clr2(i,:),'linewidth',7); hold on;
    loglog(x(i),y2r(i),'v','MarkerSize',11,'color','k'); hold on;
    plot(xe(:,i),ye2r(:,i),'k-');
end


for s=1:3
    subplot(1,3,s);
    box off;
    xlim([4,2e3]);  ylim([4,.3e4]);
    plot([tb,tb],[1,5e6],'--','color',tc);
    xlabel('genome separation (kb)');
    ylabel('search time (s)');
    title('search time')
end
subplot(1,3,1); ylim([4,.2e3]);

%% legend
figure(3); clf;
for i=1:11
    plot(i,0,'o','color',cmap(i,:),'MarkerSize',2,'linewidth',10); hold on;
     plot(i,0,'o','color','k','MarkerSize',11); hold on;
    
    plot(i,1,'^','color',cmap_dtag(i,:),'MarkerSize',2,'linewidth',7); hold on;
     plot(i,1,'^','color','k','MarkerSize',11); hold on;

       plot(i,2,'s','color',cmap(i,:),'MarkerSize',2,'linewidth',10); hold on;
     plot(i,2,'s','color','k','MarkerSize',13); hold on;

         plot(i,3,'v','color',cmap_dtag(i,:),'MarkerSize',2,'linewidth',7); hold on;
     plot(i,3,'v','color','k','MarkerSize',11); hold on;
end
ylim([-1 4]); xlim([0,12]); box off;
ylim([-1 4]); xlim([0,12]); box off;

%% domain 2 computations
theta = 50;

[nt2_wt_50,   nt2_ci_wt50] = CompContactWaitingKM(data3D_2(1,:),'d_contact',theta);
[dTag2_wt_50, dTag2_ci_wt50] = CompContactWaitingKM(data3D_2(2,:),'d_contact',theta);

[pt2_freq_nt,pt2_freq_nt_ci] = CompPassageFreq(data3D_2(1,:),'d_contact',theta);
[pt2_freq_dTag,pt2_freq_dTag_ci] = CompPassageFreq(data3D_2(2,:),'d_contact',theta);

[cd2_nt,cd2_nt_ci] = ContactDurationKM(data3D_2(1,:),'theta',theta);
[cd2_dTag,cd2_dTag_ci] = ContactDurationKM(data3D_2(2,:),'theta',theta);



%% plot results
figure(6); clf;
subplot(1,3,1); 
semilogy(dis_kb2_sign,spf*nt2_wt_50(:,3),'ro','MarkerSize',8); hold on;
plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_wt50(:,:,3)','k-'); hold on;
plot([-50,-50],[10,1e3],'k--');
box off; xlim([-500 500]); ylim([10,1e3]); ylabel('search time (s)');

subplot(1,3,2); 
semilogy(dis_kb2_sign,fph*pt2_freq_nt,'ro','MarkerSize',8); hold on;
plot([dis_kb2_sign;dis_kb2_sign],fph*pt2_freq_nt_ci','k-'); hold on;
plot([-50,-50],[1,20],'k--');
box off; xlim([-500 500]);   ylim([3.,15]);
ylabel('passage freq. (events/hr)');

subplot(1,3,3); 
semilogy(dis_kb2_sign,spf*cd2_nt(:,:,2),'ro','MarkerSize',8); hold on;
plot([dis_kb2_sign;dis_kb2_sign],spf*cd2_nt_ci(:,:,:,2),'k-'); hold on;
plot([-50,-50],[.3,1.4],'k--');
box off; xlim([-500 500]);
ylim([.4,.8])
ylabel('ave. contact duration (s)');

%%  SUPPLEMENT, alternate thresholds

figure(10); clf;
thetas = [30, 50, 100, 200, 350 500];
k=0;
N=length(thetas);
maxSamples =1e3;
nt_wt_T  =cell(N,1);  
nt_ci_wt_T =cell(N,1);
dTag_wt_T =cell(N,1); 
dTag_ci_wt_T =cell(N,1);   
nt_pt_T =cell(N,1);   
nt_ci_pt_T =cell(N,1); 
dTag_pt_T =cell(N,1); 
dTag_ci_pt_T =cell(N,1);    
p3_freq_nt_T =cell(N,1);  
p3_freq_nt_ci_T   =cell(N,1);  
p3_freq_dTag_T =cell(N,1);  
p3_freq_dTag_ci_T =cell(N,1);      
cd_nt_T =cell(N,1);    
cd_nt_ci_T =cell(N,1);  
cd_dTag_T =cell(N,1);   
cd_dTag_ci_T =cell(N,1);  
for t=1:N
    theta = thetas(t);    
    [nt_wt_T{t},   nt_ci_wt_T{t}] = CompContactWaitingKM(data3D(1,:),'d_contact',theta,'maxSamples',maxSamples);
    [dTag_wt_T{t}, dTag_ci_wt_T{t}] = CompContactWaitingKM(data3D(2,:),'d_contact',theta,'maxSamples',maxSamples);
    
    [nt_pt_T{t},   nt_ci_pt_T{t}] = CompFirstPassageKM(data3D(1,:),'d_contact',theta,'maxSamples',maxSamples);
    [dTag_pt_T{t}, dTag_ci_pt_T{t}] = CompFirstPassageKM(data3D(2,:),'d_contact',theta,'maxSamples',maxSamples);
        
    [cd_nt_T{t},cd_nt_ci_T{t}] = ContactDurationKM(data3D(1,:),'theta',theta);
    [cd_dTag_T{t},cd_dTag_ci_T{t}] = ContactDurationKM(data3D(2,:),'theta',theta);   
end

%% Plot threshold results
figure(10); clf;
for t=1:N
    gry1 = [.4 .4 .4];
    gry2 = [.8 .8 .8];
    sel = 1:11;
    spf = 0.5; % seconds per frame 
    
    tb = 600;
    if thetas(t)<300; k=3; else k=1; end
    subplot(4,N,t);
    loglog(dis_kb(sel),spf*nt_wt_T{t}(sel,k),'r^','MarkerSize',4); hold on;
    loglog(dis_kb(sel),spf*dTag_wt_T{t}(sel,k),'bo','MarkerSize',4); hold on;
    if k==3
    plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_wt_T{t}(sel,:)','r-');
    plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_wt_T{t}(sel,:)','b-');
    end
    plot([tb,tb],[1,5e4],'k--');
    max(spf*dTag_wt_T{t}(sel,3))
    box off; 
    axis tight;
     xlim([4,1e5]);   ylim([1,.5e4]); set(gca,'YTick',logspace(0,3,4),'YTickLabel',logspace(0,3,4))
    xlabel('genome separation (kb)');
    ylabel(['waiting time (s)']);
    title(['search time ',num2str(thetas(t)), ' nm'])
    % 

    k=k+1; subplot(4,N,1*N+t);
    fph = 2*60*60;  % frames to hours
    loglog(dis_kb',fph*p3_freq_nt_T{t},'r^','MarkerSize',4); hold on;
    plot([dis_kb;dis_kb],fph*p3_freq_nt_ci_T{t}','r-'); hold on;
    plot(dis_kb',fph*p3_freq_dTag_T{t},'bo','MarkerSize',4); hold on;
    plot([dis_kb;dis_kb],fph*p3_freq_dTag_ci_T{t}','b-'); hold on;
     plot([tb,tb],[1e-2,1e2],'k--');
     xlim([4,1e5]);
    box off;
    xlabel('genomic separation (kb)');
    ylabel('passage frequency (events/hr)');
    set(gcf,'color','w');
    title('passage frequency')
    % title(['contact threshold=',num2str(theta)]);
    
     k=k+1; subplot(4,N,2*N+t);
    loglog(dis_kb',spf*cd_nt_T{t},'r^','MarkerSize',4); hold on;
    plot([dis_kb;dis_kb],spf*cd_nt_ci_T{t},'r-'); hold on;
    plot(dis_kb',spf*cd_dTag_T{t},'bo','MarkerSize',4); hold on;
    plot([dis_kb;dis_kb],spf*cd_dTag_ci_T{t},'b-'); hold on;
    plot([tb,tb],[.1,100],'k--');
    box off;
    xlabel('genomic separation (kb)');
    ylabel('contact duration (s)');
    ylim([.7 2*t^2]);
    title('contact duration')
    % title(['contact threshold=',num2str(theta)]);
    
    subplot(4,N,3*N+t);
    if thetas(t)<300; k=3; else k=1; end;
    loglog(dis_kb(sel),spf*nt_pt_T{t}(sel,3),'r^','MarkerSize',4); hold on;
    loglog(dis_kb(sel),spf*dTag_pt_T{t}(sel,3),'bo','MarkerSize',4); hold on;
    if k==3
        plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_pt_T{t}(sel,:)','r-');
        plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_pt_T{t}(sel,:)','b-');
    end
    plot([tb,tb],[1,5e4],'k--');
    box off; 
    xlim([4,1e5]);   ylim([1,.5e4]); set(gca,'YTick',logspace(0,3,4),'YTickLabel',logspace(0,3,4))
    xlabel('genome distance (kb)');
    ylabel('passage time (s)')
    title('passage time')
end

%% G1 vs G2   



data3D_g1 = data3D;
data3D_g2 = data3D;
for i=1:numel(data3D_g1)
    g1 = data3D_g1{i}(:,end) < 5;
    data3D_g1{i} = data3D_g1{i}(g1,:);
    data3D_g2{i} = data3D_g2{i}(~g1,:);
end

%% 
theta = 50;

[nt_wt_50_g1,   nt_ci_wt50_g1] = CompContactWaitingKM(data3D_g1(1,:),'d_contact',theta);
[dTag_wt_50_g1, dTag_ci_wt50_g1] = CompContactWaitingKM(data3D_g1(2,:),'d_contact',theta);

[nt_pt_50_g1,   nt_ci_pt50_g1] = CompFirstPassageKM(data3D_g1(1,:),'d_contact',theta);
[dTag_pt_50_g1, dTag_ci_pt50_g1] = CompFirstPassageKM(data3D_g1(2,:),'d_contact',theta);

[p3_freq_nt_g1,p3_freq_nt_ci_g1] = CompPassageFreq(data3D_g1(1,:),'d_contact',theta);
[p3_freq_dTag_g1,p3_freq_dTag_ci_g1] = CompPassageFreq(data3D_g1(2,:),'d_contact',theta);

[cd_nt_g1,cd_nt_ci_g1] = ContactDuration(data3D_g1(1,:),'theta',theta);
[cd_dTag_g1,cd_dTag_ci_g1] = ContactDuration(data3D_g1(2,:),'theta',theta);

[nt_wt_500_g1,   nt_ci_wt500_g1] = CompContactWaitingKM(data3D_g1(1,:),'d_contact',500);
[dTag_wt_500_g1, dTag_ci_wt500_g1] = CompContactWaitingKM(data3D_g1(2,:),'d_contact',500);

[nt_pt_500_g1,   nt_ci_pt500_g1] = CompFirstPassageKM(data3D_g1(1,:),'d_contact',500);
[dTag_pt_500_g1, dTag_ci_pt500_g1] = CompFirstPassageKM(data3D_g1(2,:),'d_contact',500);

%%
gry1 = [.4 .4 .4];
gry2 = [.8 .8 .8];
sel = 1:8;
spf = 0.5; % seconds per frame 
tb = 600; % TAD border
tc = [.75 0 0]; % TAD border line color
figure(1); clf;

subplot(2,4,1);
loglog(dis_kb(sel),spf*nt_wt_50_g1(sel,3),'r^','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_wt50_g1(sel,:,3)','r-');
loglog(dis_kb(sel),spf*dTag_wt_50_g1(sel,3),'bo','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_wt50_g1(sel,:,3)','b-');
sc = max([nt_wt_50_g1(8,3)./nt_wt_500_g1(8,1),dTag_wt_50_g1(8,3)./dTag_wt_500_g1(8,1)]);
sel = 8:11; 
loglog(dis_kb(sel),sc*spf*nt_wt_500_g1(sel,1),'r>','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],sc*spf*nt_ci_wt500_g1(sel,:,1)','r-');
loglog(dis_kb(sel),sc*spf*dTag_wt_500_g1(sel,1),'bs','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],sc*spf*dTag_ci_wt500_g1(sel,:,1)','b-');
plot([tb,tb],[1,5e6],'--','color',tc);
xx = logspace(1,5,5);
rectangle('Position',[14,6,1.2e3,.2e4],'edgecolor',.7*ones(1,3));
plot(xx,.064*xx.^(2),'--','color',gry2);
plot(xx,.19*xx.^(5/3),'--','color',gry1);
box off; 
xlim([4,1e5]);  ylim([4,.8e6]);
xlabel('genome distance (kb)');
ylabel('search time (s)');
title('G1 search time')

sel = 1:8;
tb = 600; % TAD border
subplot(2,4,2);
loglog(dis_kb(sel),spf*nt_wt_50_g1(sel,3),'r^','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_wt50_g1(sel,:,3)','r-');
loglog(dis_kb(sel),spf*dTag_wt_50_g1(sel,3),'bo','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_wt50_g1(sel,:,3)','b-');
plot([tb,tb],[1,5e6],'--','color',tc);
xx = logspace(1,5,5);
plot(xx,.064*xx.^(2),'--','color',gry2);
plot(xx,.19*xx.^(5/3),'--','color',gry1);
box off; 
xlim([14,1.2e3]);  ylim([6,.2e4]);
xlabel('genome distance (kb)');
ylabel('search time (s)');
title('G1 search time')

subplot(2,4,3);
fph = 2*60*60;  % frames to hours
loglog(dis_kb',fph*p3_freq_nt_g1,'r^','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],fph*p3_freq_nt_ci_g1','r-'); hold on;
plot(dis_kb',fph*p3_freq_dTag_g1,'bo','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],fph*p3_freq_dTag_ci_g1','b-'); hold on;
 plot([tb,tb],[1e-2,1e2],'--','color',tc);
 xlim([4,1e5]);
box off;
xlabel('genomic separation (kb)');
ylabel('passage frequency (events/hr)');
set(gcf,'color','w');
title('G1 passage frequency')
% title(['contact threshold=',num2str(theta)]);

subplot(2,4,4);
loglog(dis_kb',cd_nt_g1,'r^','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],cd_nt_ci_g1,'r-'); hold on;
plot(dis_kb',cd_dTag_g1,'bo','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],cd_dTag_ci_g1,'b-'); hold on;
plot([tb,tb],[.1,10],'--','color',tc);
box off;
xlabel('genomic separation (kb)');
ylabel('contact duration (s)');
ylim([.9 1.4]);
title('G1 contact duration')
% title(['contact threshold=',num2str(theta)]);




sel =8:11;
figure(2); clf; 
loglog(dis_kb(sel),spf*nt_wt_500_g1(sel,1),'r^','MarkerSize',4); hold on
% plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_wt500(sel,:)','r-');
loglog(dis_kb(sel),spf*dTag_wt_500_g1(sel,1),'bo','MarkerSize',4); hold on;
% plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_wt500(sel,:)','b-');
xx = logspace(1,5,5);
plot(xx,5e-5*xx.^(5/3),'--','color',gry1);
xlim([600,1e5]); box off;
ylim([.9,1e3]);
xlabel('genomic separation (kb)');
ylabel('G1 search time to 500 nm')
% QuickLogFit( dis_kb(sel),spf*dTag_wt_500(sel,1) ); box off;

figure(3); clf;
subplot(1,2,1); sel=1:8;
loglog(dis_kb(sel),spf*nt_pt_50_g1(sel,3),'r^','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_pt50_g1(sel,:,3)','r-');
loglog(dis_kb(sel),spf*dTag_pt_50_g1(sel,3),'bo','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_pt50_g1(sel,:,3)','b-');
plot([tb,tb],[1,5e4],'--','color',tc);
xx = logspace(1,5,5);
plot(xx,.02*xx.^(2),'--','color',gry2);
plot(xx,.06*xx.^(5/3),'--','color',gry1);
box off; 
xlim([4,1e5]);  ylim([1,.5e4]);
xlabel('genome distance (kb)');
ylabel('passage time (s)')
title('G1 passage time')

subplot(1,2,2);
sel =8:11;
loglog(dis_kb(sel),spf*nt_pt_500_g1(sel,1),'r^','MarkerSize',4); hold on
plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_pt500_g1(sel,:,1)','r-');
loglog(dis_kb(sel),spf*dTag_pt_500_g1(sel,1),'bo','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_pt500_g1(sel,:,1)','b-');
xx = logspace(1,5,5);
plot(xx,5e-5*xx.^(5/3),'--','color',gry1);
xlim([600,1e5]); box off;
ylim([.1,1e3]);
xlabel('genomic separation (kb)');
ylabel('passage time to 500 nm')
% QuickLogFit( dis_kb(sel),spf*dTag_wt_500(sel,1) ); box off;

figure(1); subplot(2,4,1);
QuickLogFit(dis_kb(2:10),[spf*dTag_wt_50_g1(2:8,3); sc*spf*dTag_wt_500_g1(9:10,1)])
QuickLogFit(dis_kb(2:7),[spf*nt_wt_50_g1(2:7,3)],'plotOpts',{'r-'})



%% 
theta = 50;

[nt_wt_50_g2,   nt_ci_wt50_g2] = CompContactWaitingKM(data3D_g2(1,:),'d_contact',theta);
[dTag_wt_50_g2, dTag_ci_wt50_g2] = CompContactWaitingKM(data3D_g2(2,:),'d_contact',theta);

[nt_pt_50_g2,   nt_ci_pt50_g2] = CompFirstPassageKM(data3D_g2(1,:),'d_contact',theta);
[dTag_pt_50_g2, dTag_ci_pt50_g2] = CompFirstPassageKM(data3D_g2(2,:),'d_contact',theta);

[p3_freq_nt_g2,p3_freq_nt_ci_g2] = CompPassageFreq(data3D_g2(1,:),'d_contact',theta);
[p3_freq_dTag_g2,p3_freq_dTag_ci_g2] = CompPassageFreq(data3D_g2(2,:),'d_contact',theta);

[cd_nt_g2,cd_nt_ci_g2] = ContactDuration(data3D_g2(1,:),'theta',theta);
[cd_dTag_g2,cd_dTag_ci_g2] = ContactDuration(data3D_g2(2,:),'theta',theta);

[nt_wt_500_g2,   nt_ci_wt500_g2] = CompContactWaitingKM(data3D_g2(1,:),'d_contact',500);
[dTag_wt_500_g2, dTag_ci_wt500_g2] = CompContactWaitingKM(data3D_g2(2,:),'d_contact',500);

[nt_pt_500_g2,   nt_ci_pt500_g2] = CompFirstPassageKM(data3D_g2(1,:),'d_contact',500);
[dTag_pt_500_g2, dTag_ci_pt500_g2] = CompFirstPassageKM(data3D_g2(2,:),'d_contact',500);

%%
gry1 = [.4 .4 .4];
gry2 = [.8 .8 .8];
sel = 1:8;
spf = 0.5; % seconds per frame 
tb = 600; % TAD border
tc = [.75 0 0]; % TAD border line color
figure(1);

subplot(2,4,5);
loglog(dis_kb(sel),spf*nt_wt_50_g2(sel,3),'r^','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_wt50_g2(sel,:,3)','r-');
loglog(dis_kb(sel),spf*dTag_wt_50_g2(sel,3),'bo','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_wt50_g2(sel,:,3)','b-');
sc = max([nt_wt_50_g2(8,3)./nt_wt_500_g2(8,1),dTag_wt_50_g2(8,3)./dTag_wt_500_g2(8,1)]);
sel = 8:11; 
loglog(dis_kb(sel),sc*spf*nt_wt_500_g2(sel,1),'r>','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],sc*spf*nt_ci_wt500_g2(sel,:,1)','r-');
loglog(dis_kb(sel),sc*spf*dTag_wt_500_g2(sel,1),'bs','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],sc*spf*dTag_ci_wt500_g2(sel,:,1)','b-');
plot([tb,tb],[1,5e6],'--','color',tc);
xx = logspace(1,5,5);
rectangle('Position',[14,6,1.2e3,.2e4],'edgecolor',.7*ones(1,3));
plot(xx,.064*xx.^(2),'--','color',gry2);
plot(xx,.19*xx.^(5/3),'--','color',gry1);
box off; 
xlim([4,1e5]);  ylim([4,.8e6]);
xlabel('genome distance (kb)');
ylabel('search time (s)');
title('g2 search time')

sel = 1:8;
subplot(2,4,6);
loglog(dis_kb(sel),spf*nt_wt_50_g2(sel,3),'r^','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*nt_ci_wt50_g2(sel,:,3)','r-');
loglog(dis_kb(sel),spf*dTag_wt_50_g2(sel,3),'bo','MarkerSize',4); hold on;
plot([dis_kb(sel); dis_kb(sel)],spf*dTag_ci_wt50_g2(sel,:,3)','b-');
plot([tb,tb],[1,5e6],'--','color',tc);
xx = logspace(1,5,5);
plot(xx,.064*xx.^(2),'--','color',gry2);
plot(xx,.19*xx.^(5/3),'--','color',gry1);
box off; 
xlim([14,1.2e3]);  ylim([6,.2e4]);
xlabel('genome distance (kb)');
ylabel('search time (s)');
title('G2 search time')

subplot(2,4,7);
fph = 2*60*60;  % frames to hours
loglog(dis_kb',fph*p3_freq_nt_g2,'r^','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],fph*p3_freq_nt_ci_g2','r-'); hold on;
plot(dis_kb',fph*p3_freq_dTag_g2,'bo','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],fph*p3_freq_dTag_ci_g2','b-'); hold on;
 plot([tb,tb],[1e-2,1e2],'--','color',tc);
 xlim([4,1e5]);
box off;
xlabel('genomic separation (kb)');
ylabel('passage frequency (events/hr)');
set(gcf,'color','w');
title('G2 passage frequency')
% title(['contact threshold=',num2str(theta)]);

subplot(2,4,8);
loglog(dis_kb',cd_nt_g2,'r^','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],cd_nt_ci_g2,'r-'); hold on;
plot(dis_kb',cd_dTag_g2,'bo','MarkerSize',4); hold on;
plot([dis_kb;dis_kb],cd_dTag_ci_g2,'b-'); hold on;
plot([tb,tb],[.1,10],'--','color',tc);
box off;
xlabel('genomic separation (kb)');
ylabel('contact duration (s)');
ylim([.9 1.4]);
title('G2 contact duration')
% title(['contact threshold=',num2str(theta)]);


figure(1); subplot(2,4,5);
QuickLogFit(dis_kb(2:10),[spf*dTag_wt_50_g2(2:8,3); sc*spf*dTag_wt_500_g2(9:10,1)])
QuickLogFit(dis_kb(2:7),[spf*nt_wt_50_g2(2:7,3)],'plotOpts',{'r-'})



%% domain 2 cell cycle


%% domain 2 computations g1/g2

data3D_2_g1 = data3D_2;
data3D_2_g2 = data3D_2;
for i=1:numel(data3D_2_g1)
    g1 = data3D_2_g1{i}(:,end) < 5;
    data3D_2_g1{i} = data3D_2_g1{i}(g1,:);
    data3D_2_g2{i} = data3D_2_g2{i}(~g1,:);
end

%% compute domain 2, G1

theta = 50;

[nt2_wt_50_g1,   nt2_ci_wt50_g1] = CompContactWaitingKM(data3D_2_g1(1,:),'d_contact',theta);
[dTag2_wt_50_g1, dTag2_ci_wt50_g1] = CompContactWaitingKM(data3D_2_g1(2,:),'d_contact',theta);

[nt2_fp_50_g1,   nt2_ci_fp50_g1] = CompFirstPassageKM(data3D_2_g1(1,:),'d_contact',theta);
[dTag2_fp_50_g1, dTag2_ci_fp50_g1] = CompFirstPassageKM(data3D_2_g1(2,:),'d_contact',theta);

[pt2_freq_nt_g1,pt2_freq_nt_ci_g1] = CompPassageFreq(data3D_2_g1(1,:),'d_contact',theta);
[pt2_freq_dTag_g1,pt2_freq_dTag_ci_g1] = CompPassageFreq(data3D_2_g1(2,:),'d_contact',theta);

[cd2_nt_g1,cd2_nt_ci_g1] = ContactDuration(data3D_2_g1(1,:),'theta',theta);
[cd2_dTag_g1,cd2_dTag_ci_g1] = ContactDuration(data3D_2_g1(2,:),'theta',theta);



%% plot results domain 2, G1
figure(6); clf;
subplot(2,4,1); 
semilogy(dis_kb2_sign,spf*nt2_fp_50_g1(:,3),'r^','MarkerSize',4); hold on;
plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_fp50_g1(:,:,3)','r-'); hold on;
plot([-50,-50],[4,5e3],'k--');
box off; xlim([-500 500]); ylim([4,800]);  ylabel('G1  median passage time (s)')

subplot(2,4,2); 
semilogy(dis_kb2_sign,fph*pt2_freq_nt_g1,'r^','MarkerSize',4); hold on;
plot([dis_kb2_sign;dis_kb2_sign],fph*pt2_freq_nt_ci_g1','r-'); hold on;
plot([-50,-50],[1,20],'k--');
box off; xlim([-500 500]);   ylim([3.,15]);
ylabel('passage freq. (events/hr)');

subplot(2,4,3); 
semilogy(dis_kb2_sign,spf*nt2_wt_50_g1(:,3),'r^','MarkerSize',4); hold on;
plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_wt50_g1(:,:,3)','r-'); hold on;
plot([-50,-50],[10,1e3],'k--');
box off; xlim([-500 500]); ylim([10,1e3]); ylabel('med. search time (s)');

subplot(2,4,4); 
semilogy(dis_kb2_sign,cd2_nt_g1,'r^'); hold on;
plot([dis_kb2_sign;dis_kb2_sign],cd2_nt_ci_g1,'r-'); hold on;
plot([-50,-50],[.5,1.4],'k--');
box off; xlim([-500 500]);
ylim([.9,1.3])
ylabel('ave. contact duration (s)');

%% compute domain 2 G2
theta = 50;

[nt2_wt_50_g2,   nt2_ci_wt50_g2] = CompContactWaitingKM(data3D_2_g2(1,:),'d_contact',theta);
[dTag2_wt_50_g2, dTag2_ci_wt50_g2] = CompContactWaitingKM(data3D_2_g2(2,:),'d_contact',theta);

[nt2_fp_50_g2,   nt2_ci_fp50_g2] = CompFirstPassageKM(data3D_2_g2(1,:),'d_contact',theta);
[dTag2_fp_50_g2, dTag2_ci_fp50_g2] = CompFirstPassageKM(data3D_2_g2(2,:),'d_contact',theta);

[pt2_freq_nt_g2,pt2_freq_nt_ci_g2] = CompPassageFreq(data3D_2_g2(1,:),'d_contact',theta);
[pt2_freq_dTag_g2,pt2_freq_dTag_ci_g2] = CompPassageFreq(data3D_2_g2(2,:),'d_contact',theta);

[cd2_nt_g2,cd2_nt_ci_g2] = ContactDuration(data3D_2_g2(1,:),'theta',theta);
[cd2_dTag_g2,cd2_dTag_ci_g2] = ContactDuration(data3D_2_g2(2,:),'theta',theta);



%% plot results (domain 2, G2)
figure(6); 
subplot(2,4,5); 
semilogy(dis_kb2_sign,spf*nt2_fp_50_g2(:,3),'r^','MarkerSize',4); hold on;
plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_fp50_g2(:,:,3)','r-'); hold on;
plot([-50,-50],[4,1e3],'k--');
box off; xlim([-500 500]); ylim([4,130]);  ylabel('G2  median passage time (s)')

subplot(2,4,6); 
semilogy(dis_kb2_sign,fph*pt2_freq_nt_g2,'r^','MarkerSize',4); hold on;
plot([dis_kb2_sign;dis_kb2_sign],fph*pt2_freq_nt_ci_g2','r-'); hold on;
plot([-50,-50],[1,20],'k--');
box off; xlim([-500 500]);   ylim([3.,15]);
ylabel('passage freq. (events/hr)');

subplot(2,4,7); 
semilogy(dis_kb2_sign,spf*nt2_wt_50_g2(:,3),'r^','MarkerSize',4); hold on;
plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_wt50_g2(:,:,3)','r-'); hold on;
plot([-50,-50],[10,1e3],'k--');
box off; xlim([-500 500]); ylim([10,1e3]); ylabel('med. search time (s)');

subplot(2,4,8); 
semilogy(dis_kb2_sign,cd2_nt_g2,'r^'); hold on;
plot([dis_kb2_sign;dis_kb2_sign],cd2_nt_ci_g2,'r-'); hold on;
plot([-50,-50],[.5,1.4],'k--');
box off; xlim([-500 500]);
ylim([.9,1.3])
ylabel('ave. contact duration (s)');




%% Domain 2, threshold effects
ths = 4; %  length(thetas); % above 350 and 500 nm search times are often 1-2 frames and hence unreliable. 
figure(6); clf;
for t=1:ths   
    theta = thetas(t); % 50;   
    [nt2_wt_50_th,   nt2_ci_wt50_th] = CompContactWaitingKM(data3D_2(1,:),'d_contact',theta);
    [dTag2_wt_50_th, dTag2_ci_wt50_th] = CompContactWaitingKM(data3D_2(2,:),'d_contact',theta);
    [nt2_fp_50_th,   nt2_ci_fp50_th] = CompFirstPassageKM(data3D_2(1,:),'d_contact',theta);
    [dTag2_fp_50_th, dTag2_ci_fp50_th] = CompFirstPassageKM(data3D_2(2,:),'d_contact',theta);
    [pt2_freq_nt_th,pt2_freq_nt_ci_th] = CompPassageFreq(data3D_2(1,:),'d_contact',theta);
    [pt2_freq_dTag_th,pt2_freq_dTag_ci_th] = CompPassageFreq(data3D_2(2,:),'d_contact',theta);
    [cd2_nt_th,cd2_nt_ci_th] = ContactDuration(data3D_2(1,:),'theta',theta);
    [cd2_dTag_th,cd2_dTag_ci_th] = ContactDuration(data3D_2(2,:),'theta',theta);    
    % plot results
    subplot(ths,4,1+(t-1)*4); 
        semilogy(dis_kb2_sign,spf*nt2_fp_50_th(:,3),'r^','MarkerSize',4); hold on;
        plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_fp50_th(:,:,3)','r-'); hold on;
        plot([-50,-50],[4,5e3],'k--');
        box off;  xlim([-500 500]); %  ylim([4,800]);  
        ylabel(['theta = ',num2str(theta) ,'  median passage time (s)'])  
    subplot(ths,4,2+(t-1)*4); 
        semilogy(dis_kb2_sign,fph*pt2_freq_nt_th,'r^','MarkerSize',4); hold on;
        plot([dis_kb2_sign;dis_kb2_sign],fph*pt2_freq_nt_ci_th','r-'); hold on;
        plot([-50,-50],[1,200],'k--');
        box off; xlim([-500 500]);  %  ylim([3.,15]);
        ylabel('passage freq. (events/hr)');    
    subplot(ths,4,3+(t-1)*4); 
        semilogy(dis_kb2_sign,spf*nt2_wt_50_th(:,3),'r^','MarkerSize',4); hold on;
        plot([dis_kb2_sign;dis_kb2_sign],spf*nt2_ci_wt50_th(:,:,3)','r-'); hold on;
        plot([-50,-50],[1,1e3],'k--');
        box off; xlim([-500 500]); % ylim([10,1e3]); 
        ylabel('med. search time (s)');
    subplot(ths,4,4+(t-1)*4); 
        semilogy(dis_kb2_sign,cd2_nt_th,'r^'); hold on;
        plot([dis_kb2_sign;dis_kb2_sign],cd2_nt_ci_th,'r-'); hold on;
        plot([-50,-50],[.5,10],'k--');
        box off; xlim([-500 500]);
       %  ylim([.9,1.3])
        ylabel('ave. contact duration (s)');
end

