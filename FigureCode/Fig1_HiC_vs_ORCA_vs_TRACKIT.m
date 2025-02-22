
tic
dataFolder = 'U:\Manuscripts\Jude Live Imaging\Data\traces_v13\';
dis_kb = [5,20,55,70,134,260,407,799,2030,12795,73567];
nE = length(dis_kb);

data3D = cell(11,2);
for e=1:nE    
    datafile = [dataFolder,'Dis3D_2Hz_noTreatment_',num2str(dis_kb(e)),'kb.csv'];
    data3D{1,e} = readmatrix(datafile);
    datafile = [dataFolder,'Dis3D_2Hz_dTag_',num2str(dis_kb(e)),'kb.csv'];
    data3D{2,e} = readmatrix(datafile);
end
toc

% genomic positions of the insertions
tetO = 51320704;
cuO= [51321893
51336653
51371176
51387332
51451617
51576421
51724680
50525643
49294139
38529185
124883209];




%% Load Hi-C data
% Here we also convert between the ORCA-barcodes and the insertion
% positions used for live imaging

if 1 % ~exist('hox5kb','var')

jboxFolder = 'U:\GenomeData\JuiceboxExport\';
jb = 'U:\GenomeData\JuiceboxExport\';

% Bonev_mESC_hoxA_6Mb_balanced_5kb.txt
% Bonev_mESC_hoxA_30Mb_balanced_25kb.txt
% Bonev_mESC_hoxA_120Mb_balanced_100kb.txt

hox5kb = ReadJuiceboxMatrix([jb,'Bonev_mESC_hoxA_6Mb_balanced_5kb.txt'],'locus','6:51010704-52500591');
p1_5 = round((tetO - 51010704)/5e3);
p2_5 = round( (cuO(1:7) - 51010704)/5e3);
p1_kb = 51010704 + 5e3*repmat(p1_5,length(p2_5),1);
p2_kb = 51010704 + 5e3*p2_5;
hic_5 = hox5kb(p1_5,p2_5)';

hox25kb = ReadJuiceboxMatrix([jb,'Bonev_mESC_hoxA_30Mb_balanced_25kb.txt'],'locus','6:36000000-53000591');
p1_25 = round((tetO - 36000000)/25e3);
p2_25 = round( (cuO(8:10) - 36000000)/25e3);
hic_25 = hox25kb(p1_25,p2_25)'/( (25/5)^2  );
p1_kb_25 = 36000000 + 25e3*repmat(p1_25,length(p2_25),1);
p2_kb_25 = 36000000 + 25e3*p2_25;

hox100kb = ReadJuiceboxMatrix([jb,'Bonev_mESC_hoxA_120Mb_balanced_100kb.txt'],'locus','6:20,000,000-130,000,000');
p1_100 = round((tetO - 20e6)/100e3);
p2_100 = round( (cuO(11) - 20e6)/100e3);
hic_100 = hox100kb(p1_100,p2_100)'/( (100/5)^2  );
p1_kb_100 = 20e6 + 100e3*repmat(p1_100,length(p2_100),1);
p2_kb_100 = 20e6 + 100e3*p2_100;

ps = [p1_kb, p2_kb;
      p1_kb_25,p2_kb_25;
      p1_kb_100,p2_kb_100];
hic =[ hic_5; hic_25; hic_100];

figure(1); clf; loglog(abs(ps(:,2)-ps(:,1)),hic,'.-');
xlabel('separation (bp)'); ylabel('contact frequency')
set(gcf,'color','w'); 

% view Hi-C data
figure(10); clf;
subplot(1,3,1);
gx_5 = 51010704:5e3:52500591;
imagesc(gx_5,gx_5,log2(hox5kb)); colorbar; GetColorMap('WORBK');
hold on; plot(p1_kb,p2_kb,'bo');  clim([0,7]);
title('1.5 Mb, 5 kb res')

subplot(1,3,2);
gx =36000000:25e3:53000591;
imagesc(gx,gx,log2(hox25kb)); colorbar; GetColorMap('WORBK');
hold on; plot(ps(1:10,1),ps(1:10,2),'bo'); clim([0,10])
rectangle('Position',[51010704,51010704,52500591-51010704,52500591-51010704],'EdgeColor',[0,0,.5]);
title('17 Mb, 25 kb res')

subplot(1,3,3);
gx =  20e6:100e3:130e6;
imagesc(gx,gx,log2(hox100kb)); colorbar; GetColorMap('WORBK');
hold on; plot(ps(:,1),ps(:,2),'bo'); clim([3,10]);
rectangle('Position',[36000000,36000000,53000591-36000000,53000591-36000000],'EdgeColor',[0,0,.5]);
title('110 Mb, 100 kb res')
set(gcf,'color','w');

end
%% overlay the hic and live data

cmap = hsv(nE+1);
thetas = 100; % [30,50,100,150,200,250,350];
for t=1:length(thetas)
    theta = thetas(t);
    nD = length(dis_kb);
    med_3D_nt = zeros(nD,1);
    med_3D_dTag = zeros(nD,1);
    
    mean_3D_nt = zeros(nD,1);
    mean_3D_dTag = zeros(nD,1);
    contact_3D_nt = nan(nD,1);
    contact_3D_nt_ci = nan(nD,2);
    contact_3D_dTag = nan(nD,1);
    for d=1:nD
        disNT = data3D{1,d};
        nObs = sum(~isnan(disNT),2);
        [~,idx] = sort(nObs,'descend');
        disNT = disNT(idx,:);
        med_3D_nt(d) = nanmedian(disNT(:));
        contact_3D_nt(d) =  sum(disNT(:)<theta)./sum(disNT(:)<inf);
        disNTlin = disNT(:);
   
        % % supress for speed (error bars are invisbly small anyway)
        disNTlin(isnan(disNTlin)) = []; % seeds up bootstrapping
        boot = randi(length(disNTlin),length(disNTlin),10);
        contact_3D_nt_ci(d,1) =  quantile(sum(disNTlin(boot)<theta)./sum(disNTlin(boot)<inf),.05,2);
        contact_3D_nt_ci(d,2) =  quantile(sum(disNTlin(boot)<theta)./sum(disNTlin(boot)<inf),.95,2);
    end
    dg = abs(ps(:,2)-ps(:,1)); dg(1) = 5e3; % genome dist;  correct 0 to 5kb 
    
    figure(20+t); clf;  loglog(hic,contact_3D_nt,'.','MarkerSize',1); box off; hold on;
    plot([hic,hic]',contact_3D_nt_ci','k-'); box off; hold on;
    r = CorCoef(log10(hic),log10(contact_3D_nt))
    for d=1:11
        plot(hic(d),contact_3D_nt(d),'.','MarkerSize',15,'color',cmap(d,:));
    end
    title(['Pearsons R= ',num2str(r,3)]);
    xlabel('contact freq. Hi-C (Bonev)');
    ylabel('contact freq.  TRACK-IT')
end

%%