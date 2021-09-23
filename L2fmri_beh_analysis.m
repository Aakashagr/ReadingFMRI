allclear
load('L2fmri_bigram.mat')

L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.

images = L2_str.images;
image_pairs = L2_str.img_pairs;

RT_mal = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:);
RT_tel = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:);

mRT_mal = reshape(RT_mal,[size(RT_mal,1) size(RT_mal,2)*size(RT_mal,3)]);
mRT_tel = reshape(RT_tel,[size(RT_tel,1) size(RT_tel,2)*size(RT_tel,3)]);

% Malayalam and telugu pairs
mal = 1:276; tel = 277:552;  

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);
RTmm = RT_mal_mean(mal);
RTmt = RT_mal_mean(tel);
RTtt = RT_tel_mean(tel);
RTtm = RT_tel_mean(mal);

splithalfcorr(mean(RT_tel,3)')
splithalfcorr(mean(RT_mal,3)')

%% Correlating with predicted data
load pred_2RT_fmri
figure; 
% 1./RT plots 
subplot(221); corrplot(1./RTtt, pred_tt,'Tel sub on Tel pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')
subplot(224); corrplot(1./RTtm, pred_tm,'Tel sub on Mal pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')
subplot(223); corrplot(1./RTmm, pred_mm,'Mal sub on Mal pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')
subplot(222); corrplot(1./RTmt, pred_mt,'Mal sub on Tel pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')


%% correlation across subjects
[R, P] = corr(nanmean(L2_str.RT,3)); R(P==1) = NaN;
figure;imagesc(R); colorbar;  caxis([min(R(:))-.2 max(R(:))])
set(gca,'Xtick',1:16); set(gca,'Ytick',1:16); 
for i = 1:7; l{i} = ['m',num2str(i)]; end
for i = 1:8; l{end+1} = ['t',num2str(i)]; end
set(gca,'Yticklabel',l);
set(gca,'Xticklabel',l);
title('Reaction time correlation across subjects')
% Corrected Split half correlation
splithalfcorr(mean(RT_tel,3)')

%% Summary plot

figure; barweb([mean(RTtt), mean(RTmt); mean(RTtm), mean(RTmm)], ...
    [nansem(RTtt), nansem(RTmt); nansem(RTtm), nansem(RTmm)]);
ylabel('Mean reaction time, seconds');   colormap(get(gca,'ColorOrder'));
set(gca,'Xticklabel',{'Telugu Bigrams','Malayalam Bigrams'});
legend({'Telugu readers','Malayalam readers'});

figure; subplot(221); corrplot(RTmt, RTtt,'Telugu-pairs',1); 
xlabel('Malayalam group RT, seconds'); ylabel('Telugu group RT, seconds'); 
subplot(222); corrplot(RTtm, RTmm,'Malayalam-pairs',1);
xlabel('Telugu group RT, seconds'); ylabel('Malayalam group RT, seconds');

% 1./RT plots 
subplot(223); corrplot(1./RTmt, 1./RTtt,'Telugu-pairs',1); 
ylabel('Telugu group dissimilarity,1/s'); xlabel('Malayalam group dissimilarity,1/s');
subplot(224); corrplot(1./RTtm, 1./RTmm,'Malayalam-pairs',1); 
xlabel('Telugu group dissimilarity,1/s'); ylabel('Malayalam group dissimilarity,1/s');

%% Comparing subject groups

% Identifying pairs that are significantly different between two groups
h = zeros(552,1);
for i = 1:552 ;  p(i) = ranksum(mRT_mal(i,:),mRT_tel(i,:)); end
h(p<.05) = 1;
sig_mal = find(h(mal)==1); %sig_mal = sig_mal + min(mal)-1; sig_mal = sig_mal - 300;
sig_tel = find(h(tel)==1); 

% RT plots
figure; subplot(221); corrplot(RTmt, RTtt,'Telugu-pairs',1); hold on
plot(RTmt(sig_tel), RTtt(sig_tel),'r.'); %lsline
ylabel('Telugu group RT, seconds'); xlabel('Malayalam group RT, seconds')

subplot(222); corrplot(RTtm, RTmm,'Malayalam-pairs',1);
plot(RTtm(sig_mal), RTmm(sig_mal),'r.'); %lsline
xlabel('Telugu group RT, seconds'); ylabel('Malayalam group RT, seconds')

% 1./RT plots 
subplot(223); corrplot(1./RTmt, 1./RTtt,'Telugu-pairs',1); xx = find(1./RTmt(sig_tel) < 1./RTtt(sig_tel));
plot(1./RTmt(sig_tel), 1./RTtt(sig_tel),'r.'); 
plot(1./RTmt(sig_tel(xx)),polyval(polyfit(1./RTmt(sig_tel(xx)), 1./RTtt(sig_tel(xx)),1),1./RTmt(sig_tel(xx))),'r')
ylabel('Telugu group dissimilarity,1/s'); xlabel('Malayalam group dissimilarity,1/s')

subplot(224); corrplot(1./RTtm, 1./RTmm,'Malayalam-pairs',1); xx = find(1./RTtm(sig_mal) < 1./RTmm(sig_mal));
plot(1./RTtm(sig_mal), 1./RTmm(sig_mal),'r.'); 
plot(1./RTtm(sig_mal(xx)),polyval(polyfit(1./RTtm(sig_mal(xx)), 1./RTmm(sig_mal(xx)),1),1./RTtm(sig_mal(xx))),'r')
xlabel('Telugu group dissimilarity,1/s'); ylabel('Malayalam group dissimilarity,1/s')


% slope of significant pairs
idm = find(1./RTmm(sig_mal) < 1./RTtm(sig_mal)); sig_mal(idm) =  [];
[w_m, wint_m] = regress(1./RTmm(sig_mal),[1./RTtm(sig_mal),ones(numel(sig_mal),1)])

idt = find(1./RTtt(sig_tel) < 1./RTmt(sig_tel)  | 1./RTtt(sig_tel) < .5); sig_tel(idt) =  [];
[w_t, wint_t] = regress(1./RTtt(sig_tel),[1./RTmt(sig_tel),ones(numel(sig_tel),1)])

%% MDS plots

% Telugu letters
figure;
dtt = mdscale(squareform(1./RTtt),2);
dmt = mdscale(squareform(1./RTmt),2);

[~,Dtt] = procrustes(dmt,dtt,'scaling',0);

h1 = img_scatterplot(dmt(:,1),dmt(:,2),images(1:24),.05); hold all
h = plot(Dtt(:,1),Dtt(:,2),'ks');
for i = 1:size(Dtt,1);  plot([dmt(i,1) Dtt(i,1)],[dmt(i,2) Dtt(i,2)],'k:');  end
legend(h,'Telugu readers')
title('Change in representation w.r.t. Malayalam readers')

% Malayalam letters
figure;
dmm = mdscale(squareform(1./RTmm),2);
dtm = mdscale(squareform(1./RTtm),2);

[~,Dmm] = procrustes(dtm,dmm,'scaling',0);

h1 = img_scatterplot(dtm(:,1),dtm(:,2),images(25:48),.05); hold all
h = plot(Dmm(:,1),Dmm(:,2),'ks');
for i = 1:size(Dmm,1);  plot([dtm(i,1) Dmm(i,1)],[dtm(i,2) Dmm(i,2)],'k:');  end
legend(h,'Malayalam readers')
title('Change in representation w.r.t. Telugu readers')



