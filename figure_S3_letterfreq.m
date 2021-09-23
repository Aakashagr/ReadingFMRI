%% Figure 2: Average RT and dissimilarity correlation plot.
allclear
load L2_letters; load tel_bigram; load mal_bigram; load tel_lf; load mal_lf;

% removing outliers
L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.
mat = zeros(size(L2_str.RT)); mat(L2_str.RT>5) = 1; % identifying harder search pairs
temp = sum(sum(mat,2),3);
% figure; stem(temp) % identifying mistakes.
id = find(temp > 0 & temp < 7);

% removing higher RT outlier
for i = 1:length(id)
    [c, v] =ind2sub([size(L2_str.RT,2),2],find(mat(id(i),:,:)==1));
    L2_str.RT(id(i),c,v) = NaN;
end

% % subjects to remove
L2_str.RT(:,[29 28 12],:) = [];
L2_str.subjinfo.ismalayalam([29 28 12]) = [];
L2_str.fluencytest.reading_time([29 28 12]) = [];

images = L2_str.images;
image_pairs = L2_str.img_pairs(1:1260,:);
RT = L2_str.RT;
ismal = L2_str.subjinfo.ismalayalam;
mal = 1:630; tel = 631:1260;

% Search time
RTmm = nanmean(nanmean(RT(mal,ismal == 1,:),2),3);
RTmt = nanmean(nanmean(RT(tel,ismal == 1,:),2),3);
RTtt = nanmean(nanmean(RT(tel,ismal == 0,:),2),3);
RTtm = nanmean(nanmean(RT(mal,ismal == 0,:),2),3);

tel_lf = tel_lf./sum(tel_lf);  mal_lf = mal_lf./sum(mal_lf); 
tbf = tel_lf(nchoosek(1:36,2));
mbf = mal_lf(nchoosek(1:36,2));

%%
TB = nanmean(tel_bigram,2); TB = TB/nansum(TB)*100;
MB = nanmean(mal_bigram,2); MB = MB/nansum(MB)*100;
TB(TB >= 1) = NaN; MB(MB >= 1) = NaN;
figure(1); corrplot(MB,(1./RTmm - 1./RTtm)); xlabel('% Bigram frequency'); ylabel('d(Malayalam) - d(Telugu) readers');hold on;
figure(2); corrplot(TB,(1./RTtt - 1./RTmt)); xlabel('% Bigram frequency'); ylabel('d(Telugu) - d(Malayalam) readers');hold on;
% [r p] = partialcorri((1./RTtt - 1./RTmt),[TB,tbf]); 

TB = nanmean(tel_bigram,2); TB = TB/nansum(TB)*100;
MB = nanmean(mal_bigram,2); MB = MB/nansum(MB)*100;
TB(TB < 1) = NaN; MB(MB < 1) = NaN;
figure(1); corrplot(MB,(1./RTmm - 1./RTtm)); xlabel('% Bigram frequency'); ylabel('d(Malayalam) - d(Telugu) readers'); lsline
figure(2); corrplot(TB,(1./RTtt - 1./RTmt)); xlabel('% Bigram frequency'); ylabel('d(Telugu) - d(Malayalam) readers'); lsline




