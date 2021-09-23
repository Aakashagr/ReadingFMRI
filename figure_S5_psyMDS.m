allclear
load('L2_telmal_bigram.mat')

L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.
mat = zeros(size(L2_str.RT)); mat(L2_str.RT>5) = 1; % identifying harder search pairs
temp = sum(sum(mat,2),3);
% figure; stem(temp) % identifying mistakes.
id = find(temp > 0 & temp < 2);

% removing higher RT outlier
for i = 1:length(id)
    [c, v] =ind2sub([size(L2_str.RT,2),2],find(mat(id(i),:,:)==1));
    L2_str.RT(id(i),c,v) = NaN;
end

for i = 1:numel(L2_str.images); images{i} = imcomplement(L2_str.images{i}); end

RT_mal = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:);
RT_tel = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:);

mRT_mal = nanmean(RT_mal,3);
mRT_tel = nanmean(RT_tel,3);

% Malayalam and telugu pairs
tel = 1:300; mal = 301:600;

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);
RTmm = RT_mal_mean(mal);
RTtt = RT_tel_mean(tel);

%% Telugu bigrams
figure; 
d_mal = mdscale(squareform(1./RT_mal_mean(tel)),2);
d_tel = mdscale(squareform(1./RT_tel_mean(tel)),2);

[~,D_tel] = procrustes(d_mal,d_tel,'scaling',0);

h1 = img_scatterplot(d_mal(:,1),d_mal(:,2),images(1:25),.07); hold all
h = plot(D_tel(:,1),D_tel(:,2),'r+'); 
for i = 1:size(D_tel,1);  plot([d_mal(i,1) D_tel(i,1)],[d_mal(i,2) D_tel(i,2)],'r:');  end
legend(h, {'Telugu readers'})
title('Perceptual space (Telugu Bigrams)')


%% Malayalam bigrams
figure;
d_mal = mdscale(squareform(1./RT_mal_mean(mal)),2);
d_tel = mdscale(squareform(1./RT_tel_mean(mal)),2);

[~,D_mal] = procrustes(d_tel,d_mal,'scaling',0);

h1 = img_scatterplot(d_tel(:,1),d_tel(:,2),images(26:50),.07); hold all
h = plot(D_mal(:,1),D_mal(:,2),'r+');
for i = 1:size(D_mal,1);  plot([d_tel(i,1) D_mal(i,1)],[d_tel(i,2) D_mal(i,2)],'r:');  end
legend(h, {'Malayalam readers'})
title('Perceptual space (Malayalam Bigrams)')
