% MDS plots across all subjects
allclear; 
load L2_letters.mat

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

normalize = 0; % normalize the data by subtracting the mean base line reaction 
if normalize
    for i = 1:numel(L2_str.BRT)
%         brt_mean(i) = mean(L2_str.BRT{i}(4:10));
%         L2_str.RT(:,i,:) =  L2_str.RT(:,i,:) - brt_mean(i);
        prac_mean(i) = mean(L2_str.Practise{i});
        L2_str.RT(:,i,:) =  L2_str.RT(:,i,:) - prac_mean(i);
    end    
end

% % subjects to remove
L2_str.RT(:,[29 28 12],:) = [];
L2_str.subjinfo.ismalayalam([29 28 12]) = [];
 

% extracting variables of interest from L2 structure
images = L2_str.images;
image_pairs = L2_str.img_pairs(1:1260,:);
RT_mal = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:);
RT_tel = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:);
PC_mal = L2_str.PC(:,L2_str.subjinfo.ismalayalam == 1,:);
PC_tel = L2_str.PC(:,L2_str.subjinfo.ismalayalam == 0,:);

% mRT_both = mean(L2_str.RT(:,L2_str.subjinfo.ismalayalam == 2,:),3);
mRT_mal = reshape(RT_mal,[size(RT_mal,1) size(RT_mal,2)*size(RT_mal,3)]);
mRT_tel = reshape(RT_tel,[size(RT_tel,1) size(RT_tel,2)*size(RT_tel,3)]);

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);

% Malayalam and telugu pairs
mal = 1:630;
tel = 631:1260;

%%
% Telugu letters
figure; 
d_mal = mdscale(squareform(1./RT_mal_mean(tel)),2);
d_tel = mdscale(squareform(1./RT_tel_mean(tel)),2);

[~,D_tel] = procrustes(d_mal,d_tel,'scaling',0);

h1 = img_scatterplot(d_mal(:,1),d_mal(:,2),images(37:72),.07); hold all
h = plot(D_tel(:,1),D_tel(:,2),'r+'); 
for i = 1:size(D_tel,1);  plot([d_mal(i,1) D_tel(i,1)],[d_mal(i,2) D_tel(i,2)],'r:');  end
legend(h, {'Telugu readers'})
title('Perceptual space (Telugu letters)')


%% Malayalam letters
figure;

d_mal = mdscale(squareform(1./RT_mal_mean(mal)),2);
d_tel = mdscale(squareform(1./RT_tel_mean(mal)),2);

[~,D_mal] = procrustes(d_tel,d_mal,'scaling',0);

h1 = img_scatterplot(d_tel(:,1),d_tel(:,2),images(1:36),.07); hold all
h = plot(D_mal(:,1),D_mal(:,2),'r+');
for i = 1:size(D_mal,1);  plot([d_tel(i,1) D_mal(i,1)],[d_tel(i,2) D_mal(i,2)],'r:');  end
legend(h, {'Malayalam readers'})
title('Perceptual space (Malayalam letters)')
