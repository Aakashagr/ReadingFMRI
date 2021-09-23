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

images = L2_str.images;
image_pairs = L2_str.img_pairs;
RT = L2_str.RT;
part_id = L2_str.partsTD;

ismal = L2_str.subjinfo.ismalayalam;
RT_mal = RT(:,ismal == 1,:); 
RT_tel = RT(:,ismal == 0,:);
PC_mal = L2_str.PC(:,ismal == 1,:);
PC_tel = L2_str.PC(:,ismal == 0,:);

mRT_mal = nanmean(RT_mal,3);
mRT_tel = nanmean(RT_tel,3);

% Malayalam and telugu pairs
tel = 1:300; mal = 301:600;

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);
RTmm = RT_mal_mean(mal);
RTmt = RT_mal_mean(tel);
RTtt = RT_tel_mean(tel);
RTtm = RT_tel_mean(mal);

%% Baton model

[y_pred_tt,w_tt,~,mat] = Batonmodel(1./RTtt, 5, 2,1,0);
[y_pred_mt,w_mt] = Batonmodel(1./RTmt, 5, 2,1,0);
[y_pred_mm,w_mm] = Batonmodel(1./RTmm, 5, 2,1,0);
[y_pred_tm,w_tm] = Batonmodel(1./RTtm, 5, 2,1,0);

% sigmoid fit
[~,y_pred_new_tt] = fitsigmoid(y_pred_tt,1./RTtt);
[~,y_pred_new_mt] = fitsigmoid(y_pred_mt,1./RTmt);
[~,y_pred_new_tm] = fitsigmoid(y_pred_tm,1./RTtm);
[~,y_pred_new_mm] = fitsigmoid(y_pred_mm,1./RTmm);

%%
% Calculating errors
ett = 1./RTtt - y_pred_new_tt; etm = 1./RTtm - y_pred_new_tm;
emm = 1./RTmm - y_pred_new_mm; emt = 1./RTmt - y_pred_new_mt;

% Loading bigram frequency data
load('tel_25bigram.mat')
load('mal_25bigram.mat')

% fitting bigram frequency model to the residual data
pairs = nchoosek(1:25,2);
Xt = [tel_bigram(pairs(:,1)) tel_bigram(pairs(:,2)) ones(size(ett,1),1)];
Xm = [mal_bigram(pairs(:,1)) mal_bigram(pairs(:,2)) ones(size(ett,1),1)];

% Telugu stimuli
w = regress(ett,Xt); subplot(221); corrplot(Xt*w,ett);
w = regress(emt,Xt); subplot(222); corrplot(Xt*w,emt);

% Malayalam stimuli
w = regress(emm,Xm); subplot(223); corrplot(Xm*w,emm);
w = regress(etm,Xm); subplot(224); corrplot(Xm*w,etm);

