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


% First letter different
q = [1:5:25];
pairidx = q(nchoosek(1:5,2));
pairidx = [pairidx;pairidx+1;pairidx+2;pairidx+3;pairidx+4];
qfirst = find(ismember(nchoosek(1:25,2),pairidx,'rows'));

% second letter different
pairidx = nchoosek(1:5,2);
pairidx = [pairidx;pairidx+5;pairidx+10;pairidx+15;pairidx+20];
qsecond = find(ismember(nchoosek(1:25,2),pairidx,'rows'));

%% first letter advantage raw data
figure;subplot(121); barweb([mean(1./RTtt([qfirst qsecond]));mean(1./RTmt([qfirst qsecond]))],[nansem(1./RTtt([qfirst qsecond])); nansem(1./RTmt([qfirst qsecond]))]);
legend('First letter different','Second letter different');  title('Telugu bigram');
 ylabel('mean dissimilarity'); ylim([0 0.6])
set(gca,'XtickLabel',{'Telugu readers','Malayalam readers'});

subplot(122);  barweb([mean(1./RTmm([qfirst qsecond]));mean(1./RTtm([qfirst qsecond]))],[nansem(1./RTmm([qfirst qsecond])); nansem(1./RTtm([qfirst qsecond]))]);
legend('First letter different','Second letter different'); title('Malayalam bigram');
ylabel('mean dissimilarity');ylim([0 0.6])
set(gca,'XtickLabel',{'Telugu readers','Malayalam readers'});

% STATS
ranksum(1./RTmm(qfirst),1./RTmm(qsecond))

%% Baton model

[y_pred_tt,w_tt,~,mat] = Batonmodel(1./RTtt, 5, 2,1,1);
[y_pred_mt,w_mt] = Batonmodel(1./RTmt, 5, 2,1,1);
[y_pred_mm,w_mm] = Batonmodel(1./RTmm, 5, 2,1,1);
[y_pred_tm,w_tm] = Batonmodel(1./RTtm, 5, 2,1,1);

%% Panel E and F
xx = reshape(w_tt(1:40),[10,4]);
figure;barweb([mean(reshape(w_tt(1:40),[10,4]));mean(reshape(w_mt(1:40),[10,4]))]',[nansem(reshape(w_tt(1:40),[10,4])); nansem(reshape(w_mt(1:40),[10,4]))]');
legend('Readers','Nonreaders'); title('Telugu bigram');
xlabel('Parameters'); ylabel('Magnitude');ylim([-.25 .25])
set(gca,'XtickLabel',{'Corresponding 1','Corresponding 2', 'Across part terms','Within part terms'},'XTickLabelRotation',45);

figure; barweb([mean(reshape(w_mm(1:40),[10,4]));mean(reshape(w_tm(1:40),[10,4]))]',[nansem(reshape(w_mm(1:40),[10,4])); nansem(reshape(w_tm(1:40),[10,4]))]');
legend('Readers','Nonreaders'); title('Malayalam bigram');
xlabel('Parameters'); ylabel('Magnitude'); ylim([-.25 .25])
set(gca,'XtickLabel',{'Corresponding 1','Corresponding 2','Across part terms','Within part terms'},'XTickLabelRotation',45);

% STATS
xt = reshape(w_tt(1:40),[10,4]); 
xm = reshape(w_mt(1:40),[10,4]);
yt = reshape(w_tm(1:40),[10,4]);
ym = reshape(w_mm(1:40),[10,4]);
p = statcomparemean(xt(:,1), xm(:,1));

%%
clear
load L2_telmal_3baton.mat

L2_str.RT(:,[2 25],:) = []; L2_str.PC(:,[2 25],:) = [];
L2_str.subjinfo.ismalayalam([2 25]) = []; L2_str.BRT([2 25]) = [];

tr = [2 6:8 10 13 14 16 17 20 22:24]; % Subjects to remove
L2_str.RT(:,tr,:) = []; L2_str.PC(:,tr,:) = [];
L2_str.subjinfo.ismalayalam(tr) = []; L2_str.BRT(tr) = [];

ismal = L2_str.subjinfo.ismalayalam;
    
L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.
mat = zeros(size(L2_str.RT)); mat(L2_str.RT>5) = 1; % identifying harder search pairs
temp = sum(sum(mat,2),3);
% figure; stem(temp) % identifying mistakes.
id = find(temp > 0 & temp <7);

% removing higher RT outlier
for i = 1:length(id)
    [c, v] =ind2sub([size(L2_str.RT,2),2],find(mat(id(i),:,:)==1));
    L2_str.RT(id(i),c,v) = NaN;
end

images = L2_str.images;
mal = 1:325; tel = 326:650;

m_image_pairs = L2_str.img_pairs(mal,:);
t_image_pairs = L2_str.img_pairs(tel,:);
nimg = 6; img_len = 3;
obj_parts = permn(1:nimg,img_len);

RT = L2_str.RT; PC = L2_str.PC;

RT_mal = RT(:,L2_str.subjinfo.ismalayalam == 1,:);  
RT_tel = RT(:,L2_str.subjinfo.ismalayalam == 0,:);
mRT_mal = nanmean(RT_mal,3);
mRT_tel = nanmean(RT_tel,3);

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);
RTmm = RT_mal_mean(mal);
RTmt = RT_mal_mean(tel);
RTtt = RT_tel_mean(tel);
RTtm = RT_tel_mean(mal);
%% Baton model

[y_pred_tt,w_tt] = Batonmodel(1./RTtt, nimg, img_len,1,1,[],t_image_pairs-216);
[y_pred_mt,w_mt] = Batonmodel(1./RTmt, nimg, img_len,1,1,[],t_image_pairs-216);
[y_pred_mm,w_mm] = Batonmodel(1./RTmm, nimg, img_len,1,1,[],m_image_pairs);
[y_pred_tm,w_tm] = Batonmodel(1./RTtm, nimg, img_len,1,1,[],m_image_pairs);

xx = reshape(w_tt(1:135),[15,9]);
figure;barweb([mean(reshape(w_tt(1:135),[15,9]));mean(reshape(w_mt(1:135),[15,9]))]',[nansem(reshape(w_tt(1:135),[15,9])); nansem(reshape(w_mt(1:135),[15,9]))]');
legend('Readers','Nonreaders'); title('Telugu Trigram');
xlabel('Parameters'); ylabel('Magnitude');
set(gca,'XtickLabel',{'C1','C2','C3','An1','An2','Af','Wn1','Wn2','Wf'},'XTickLabelRotation',45);

figure; barweb([mean(reshape(w_mm(1:135),[15,9]));mean(reshape(w_tm(1:135),[15,9]))]',[nansem(reshape(w_mm(1:135),[15,9])); nansem(reshape(w_tm(1:135),[15,9]))]');
legend('Readers','Nonreaders'); title('Malayalam Trigram');
xlabel('Parameters'); ylabel('Magnitude');
set(gca,'XtickLabel',{'C1','C2','C3','An1','An2','Af','Wn1','Wn2','Wf'},'XTickLabelRotation',45);

