%% Using reduced model

allclear
% Upright inverted stimuli
load('L2_telmal_bigram.mat')
L2_str.RT = rmRToutlier(L2_str.RT,4,3);

RT_mal = nanmean(L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:),3);
RT_tel = nanmean(L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:),3);
tel = 1:300; mal = 301:600;

RTmm = RT_mal(mal,:); RTmt = RT_mal(tel,:);
RTtt = RT_tel(tel,:); RTtm = RT_tel(mal,:);

% Fitting baton model
load('pair_RT_2baton.mat');
for sub = 1:8
    [~,bmm(sub,:)] = Batonmodel(1./RTmm(:,sub), 5, 2,3,0,[],[],1./rt_mm);
   
    % for telugu subjects
    [~,btt(sub,:)] = Batonmodel(1./RTtt(:,sub), 5, 2,3,0,[],[],1./rt_tt);
end
fluency = [51.19 67.59 133.42 77 66.35 166.77 67.1 298.8 71.36 119.52 108.3 148.34 156.03 116.8 143.24 68.44]; % first 8 are malayalam
baseline1 = cellfun(@(x) nanmean(x(5:10)), L2_str.BRT);
%% 
load('L2_EMTbigram_searches.mat')
qt_sub = L2_str.subjinfo.ismalayalam==0;
qm_sub = L2_str.subjinfo.ismalayalam==1;
qe_sub = L2_str.subjinfo.ismalayalam==2;

RTt = rmRToutlier(L2_str.RT(:,qt_sub,:),2,2);
RTm = rmRToutlier(L2_str.RT(:,qm_sub,:),2,2);
RTe = rmRToutlier(L2_str.RT(:,qe_sub,:),3,2);

RT_mal = nanmean(RTm,3);
RT_tel = nanmean(RTt,3);
RT_eng = nanmean(RTe,3);

for sub = 1:size(RT_mal,2); [~,bm(sub,:)] = Batonmodel(1./RT_mal(:,sub), 5, 2,3,0,[],[],1./rt_mm); end
fluency_mal = L2_str.subjinfo.reading_time(qm_sub);

% Telugu
for sub = 1:size(RT_tel,2); [~,bt(sub,:)] = Batonmodel(1./RT_tel(:,sub), 5, 2,3,0,[],[],1./rt_tt); end
fluency_tel = L2_str.subjinfo.reading_time(qt_sub);

% English
load dis_letter; q1 = [1 12 18 19 20];
qidx = find(ismember(nchoosek(1:26,2),q1(nchoosek(1:5,2)),'rows'));
for sub = 1:size(RT_eng,2); [~,be(sub,:)] = Batonmodel(1./RT_eng(:,sub), 5, 2,3,0,[],[],dis_uletter(qidx)); end
fluency_eng = L2_str.subjinfo.reading_time(qe_sub);

baseline2t = cellfun(@(x) nanmean(x(5:10)), L2_str.BRT(qt_sub));
baseline2m = cellfun(@(x) nanmean(x(5:10)), L2_str.BRT(qm_sub));

%% fMRI bigrams
load('L2fmri_bigram.mat')
L2_str.RT = rmRToutlier(L2_str.RT,4,5);
RT_mal = nanmean(L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:),3);
RT_tel = nanmean(L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:),3);
mal = 1:276; tel = 277:552;  

RTmm = RT_mal(mal,:); RTmt = RT_mal(tel,:);
RTtt = RT_tel(tel,:); RTtm = RT_tel(mal,:);
load('pair_RT_fmri.mat')
rt_mm = nanmean(nanmean(rt_mm,3),2); rt_mt = nanmean(nanmean(rt_mt,3),2); 
rt_tt = nanmean(nanmean(rt_tt,3),2); rt_tm = nanmean(nanmean(rt_tm,3),2);

img_pair = sortrows([2 4; 6 2; 6 4; 8 2; 8 4; 8 6; 10 2; 10 6; 10 8; 5 2; 1 2; 1 8;...
            4 10; 7 6; 5 6; 10 7; 3 9; 9 5; 9 1; 9 7; 3 1; 3 5; 5 1; 3 7]);

for sub = 1:7
    [~,Fbmm(sub,:)] = Batonmodel(1./RTmm(:,sub), 10, 2,3,0,img_pair,[],1./rt_mm);
end
for sub = 1:4
    % for telugu subjects
    [~,Fbtt(sub,:)] = Batonmodel(1./RTtt(:,sub), 10, 2,3,0,img_pair,[],1./rt_tt);
end
fluency_fmri = [142.2 61.01 119.07 44.82 35.5 30.89 195.36 47.28 144.94 75.24 58.21]; % first 7 are malayalam
baseline3 = cellfun(@(x) nanmean(x(5:10)), L2_str.BRT);


%%
% Zscored measure
fluency(fluency > 250) = nan;  fluency_mal(fluency_mal > 250) = nan;  fluency_tel(fluency_tel > 250) = nan; fluency_eng(fluency_eng > 200) = nan;
fluency = 300-fluency; fluency_mal = 300-fluency_mal; fluency_tel = 300-fluency_tel; 
fluency_fmri = 300-fluency_fmri; 
f1 = (fluency(1:8)  - nanmean(fluency(1:8)))/nanstd(fluency(1:8));
f2 = (fluency(9:16) - nanmean(fluency(9:16)))/nanstd(fluency(9:16));
f3 = (fluency_mal   - nanmean(fluency_mal))/nanstd(fluency_mal);
f4 = (fluency_tel   - nanmean(fluency_tel))/nanstd(fluency_tel);
% fe = (fluency_eng   - nanmean(fluency_eng))/nanstd(fluency_eng);
f5 = (fluency_fmri(1:7)  - nanmean(fluency_fmri(1:7)))/nanstd(fluency_fmri(1:7));
f6 = (fluency_fmri(8:11)  - nanmean(fluency_fmri(8:11)))/nanstd(fluency_fmri(8:11));

% fluency_all = [f1 f2 f3' f4' fe'];
% interaction = [zscore(bmm); zscore(btt); zscore(bm); zscore(bt); zscore(be)];
% [r, p] = partialcorri(fluency_all',interaction,'rows','complete');

fluency_all = [f1 f2 f3' f4' f5 f6];
interaction = [zscore(bmm); zscore(btt); zscore(bm); zscore(bt); zscore(Fbmm); zscore(Fbtt)];
[r,p] = partialcorri(fluency_all',interaction,'rows','complete','type','spearman');

[b, bint] = regress(fluency_all', interaction);
pred = interaction*b;
figure; corrplot(pred,fluency_all',[],1); xlabel('Predicted fluency'); ylabel('Observed fluency'); 
[r(5) p(5)] = corr(fluency_all',pred,'rows','complete','type','spearman');

figure; bar(r)
return

corrplot([f1 f2 f3' f4' f5 f6],[zscore(baseline1) zscore(baseline2m') zscore(baseline2t') zscore(baseline3)])
% fluency_all = [f1 f3' f5] ;
% interaction = [zscore(bmm);zscore(bm); zscore(Fbmm)];
% [r, p] = partialcorri(fluency_all',interaction,'rows','complete','type','spearman');
% 
% 
% fluency_all = [f2 f4' f6];
% interaction = [zscore(btt); zscore(bt); zscore(Fbtt)];
% [r, p] = partialcorri(fluency_all',interaction,'rows','complete','type','spearman');

%%
% fluency_all = [fluency(1:8) fluency_mal']; fluency_all(fluency_all > 250) = nan;
% interaction = [bmm; bm];
% [rm, pm] = partialcorri(fluency_all',interaction,'rows','complete');
% 
% fluency_all = [fluency(9:16) fluency_tel']; fluency_all(fluency_all > 250) = nan;
% interaction = [btt; bt];
% [rt, pt] = partialcorri(fluency_all',interaction,'rows','complete');
% 
% fluency_all = fluency_eng; fluency_all(fluency_all > 200) = nan;
% interaction = be;
% qidx = setdiff(1:numel(fluency_all), find(isnan(fluency_all)));
% [re, pe] = partialcorri(fluency_all,interaction,'rows','complete');

fluency_all = [fluency(1:8) fluency(9:16) fluency_mal' fluency_tel']; fluency_all(fluency_all > 250) = nan;
interaction = [bmm; btt; bm; bt];
[rmt, pmt] = partialcorri(fluency_all',interaction,'rows','complete');

