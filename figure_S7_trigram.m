allclear
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
PC_mal = PC(:,L2_str.subjinfo.ismalayalam == 1,:);
PC_tel = PC(:,L2_str.subjinfo.ismalayalam == 0,:);

mRT_mal = nanmean(RT_mal,3);
mRT_tel = nanmean(RT_tel,3);

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);
RTmm = RT_mal_mean(mal);
RTmt = RT_mal_mean(tel);
RTtt = RT_tel_mean(tel);
RTtm = RT_tel_mean(mal);

% Identifying word-word search pairs
words_tel = [1 4 5; 1 4 3; 1 5 2; 1 3 2; 1 3 3;1 6 4; 2 4 3; 4 2 5; 4 5 2;4 6 2; 6 4 5; 6 5 6 ];
words_mal = [1 4 1; 2 4 1; 2 4 3; 2 4 4; 3 1 6; 4 6 1; 5 1 6; 5 5 3; 6 2 3  ];

% 3 letter transposition stimuli
trans1_tel = perms([6 2 1]); trans1_mal = perms([3 2 1]);
trans2_tel = perms([5 2 4]); trans2_mal = perms([5 1 6]);

% identifying pair locations in object parts list
% telugu 
words_telid = find(ismember(obj_parts,words_tel,'rows') == 1);
trans1_telid = find(ismember(obj_parts,trans1_tel,'rows') == 1);
trans2_telid = find(ismember(obj_parts,trans2_tel,'rows') == 1);

% Malayalam
words_malid = find(ismember(obj_parts,words_mal,'rows') == 1);
trans1_malid = find(ismember(obj_parts,trans1_mal,'rows') == 1);
trans2_malid = find(ismember(obj_parts,trans2_mal,'rows') == 1);

% word search pairs to be included in the stimulus set
% Telugu
words_telpairs = words_telid(nchoosek(1:numel(words_telid),2))+numel(images)/2;
trans1_telpairs = trans1_telid(nchoosek(1:numel(trans1_telid),2))+ numel(images)/2;
trans2_telpairs = trans2_telid(nchoosek(1:numel(trans2_telid),2))+ numel(images)/2;

% Malayalam
words_malpairs = words_malid(nchoosek(1:numel(words_malid),2));
trans1_malpairs = trans1_malid(nchoosek(1:numel(trans1_malid),2));
trans2_malpairs = trans2_malid(nchoosek(1:numel(trans2_malid),2));

% identifying whether seach pairs are already included in the list
% Telugu
words_tpairs_id = find(ismember(t_image_pairs, words_telpairs,'rows') == 1); 
trans1_tpairs_id = find(ismember(t_image_pairs, trans1_telpairs,'rows') == 1); 
trans2_tpairs_id = find(ismember(t_image_pairs, trans2_telpairs,'rows') == 1); 
nw_tpairs_id = setdiff(1:325, words_tpairs_id);

% Malayalam
words_mpairs_id = find(ismember(m_image_pairs, words_malpairs,'rows') == 1); 
trans1_mpairs_id = find(ismember(m_image_pairs, trans1_malpairs,'rows') == 1); 
trans2_mpairs_id = find(ismember(m_image_pairs, trans2_malpairs,'rows') == 1); 
nw_mpairs_id = setdiff(1:325, words_mpairs_id);

% Baseline reaction time
for i = 1:numel(L2_str.BRT)
    brtmean(i,1) = mean(L2_str.BRT{i}(4:10));
end

%% panel B- Summary plot

figure; barweb([nanmean(RTtt), nanmean(RTmt); nanmean(RTtm), nanmean(RTmm); nanmean(brtmean(ismal==0)) nanmean(brtmean(ismal ==1))], ...
    [nansem(nanmean(mRT_tel(tel,:))'), nansem(nanmean(mRT_mal(tel,:))'); ...
    nansem(nanmean(mRT_tel(mal,:))'), nansem(nanmean(mRT_mal(mal,:))'); nansem(brtmean(ismal==0)) nansem(brtmean(ismal ==1))])
ylabel('Avg Search Time, Seconds');   colormap(get(gca,'ColorOrder'));
set(gca,'Xticklabel',{'Telugu Bigrams','Malayalam Bigrams','Baseline stimuli'});
legend({'Telugu readers','Malayalam readers'});
%% Panel F & G -Baton model - Telugu

% full model
% part_pos = [7,9,1,4,6,8,4,2,5,6,5,3,7,9,8];

% reduced model
% part_pos = [6,7,1,4,5,6,4,2,4,5,4,3,6,7,6];

% reduced model - 2
part_pos = [4,5,1,2,3,4,2,1,2,3,2,1,4,5,4];

parm_len = max(part_pos);


comb = nchoosek(1:2*img_len,2); % number of possible part pairs
part_pairs = nchoosek(1:nimg,2); 
mat_tel = zeros(size(t_image_pairs,1),size(part_pairs,1)*parm_len+1);
mat_tel(:,end) = 1;
xx = nchoosek(nimg,2);

for  i = 1:size(t_image_pairs,1)
   
    q1 = [obj_parts(t_image_pairs(i,1)- numel(images)/2,:), obj_parts(t_image_pairs(i,2)- numel(images)/2,:)];
    
    for  id = 1:numel(part_pos)
        [~,idw1] = ismember([q1(comb(id,1)), q1(comb(id,2))], part_pairs,'rows');
        [~,idw2] = ismember([q1(comb(id,2)), q1(comb(id,1))], part_pairs,'rows');            

        if idw1>0; mat_tel(i,xx*(part_pos(id)-1)+ idw1)  = mat_tel(i,xx*(part_pos(id)-1)+ idw1) + 1; end    
        if idw2>0; mat_tel(i,xx*(part_pos(id)-1)+ idw2)  = mat_tel(i,xx*(part_pos(id)-1)+ idw2) + 1; end
 
    end

end

% for Telugu letters
w_tt = regress(1./RTtt,mat_tel); y_pred_tt = mat_tel*w_tt;
w_mt = regress(1./RTmt,mat_tel); y_pred_mt = mat_tel*w_mt;
% sigmoid fit
[~,y_pred_new_tt] = fitsigmoid(y_pred_tt,1./RTtt);
[~,y_pred_new_mt] = fitsigmoid(y_pred_mt,1./RTmt);

% summary plot
figure; barweb([mean(reshape(w_tt(1:end-1), [15 parm_len]))', mean(reshape(w_mt(1:end-1), [15 parm_len]))'],...
               [nansem(reshape(w_tt(1:end-1), [15 parm_len]))',nansem(reshape(w_mt(1:end-1), [15 parm_len]))']);  
xlabel('Parameters');  ylabel('Mean Estimated coefficient')
if parm_len == 9
    set(gca,'Xticklabel',{'C1','C2','C3','A-near1', 'A-near2', 'A-far','W-near1', 'W-near2', 'W-far'})
elseif parm_len == 5
    set(gca,'Xticklabel',{'Corr.','A-near', 'A-far','W-near', 'W-far'});  title('Telugu stimuli')
end
legend('Native readers','Non-native readers');title('Telugu stimuli')

% Baton model- Malayalam
comb = nchoosek(1:2*img_len,2); % number of possible part pairs
part_pairs = nchoosek(1:nimg,2); 
mat_mal = zeros(size(m_image_pairs,1),size(part_pairs,1)*parm_len+1);
mat_mal(:,end) = 1;
xx = nchoosek(nimg,2);

for  i = 1:size(m_image_pairs,1)
   
    q1 = [obj_parts(m_image_pairs(i,1),:), obj_parts(m_image_pairs(i,2),:)];
    
    for  id = 1:numel(part_pos)
        [~,idw1] = ismember([q1(comb(id,1)), q1(comb(id,2))], part_pairs,'rows');
        [~,idw2] = ismember([q1(comb(id,2)), q1(comb(id,1))], part_pairs,'rows');            

        if idw1>0; mat_mal(i,xx*(part_pos(id)-1)+ idw1)  = mat_mal(i,xx*(part_pos(id)-1)+ idw1) + 1; end    
        if idw2>0; mat_mal(i,xx*(part_pos(id)-1)+ idw2)  = mat_mal(i,xx*(part_pos(id)-1)+ idw2) + 1; end
 
    end

end

w_tm = regress(1./RTtm,mat_mal); y_pred_tm = mat_mal*w_tm; 
w_mm = regress(1./RTmm,mat_mal); y_pred_mm = mat_mal*w_mm;

[~,y_pred_new_tm] = fitsigmoid(y_pred_tm,1./RTtm);
[~,y_pred_new_mm] = fitsigmoid(y_pred_mm,1./RTmm);

% summary plot
figure; barweb([mean(reshape(w_mm(1:end-1), [15 parm_len]))', mean(reshape(w_tm(1:end-1), [15 parm_len]))'],...
               [nansem(reshape(w_mm(1:end-1), [15 parm_len]))',nansem(reshape(w_tm(1:end-1), [15 parm_len]))']);  
xlabel('Parameters');  ylabel('Mean Estimated coefficient')
if parm_len == 9
    set(gca,'Xticklabel',{'C1','C2','C3','A-near1', 'A-near2', 'A-far','W-near1', 'W-near2', 'W-far'})
elseif parm_len == 5
    set(gca,'Xticklabel',{'Corr.','A-near', 'A-far','W-near', 'W-far'})
end
legend('Native readers','Non-native readers'); title('Malayalam stimuli')

%%
ids = 61:75;
pt = statcomparemean(w_tt(ids),w_mt(ids));
pm = statcomparemean(w_mm(ids),w_tm(ids));

%% Panel D- model fit
% Telugu on Telugu
figure; corrplot(y_pred_new_tt,1./RTtt,'Tel-Tel',1); hold all;
h(1) = plot(y_pred_new_tt(words_tpairs_id),1./RTtt(words_tpairs_id),'rs'); 
h(2) = plot(y_pred_new_tt([trans1_tpairs_id;trans2_tpairs_id]),1./RTtt([trans1_tpairs_id;trans2_tpairs_id]),'gd'); 
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity'); 
legend(h,'Word-word searches','Transposition searches');

error1 = y_pred_new_tt(words_tpairs_id)- 1./RTtt(words_tpairs_id);
error2 = y_pred_new_tt(nw_tpairs_id)- 1./RTtt(nw_tpairs_id);
error3 = y_pred_new_tt([trans1_tpairs_id;trans2_tpairs_id])- 1./RTtt([trans1_tpairs_id;trans2_tpairs_id]);

%% Pabel E Calculating noise ceiling
[~,cmm] = splithalfcorr(nanmean(RT_mal(mal,:,:),3)'); cmm = spearmanbrowncorrection(cmm,2);
[~,ctt] = splithalfcorr(nanmean(RT_tel(tel,:,:),3)'); ctt = spearmanbrowncorrection(ctt,2);
[~,cmt] = splithalfcorr(nanmean(RT_mal(tel,:,:),3)'); cmt = spearmanbrowncorrection(cmt,2);
[~,ctm] = splithalfcorr(nanmean(RT_tel(mal,:,:),3)'); ctm = spearmanbrowncorrection(ctm,2);


r_tt = nancorrcoef(y_pred_new_tt,1./RTtt);
r_tm = nancorrcoef(y_pred_new_tm,1./RTtm);
r_mt = nancorrcoef(y_pred_new_mt,1./RTmt);
r_mm = nancorrcoef(y_pred_new_mm,1./RTmm);

figure; bar([r_tt, r_mt, r_mm, r_tm]); ylabel('Correlation Coefficient'); hold on
set(gca,'XtickLabel',{'TT', 'MT', 'MM', 'TM'}); title('Model Prediction')
shadedErrorBar([0.8,1.2],[mean(ctt) mean(ctt)], [std(ctt) std(ctt)])
shadedErrorBar([1.8,2.2],[mean(cmt) mean(cmt)], [std(cmt) std(cmt)])
shadedErrorBar([2.8,3.2],[mean(cmm) mean(cmm)], [std(cmm) std(cmm)])
shadedErrorBar([3.8,4.2],[mean(ctm) mean(ctm)], [std(ctm) std(ctm)])

%% Stats
% Telugu & Malayalam letters separately ters
cnt = 0;
for imgpairs = 1:325
    for sub = 1:size(RT,2)
        for reps = 1:size(RT,3)
            cnt = cnt + 1;
            TRT(cnt) = RT(tel(imgpairs),sub,reps);
            MRT(cnt) = RT(mal(imgpairs),sub,reps);
            grp(cnt) = ismal(sub);
            IMGpairs(cnt) = imgpairs;
        end
    end
end
            
pt = anovan(TRT,{grp,IMGpairs},'model','full'); % p < 10-23
pm = anovan(MRT,{grp,IMGpairs},'model','full'); % p < 10-12




