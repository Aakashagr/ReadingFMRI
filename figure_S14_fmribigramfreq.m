allclearL2
% Load data structure
% load L2fmri_READINGw

ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli
qmal_hf = [13,21:27,31:34]; qmal_lf = setdiff(qm(qd), qmal_hf);
qtel_hf = [47,56,57,59,60,61,65,66,68,54,45,46]; qtel_lf = setdiff(qt(qd), qtel_hf);
qmal_hf = qmal_hf - 10; qmal_lf = qmal_lf - 10;
qtel_hf = qtel_hf - 44; qtel_lf = qtel_lf - 44;

[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions

% loading letter frequency
load tel_lf; load mal_lf;
lid = [13 17 18 20 23 25 27 29 33 35];
tel_lf = tel_lf(lid); mal_lf = mal_lf(lid);
tel_lf = tel_lf./sum(tel_lf); mal_lf = mal_lf./sum(mal_lf);
m_lfpair = mal_lf(L2_str.letterpairs);
t_lfpair = tel_lf(L2_str.letterpairs);

% loading bigram frequency
load tel_bigram; load mal_bigram
lr = find(L2_str.letterpairs(:,1) < L2_str.letterpairs(:,2));
rl = find(L2_str.letterpairs(:,1) > L2_str.letterpairs(:,2));
for i = 1:24
    qfreq(i) = find(ismember(nchoosek(1:36,2),[lid(L2_str.letterpairs(i,:));fliplr(lid(L2_str.letterpairs(i,:)))],'rows'));
end
tel_bf(lr) = tel_bigram(qfreq(lr),1);  tel_bf(rl) = tel_bigram(qfreq(rl),2); 
mal_bf(lr) = mal_bigram(qfreq(lr),1);  mal_bf(rl) = mal_bigram(qfreq(rl),2);
tel_bf = tel_bf'/sum(tel_bf);
mal_bf = mal_bf'/sum(mal_bf);

% building regression matrix
tel_MAT = [t_lfpair tel_bf ones(24,1)];
mal_MAT = [m_lfpair mal_bf ones(24,1)];

%% Bigram frequency plots
cnt = 0;
for roi = 1:4
    cnt = cnt + 1;
    for sub = 1:numel(ismal)
        if ismal(sub); N = qm(qd); NN = qt(qd); Xnn = tel_MAT; Xn = mal_MAT; else, N = qt(qd); NN = qm(qd); Xn = tel_MAT; Xnn = mal_MAT; end % Identifying native (N) and non-native (NN) stimuli index
        betas = L2_str.mergedevtbeta{sub};                      % Extracting Beta (or regression weights) for each subject
        if roi == 4; max_vox = 20; elseif roi == 3, max_vox = 200; else, max_vox = Inf; end     % Restricting VWFA definitions upto top 75 voxels
        nvox = min(numel(ids{sub,roi}),max_vox);                % Total number of voxels considered in the analysis
        
        % Averaging betas across single and two letter stimuli; separately for native and nonnative stimuli
        actn(sub,:)  = nanmean(betas(ids{sub,roi}(1:nvox),N));
        actnn(sub,:) = nanmean(betas(ids{sub,roi}(1:nvox),NN));
        
        % Fitting frequncy model to predict fMRI responses
%         w = regress(vec(actn(sub,:)),Xn); pred = Xn*w; [rho_n(sub,roi) pval_n(sub,roi)] = nancorrcoef(vec(actn(sub,:)),pred);
%         w = regress(vec(actnn(sub,:)),Xnn); pred = Xnn*w; [rho_nn(sub,roi) pval_nn(sub,roi)] = nancorrcoef(vec(actnn(sub,:)),pred);
        [rho_n(sub,roi) pval_n(sub,roi)] = nancorrcoef(vec(actn(sub,:)),Xn(:,3));
        [rho_nn(sub,roi) pval_nn(sub,roi)] = nancorrcoef(vec(actnn(sub,:)),Xn(:,3));
        
        % partial correlation results
        w = regress(vec(actn(sub,:)),Xn(:,[1 2 4])); pred = Xn(:,[1 2 4])*w; [Prho_n(sub,roi,:) Ppval_n(sub,roi,:)] = partialcorri(vec(actn(sub,:)),[pred Xn(:,3)]);
        Eactn(sub,:) = actn(sub,:) - pred';
        w = regress(vec(actnn(sub,:)),Xnn(:,[1 2 4])); pred = Xnn(:,[1 2 4])*w; [Prho_nn(sub,roi,:) Ppval_nn(sub,roi,:)] = partialcorri(vec(actnn(sub,:)),[pred Xnn(:,3)]);
        Eactnn(sub,:) = actnn(sub,:) - pred';        
        
    end
%     
%     meanacthf = [nanmean(vec(actn(ismal==0,qtel_hf))) nanmean(vec(actnn(ismal==1,qtel_hf)));  nanmean(vec(actn(ismal==0,qtel_lf))) nanmean(vec(actnn(ismal==1,qtel_lf)))];
%     semacthf = [nansem(vec(actn(ismal==0,qtel_hf))) nansem(vec(actnn(ismal==1,qtel_hf))); nansem(vec(actn(ismal==0,qtel_lf))) nansem(vec(actnn(ismal==1,qtel_lf)))];
%     figure(1); subplot(2,4,roi); barweb(meanacthf',semacthf');
%     set(gca,'XTickLabel',{'Readers','Non-readers'}); title(ROIname{roi}); ylabel('Average Betas');
%     global_title([],[],'Telugu Bigrams');
%     
%     meanactlf = [nanmean(vec(actn(ismal==1,qmal_hf))) nanmean(vec(actnn(ismal==0,qmal_hf)));  nanmean(vec(actn(ismal==1,qmal_lf))) nanmean(vec(actnn(ismal==0,qmal_lf)))];
%     semactlf = [nansem(vec(actn(ismal==1,qmal_hf))) nansem(vec(actnn(ismal==0,qmal_hf)));  nansem(vec(actn(ismal==1,qmal_lf))) nansem(vec(actnn(ismal==0,qmal_lf)))];
%     figure(1); subplot(2,4,roi+4); barweb(meanactlf',semactlf');
%     set(gca,'XTickLabel',{'Readers','Non-readers'}); title(ROIname{roi});  ylabel('Average Betas');
%     global_title(.1,.45,'Malayalam Bigrams');
%     
    Tactdiff(roi,:) = [nanmean(nanmean(actn(ismal==0,qtel_lf),2)- nanmean(actn(ismal==0,qtel_hf),2)); nanmean(nanmean(actnn(ismal==1,qtel_lf),2)- nanmean(actnn(ismal==1,qtel_hf),2))];
    Tsemdiff(roi,:) = [nansem(nanmean(actn(ismal==0,qtel_lf),2)- nanmean(actn(ismal==0,qtel_hf),2)); nansem(nanmean(actnn(ismal==1,qtel_lf),2)- nanmean(actnn(ismal==1,qtel_hf),2))];
    Tp(roi) = statcomparemean(nanmean(actn(ismal==0,qtel_lf),2)- nanmean(actn(ismal==0,qtel_hf),2), nanmean(actnn(ismal==1,qtel_lf),2)- nanmean(actnn(ismal==1,qtel_hf),2));
    Mactdiff(roi,:) = [nanmean(nanmean(actn(ismal==1,qmal_lf),2)- nanmean(actn(ismal==1,qmal_hf),2)); nanmean(nanmean(actnn(ismal==0,qmal_lf),2)- nanmean(actnn(ismal==0,qmal_hf),2))];
    Msemdiff(roi,:) = [nansem(nanmean(actn(ismal==1,qmal_lf),2)- nanmean(actn(ismal==1,qmal_hf),2)); nansem(nanmean(actnn(ismal==0,qmal_lf),2)- nanmean(actnn(ismal==0,qmal_hf),2))];
    Mp(roi) = statcomparemean(nanmean(actn(ismal==1,qmal_lf),2)- nanmean(actn(ismal==1,qmal_hf),2), nanmean(actnn(ismal==0,qmal_lf),2)- nanmean(actnn(ismal==0,qmal_hf),2));
    
    
    ETactdiff(roi,:) = [nanmean(nanmean(Eactn(ismal==0,qtel_lf),2)- nanmean(Eactn(ismal==0,qtel_hf),2)); nanmean(nanmean(Eactnn(ismal==1,qtel_lf),2)- nanmean(Eactnn(ismal==1,qtel_hf),2))];
    ETsemdiff(roi,:) = [nansem(nanmean(Eactn(ismal==0,qtel_lf),2)- nanmean(Eactn(ismal==0,qtel_hf),2)); nansem(nanmean(Eactnn(ismal==1,qtel_lf),2)- nanmean(Eactnn(ismal==1,qtel_hf),2))];
    ETp(roi) = statcomparemean(nanmean(Eactn(ismal==0,qtel_lf),2)- nanmean(Eactn(ismal==0,qtel_hf),2), nanmean(Eactnn(ismal==1,qtel_lf),2)- nanmean(Eactnn(ismal==1,qtel_hf),2));
    EMactdiff(roi,:) = [nanmean(nanmean(Eactn(ismal==1,qmal_lf),2)- nanmean(Eactn(ismal==1,qmal_hf),2)); nanmean(nanmean(Eactnn(ismal==0,qmal_lf),2)- nanmean(Eactnn(ismal==0,qmal_hf),2))];
    EMsemdiff(roi,:) = [nansem(nanmean(Eactn(ismal==1,qmal_lf),2)- nanmean(Eactn(ismal==1,qmal_hf),2)); nansem(nanmean(Eactnn(ismal==0,qmal_lf),2)- nanmean(Eactnn(ismal==0,qmal_hf),2))];
    EMp(roi) = statcomparemean(nanmean(Eactn(ismal==1,qmal_lf),2)- nanmean(Eactn(ismal==1,qmal_hf),2), nanmean(Eactnn(ismal==0,qmal_lf),2)- nanmean(Eactnn(ismal==0,qmal_hf),2));

    disp(roi);
end
% figure(1); legend('High frequency bigram','Low frequency bigram','Location','Best');

figure;
subplot(211); barweb(Tactdiff,Tsemdiff); set(gca,'XTickLabel',ROIname); ylabel('(LF - HF) Betas'); title('Telugu Bigrams');
subplot(212); barweb(Mactdiff,Msemdiff); set(gca,'XTickLabel',ROIname); ylabel('(LF - HF) Betas'); title('Malayalam Bigrams');
legend('Readers','Non-readers')

figure;
subplot(211); barweb(ETactdiff,ETsemdiff); set(gca,'XTickLabel',ROIname); ylabel('(LF - HF) Betas'); title('Telugu Bigrams');
subplot(212); barweb(EMactdiff,EMsemdiff); set(gca,'XTickLabel',ROIname); ylabel('(LF - HF) Betas'); title('Malayalam Bigrams');
legend('Readers','Non-readers')


%% comparing model fits across roi
% figure; 
% barweb([nanmean(rho_n(ismal==0,:));nanmean(rho_nn(ismal==1,:))]', [nansem(rho_n(ismal==0,:));nansem(rho_nn(ismal==1,:))]'); 
% set(gca,'XTickLabel',{'V1-V3','V4','LOC','VWFA'}); legend('readers','nonreaders') 
% title('Telugu bigrams'); ylabel('Correlation Coefficient')
% median(pval_n(ismal==0,:))
% 
% 
% figure; 
% barweb([nanmean(rho_n(ismal==1,:));nanmean(rho_nn(ismal==0,:))]', [nansem(rho_n(ismal==1,:));nansem(rho_nn(ismal==0,:))]'); 
% set(gca,'XTickLabel',{'V1-V3','V4','LOC','VWFA'}); legend('readers','nonreaders') 
% title('Malayalam bigrams'); ylabel('Correlation Coefficient')
% %% comparing parital correlation after regressing out letter frequency predictions
% figure; 
% barweb([nanmean(Prho_n(ismal==0,:,2));nanmean(Prho_nn(ismal==1,:,2))]', [nansem(Prho_n(ismal==0,:,2));nansem(Prho_nn(ismal==1,:,2))]'); 
% set(gca,'XTickLabel',{'V1-V3','V4','LOC','VWFA'}); legend('readers','nonreaders') 
% title('Telugu bigrams'); ylabel('Correlation Coefficient')
% median(pval_n(ismal==0,:))
% 
% 
% figure; 
% barweb([nanmean(Prho_n(ismal==1,:,2));nanmean(Prho_nn(ismal==0,:,2))]', [nansem(Prho_n(ismal==1,:,2));nansem(Prho_nn(ismal==0,:,2))]'); 
% set(gca,'XTickLabel',{'V1-V3','V4','LOC','VWFA'}); legend('readers','nonreaders') 
% title('Malayalam bigrams'); ylabel('Correlation Coefficient')
