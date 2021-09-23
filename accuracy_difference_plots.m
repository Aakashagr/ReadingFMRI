allclearL2
% Load data structure
% load L2fmri_READINGw
ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli
[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions
load acc_diff
HD = AccDiff_idx(1:17);
LD = AccDiff_idx(18:35);
%% Plotting mean activity values for different regions
for roi = 1:5 %
    for sub = 1:numel(ismal)
        if ismal(sub); N = qm; NN = qt; else, N = qt; NN = qm; end % Identifying native (N) and non-native (NN) stimuli index
        betas = L2_str.mergedevtbeta{sub};                      % Extracting Beta (or regression weights) for each subject
        if roi == 4; max_vox = 20; elseif roi == 3, max_vox = 200; else, max_vox = Inf; end     % Restricting VWFA definitions upto top 75 voxels
        nvox = min(numel(ids{sub,roi}),max_vox);                % Total number of voxels considered in the analysis
        
        % Averaging betas across single and two letter stimuli; separately for native and nonnative stimuli
        actn(sub,:)  = nanmean(betas(ids{sub,roi}(1:nvox),N));
        actnn(sub,:) = nanmean(betas(ids{sub,roi}(1:nvox),NN));
    end
    
    % Calculating mean and standard error of mean across subjects. Grouping data according to Readers
    HDmean(roi,:) = [nanmean(vec(actn(HD,:))) nanmean(vec(actnn(HD,:)))];
    LDmean(roi,:) = [nanmean(vec(actn(LD,:))) nanmean(vec(actnn(LD,:)))];

    HDsem(roi,:) = [nansem(nanmean(actn(HD,:),2)) nansem(nanmean(actnn(HD,:),2))];
    LDsem(roi,:) = [nansem(nanmean(actn(LD,:),2)) nansem(nanmean(actnn(LD,:),2))];
    
    % pair-wise statistical test
    P(roi,:) =  [signrank(nanmean(actn(HD,:),2), nanmean(actnn(HD,:),2)); signrank(nanmean(actn(LD,:),2), nanmean(actnn(LD,:),2))];
end

figure; subplot(121); barweb(HDmean, HDsem); set(gca,'Xticklabel',{'V1-V3','V4','LOC','VWFA','TG'}); title('High accuracy difference')
subplot(122); barweb(LDmean, LDsem); set(gca,'Xticklabel',{'V1-V3','V4','LOC','VWFA','TG'});  title('Low accuracy difference')
legend('Known stimuli','Unknown stimuli','Location','Best');

%% Modelling activation of Bigrams using activations of single letters
for roi = 1:5 % 3 = LOC
    % Selecting top N voxels to model and initialising the variables
    topnvoxels = 200; N = NaN(numel(ismal),topnvoxels);
    cvn = N; pvn = N; cvnn = N; pvnn = N; bn = NaN(numel(ismal),topnvoxels,3); bnn = bn;
    cpvn = N; ppvn = N; cpvnn = N; ppvnn = N;
    
    clear bpn bpnn
    for sub = 2:35
        if ismal(sub); N = qm; NN = qt; else, N = qt; NN = qm; end
        nvox = min(numel(ids{sub,roi}),topnvoxels);                % Total number of voxels considered in the analysis
        betas = L2_str.mergedevtbeta{sub}(ids{sub,roi}(1:nvox),:);
        betas(isnan(nanmean(betas,2)),:) = [];
        
        % fit a model for every bigram relating parts to wholes
        clear rpredpn rpredpnn
        for bid = 1:24
            % model for native bigram responses
            qstim = N;
            robs2m = betas(:,qstim(qd(bid)));   % Observed activity for a bigram
            p = L2_str.letterpairs(bid,:); % Single letter ID of each bigram
            r1 = [betas(:,qstim(p(1))) betas(:,qstim(p(2)))]; %r1 = sort(r1,2); % Corresponding single letter activity
            X = [r1 ones(size(r1,1),1)]; bp = regress(robs2m,X); rpred2m = X*bp; % predicting betas for each bigram
            [cpn(sub,bid),ppn(sub,bid)] = nancorrcoef(robs2m,rpred2m);
            bpn(sub,bid,:) = bp;  rpredpn(bid,:) = rpred2m;
            
            % model for nonnative bigram responses
            qstim = NN;
            robs2m = betas(:,qstim(qd(bid)));
            p = L2_str.letterpairs(bid,:);
            r1 = [betas(:,qstim(p(1))) betas(:,qstim(p(2)))]; %r1 = sort(r1,2);
            X = [r1 ones(size(r1,1),1)]; bp = regress(robs2m,X); rpred2m = X*bp;
            [cpnn(sub,bid),ppnn(sub,bid)] = nancorrcoef(robs2m,rpred2m);
            bpnn(sub,bid,:) = bp;         rpredpnn(bid,:) = rpred2m;
        end
        
        % calculate correlation of population model predictions for each voxel
        for voxid = 1:size(betas,1)
            qstim = N;
            [c,p] = nancorrcoef(betas(voxid,qstim(qd)),rpredpn(:,voxid));
            cpvn(sub,voxid) = c; ppvn(sub,voxid) = p; % correlation of pop model on voxels for native letters
            
            qstim = NN;
            [c,p] = nancorrcoef(betas(voxid,qstim(qd)),rpredpnn(:,voxid));
            cpvnn(sub,voxid) = c; ppvnn(sub,voxid) = p; % correlation of pop model on voxels for nonnative letters
        end
        
        fprintf('Subject %d \n',sub);
    end
    R_n(roi,:) = nanmean(cpvn,2);
    R_nn(roi,:) = nanmean(cpvnn,2);    
    
end

figure; subplot(121); barweb([nanmean(R_n(:,HD),2) nanmean(R_nn(:,HD),2)], [nansem(R_n(:,HD),2) nansem(R_nn(:,HD),2)]);
set(gca,'Xticklabel',{'V1-V3','V4','LOC','VWFA','TG'}); title('High accuracy difference')
legend('Known', 'Unknown'); ylabel('Correlation Coefficient');

subplot(122); barweb([nanmean(R_n(:,LD),2) nanmean(R_nn(:,LD),2)], [nansem(R_n(:,LD),2) nansem(R_nn(:,LD),2)]);
set(gca,'Xticklabel',{'V1-V3','V4','LOC','VWFA','TG'});  title('Low accuracy difference')
legend('Known', 'Unknown'); ylabel('Correlation Coefficient');

for i = 1:5
    HPmodel(i) = signrank(R_n(i,HD),R_nn(i,HD));
    LPmodel(i) = signrank(R_n(i,LD),R_nn(i,LD));
    LPmodel(i) = signrank(R_n(i,:),R_nn(i,:));

end
%% Match with behaviour
allclearL2;
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli
ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions
load acc_diff
HD = AccDiff_idx(1:17);
LD = AccDiff_idx(18:35);

distmsr = 'spearman'; % 1- r (correlation coefficient) is used as a dissimilarity measure.
% Visual dissimilarity based on model prediction
dpsyTT = L2_str.search.pdouble.dTT; % Telugu    readers on Telugu    bigrams
dpsyMT = L2_str.search.pdouble.dMT; % Malayalam readers on Telugu    bigrams
dpsyTM = L2_str.search.pdouble.dTM; % Telugu    readers on Malayalam bigrams
dpsyMM = L2_str.search.pdouble.dMM; % Malayalam readers on Malayalam bigrams

% Extracting dissimilarity for each subject and each ROI
for rep = 1:1000
    Subjsample = randi(numel(ismal),[numel(ismal),1])';
    for roi = 1:5
        Mcnt = 0; Tcnt = 0; Tdisn =[]; Tdisnn = []; Mdisn =[]; Mdisnn = []; isTHD =[]; isMHD = [];
        for sub = Subjsample
            if ismal(sub); N = qm; NN = qt; else, N = qt; NN = qm; end
            betas = L2_str.mergedevtbeta{sub};                      % Extracting Beta (or regression weights) for each subject
            if roi == 4; max_vox = 20; elseif roi == 3, max_vox = 200; else, max_vox = Inf; end     % Restricting VWFA definitions upto top 75 voxels
            nvox = min(numel(ids{sub,roi}),max_vox);                % Total number of voxels considered in the analysis
            
            % Calculating pair-wise dissimilarities for single and double letter stimuli
            if ismal(sub)
                Mcnt = Mcnt + 1;
                xx = betas(ids{sub,roi}(1:nvox),N(qd))';  xx(:,isnan(mean(xx))) = []; Mdisn(:,Mcnt)  = pdist(xx,distmsr);
                xx = betas(ids{sub,roi}(1:nvox),NN(qd))'; xx(:,isnan(mean(xx))) = []; Mdisnn(:,Mcnt) = pdist(xx,distmsr);
                if sum(ismember(HD,sub)), isMHD(Mcnt) = 1; else, isMHD(Mcnt) = 0; end
            else
                Tcnt = Tcnt + 1;
                xx = betas(ids{sub,roi}(1:nvox),N(qd))';  xx(:,isnan(mean(xx))) = []; Tdisn(:,Tcnt)  = pdist(xx,distmsr);
                xx = betas(ids{sub,roi}(1:nvox),NN(qd))'; xx(:,isnan(mean(xx))) = []; Tdisnn(:,Tcnt) = pdist(xx,distmsr);
                if sum(ismember(HD,sub)), isTHD(Tcnt) = 1; else, isTHD(Tcnt) = 0; end
            end
            
        end
        
        [R_Nhd(roi,rep), P_Nhd(roi,rep)] = nancorrcoef([zscore(dpsyTT) zscore(dpsyMM)] , [zscore(nanmean(Tdisn(:,find(isTHD)),2)) zscore(nanmean(Mdisn(:,find(isMHD)),2))]);
        [R_NNhd(roi,rep), P_NNhd(roi,rep)] = nancorrcoef([zscore(dpsyMT) zscore(dpsyTM)] , [zscore(nanmean(Mdisnn(:,find(isMHD)),2))  zscore(nanmean(Tdisnn(:,find(isTHD)),2))]); 
        
        [R_Nld(roi,rep), P_Nld(roi,rep)] = nancorrcoef([zscore(dpsyTT) zscore(dpsyMM)] , [zscore(nanmean(Tdisn(:,find(~isTHD)),2)) zscore(nanmean(Mdisn(:,find(~isMHD)),2))]);
        [R_NNld(roi,rep), P_NNld(roi,rep)] = nancorrcoef([zscore(dpsyMT) zscore(dpsyTM)] , [zscore(nanmean(Mdisnn(:,find(~isMHD)),2))  zscore(nanmean(Tdisnn(:,find(~isTHD)),2))]); 
    end
    disp(rep)
end

p_N_hd = median(P_Nhd,2); p_NN_hd = median(P_NNhd,2);
p_N_ld = median(P_Nld,2); p_NN_ld = median(P_NNld,2);

figure; subplot(121); barweb([mean(R_Nhd,2), mean(R_NNhd,2)]',[std(R_Nhd,[],2), std(R_NNhd,[],2)]'); colormap(pink); title('High accuracy difference '); legend(ROIname)
set(gca,'Xticklabel',{'Known Stimuli','Unknown Stimuli'}); ylabel('Correlation with Behaviour'); 

subplot(122);  barweb([mean(R_Nld,2), mean(R_NNld,2)]',[std(R_Nld,[],2), std(R_NNld,[],2)]'); colormap(pink); title('Low accuracy difference '); legend(ROIname)
set(gca,'Xticklabel',{'Known Stimuli','Unknown Stimuli'}); ylabel('Correlation with Behaviour'); 
