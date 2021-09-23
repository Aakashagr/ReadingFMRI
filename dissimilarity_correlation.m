load L2fmri_READINGw.mat
%%
allclearL2
ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli
[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions
%%
distmsr = 'spearman'; % 1- r (correlation coefficient) is used as a dissimilarity measure.
% Visual dissimilarity based on model prediction
dpsyTT = L2_str.search.pdouble.dTT; % Telugu    readers on Telugu    bigrams
dpsyMT = L2_str.search.pdouble.dMT; % Malayalam readers on Telugu    bigrams
dpsyTM = L2_str.search.pdouble.dTM; % Telugu    readers on Malayalam bigrams
dpsyMM = L2_str.search.pdouble.dMM; % Malayalam readers on Malayalam bigrams

% Extracting dissimilarity for each subject and each ROI
for rep = 1
    Subjsample = 1:35;
    for roi = 1:5
        Mcnt = 0; Tcnt = 0;
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
            else
                Tcnt = Tcnt + 1;
                xx = betas(ids{sub,roi}(1:nvox),N(qd))';  xx(:,isnan(mean(xx))) = []; Tdisn(:,Tcnt)  = pdist(xx,distmsr);
                xx = betas(ids{sub,roi}(1:nvox),NN(qd))'; xx(:,isnan(mean(xx))) = []; Tdisnn(:,Tcnt) = pdist(xx,distmsr);
            end
        end
        tdisn(:,roi) = nanmean(Tdisn,2);  tdisnn(:,roi) = nanmean(Tdisnn,2);
        mdisn(:,roi) = nanmean(Mdisn,2);  mdisnn(:,roi) = nanmean(Mdisnn,2);
    end
end

[R_tn P_tn] = corr(tdisn); R_tn(P_tn > 0.06) = NaN;
[R_tnn P_tnn] = corr(tdisnn);  R_tnn(P_tnn > 0.06) = NaN;
[R_mn P_mn] = corr(mdisn); R_mn(P_mn > 0.06) = NaN;
[R_mnn P_mnn] = corr(mdisnn); R_mnn(P_mnn > 0.06) = NaN;

figure; subplot(221); imagesc(R_tn); colorbar
subplot(222); imagesc(R_tnn); colorbar
subplot(223); imagesc(R_mn); colorbar
subplot(224); imagesc(R_mnn); colorbar