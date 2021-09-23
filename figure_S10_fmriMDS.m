allclearL2
% Load data structure
% load L2fmri_READINGw

ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli' 
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli
[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions

%%
distmsr = 'spearman'; % 1- r (correlation coefficient) is used as a dissimilarity measure. 
for roi = 1:5
    for sub = 1:numel(ismal)
        betas = L2_str.mergedevtbeta{sub};                      % Extracting Beta (or regression weights) for each subject    
        if roi == 4; max_vox = 20; else, max_vox = Inf; end     % Restricting VWFA definitions upto top 75 voxels
        nvox = min(numel(ids{sub,roi}),max_vox);                % Total number of voxels considered in the analysis
        
        % Calculating pair-wise dissimilarities for single and double letter stimuli
        xx = betas(ids{sub,roi}(1:nvox),1:68)';  xx(:,isnan(mean(xx))) = []; dis(:,sub,roi)  = pdist(xx,distmsr);  
    end
end

%% RDM for each ROI
for roi = 1:5
    figure; imagesc(squareform(squeeze(mean(dis(:,ismal == 1,roi),2)))); colormap(jet); colorbar; title(ROIname{roi}); caxis([0,1.1]);
%     saveas(gcf,[ROIname{roi}, '.svg'])
end

%% For Telugu subjects
for roi = 1:5
    figure; imagesc(squareform(squeeze(mean(dis(:,ismal == 0,roi),2)))); colormap(jet); colorbar; title(ROIname{roi}); caxis([0,1.1]);
%     saveas(gcf,[ROIname{roi}, '.svg'])
end
