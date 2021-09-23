allclearL2
% Load data structure
% load L2fmri_READINGw
% age = std([24 24 25 22 20 24 25 24 27 29 25 24 22 23 26 25 34 26 23 25 26 25 26 26 32 32 23 30 24 21 21 30 23 23 19])
ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli
[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions
%% Extracting Peak ROI location

for roi = 1:5
    for sub = 1:35
        [x,y,z] = ind2sub([53 63 52], ids{sub,roi}(1));
        idx(roi,sub,:) = cor2mni([x,y,z]);
    end
end

std(idx(:,:,1),[],2)