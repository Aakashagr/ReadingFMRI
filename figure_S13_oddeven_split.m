% identifying ROIs using odd/even split for population model and separately
% for Telugu and Malayalam
allclear
load normalised27
pMOD_diff = Rpmod_N - Rpmod_NN;

%%
qm = find(ismal == 1);
qt = find(ismal == 0);
mpMOD_diff(mask == 0) = nan;

pMOD_diff_To = nanmean(pMOD_diff(:,:,:,qt(1:2:end)),4); pMOD_diff_To(mask == 0) = nan;
pMOD_diff_Te = nanmean(pMOD_diff(:,:,:,qt(2:2:end)),4); pMOD_diff_Te(mask == 0) = nan;
pMOD_diff_Mo = nanmean(pMOD_diff(:,:,:,qm(1:2:end)),4); pMOD_diff_Mo(mask == 0) = nan;
pMOD_diff_Me = nanmean(pMOD_diff(:,:,:,qm(2:2:end)),4); pMOD_diff_Me(mask == 0) = nan;
pMOD_diff_T = nanmean(pMOD_diff(:,:,:,qt),4); pMOD_diff_T(mask == 0) = nan;
pMOD_diff_M = nanmean(pMOD_diff(:,:,:,qm),4); pMOD_diff_M(mask == 0) = nan;
pMOD_diff_O = nanmean(pMOD_diff(:,:,:,[qt(1:2:end); qm(1:2:end)]),4); pMOD_diff_O(mask == 0) = nan;
pMOD_diff_E = nanmean(pMOD_diff(:,:,:,[qt(2:2:end); qm(2:2:end)]),4); pMOD_diff_E(mask == 0) = nan;

V = spm_vol('spmT_0001.nii');
V = rmfield(V,'pinfo');
V.fname = 'roidef_figure5\pMOD_diff_To.nii';   spm_write_vol(V,pMOD_diff_To);
V.fname = 'roidef_figure5\pMOD_diff_Te.nii';   spm_write_vol(V,pMOD_diff_Te);
V.fname = 'roidef_figure5\pMOD_diff_Mo.nii';   spm_write_vol(V,pMOD_diff_Mo);
V.fname = 'roidef_figure5\pMOD_diff_Me.nii';   spm_write_vol(V,pMOD_diff_Me);
V.fname = 'roidef_figure5\pMOD_diff_T.nii';   spm_write_vol(V,pMOD_diff_T);
V.fname = 'roidef_figure5\pMOD_diff_M.nii';   spm_write_vol(V,pMOD_diff_M);
V.fname = 'roidef_figure5\pMOD_diff_O.nii';   spm_write_vol(V,pMOD_diff_O);
V.fname = 'roidef_figure5\pMOD_diff_E.nii';   spm_write_vol(V,pMOD_diff_E);
