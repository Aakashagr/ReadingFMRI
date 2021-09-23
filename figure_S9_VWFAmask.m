allclearL2
ismal = [1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 0 0 1 1 1 0 1 0 0 0 0 0 0 1 0 1 1 1 1 1];


mask = zeros(53,63,52);
 for sub = find(ismal == 0)
     xx = spm_read_vols(spm_vol(['VWFA_mask/sub' num2str(sub,'%02d') '_vwfa.nii']));
     mask = mask + xx;     
 end
V = spm_vol('VWFA.nii'); V.fname = 'VWFA_tel.nii';
spm_write_vol(V,mask);


mask = zeros(53,63,52);
 for sub = find(ismal == 1)
     xx = spm_read_vols(spm_vol(['VWFA_mask/sub' num2str(sub,'%02d') '_vwfa.nii']));
     mask = mask + xx;     
 end
V = spm_vol('VWFA.nii'); V.fname = 'VWFA_mal.nii';
spm_write_vol(V,mask);