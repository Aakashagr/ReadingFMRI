%% Search light over normalized brain to look for correlates of behaviour.
% allclear; load L2fmri_READINGw
allclearL2
ismal = L2_str.ismal;
letterparts = L2_str.letterpairs;

qm = 1:34; qt = 35:68;
qs = 1:10; qd = 11:34;

% Visual dissimilarity
d2psyTT = L2_str.search.pdouble.dTT;  d1psyTT = 1./nanmean(nanmean(L2_str.search.single.RT_TT,2),3);
d2psyMT = L2_str.search.pdouble.dMT;  d1psyMT = 1./nanmean(nanmean(L2_str.search.single.RT_MT,2),3);
d2psyTM = L2_str.search.pdouble.dTM;  d1psyTM = 1./nanmean(nanmean(L2_str.search.single.RT_TM,2),3);
d2psyMM = L2_str.search.pdouble.dMM;  d1psyMM = 1./nanmean(nanmean(L2_str.search.single.RT_MM,2),3);


% Initialsing weights for natives and non-natives
ACT_N    = single(NaN(53,63,52,35));
ACT_N1   = single(NaN(53,63,52,35));
ACT_N2   = single(NaN(53,63,52,35));
Rpmod_N  = single(NaN(53,63,52,35));
Rvmod_N  = single(NaN(53,63,52,35));
R_beh_N1 = single(NaN(53,63,52,2));             p_beh_N1 = R_beh_N1;
R_beh_N2 = single(NaN(53,63,52,2));             p_beh_N2 = R_beh_N2;

ACT_NN    = single(NaN(53,63,52,35));
ACT_NN1   = single(NaN(53,63,52,35));
ACT_NN2   = single(NaN(53,63,52,35));
Rpmod_NN  = single(NaN(53,63,52,35));
Rvmod_NN  = single(NaN(53,63,52,35));
R_beh_NN1 = single(NaN(53,63,52,2));             p_beh_NN1 = R_beh_NN1;
R_beh_NN2 = single(NaN(53,63,52,2));             p_beh_NN2 = R_beh_NN2;
disn1 = single(NaN(53,63,52,45,2));
disnn1 = single(NaN(53,63,52,45,2));
disn2 = single(NaN(53,63,52,276,2));
disnn2 = single(NaN(53,63,52,276,2));

for sub = 1:35
    if ismal(sub) == 0
        N = qt; NN = qm; gid = 1;
    else
        N = qm; NN = qt; gid = 2;
    end
    
    beta = reshape(L2_str.mergedevtbeta{sub},[53 63 52 69]);  beta(:,:,:,69) = [];
    mask = reshape(L2_str.mergedevtbetash{sub},[53 63 52]);   mask(~isnan(mask)) = 1;
    [x,y,z] = ind2sub([53 63 52],find(mask >0));
    
    for idx = 1:numel(x)
        x_range = max(x(idx)-1,1):min(x(idx)+1,53);
        y_range = max(y(idx)-1,1):min(y(idx)+1,63);
        z_range = max(z(idx)-1,1):min(z(idx)+1,52);
        
        act = reshape(beta(x_range,y_range,z_range,:),[numel(x_range)*numel(y_range)*numel(z_range),68]);
        act(isnan(nanmean(act,2)) ,:) = [];  % removing voxels which are not included in the mask.
        vox = find(ismember(act,squeeze(beta(x(idx),y(idx),z(idx),:))','rows'));
        % Comparing ACTIVITY
        % Identifying areas where ACTIVITY for Native stimuli is greater than Non-native stimuli
        ACT_N(x(idx),y(idx),z(idx),sub)  = nanmedian(vec(beta(x(idx),y(idx),z(idx),N)));
        ACT_NN(x(idx),y(idx),z(idx),sub) = nanmedian(vec(beta(x(idx),y(idx),z(idx),NN)));
        
        % Identifying areas where ACTIVITY for single letter Native stimuli is greater than Non-native stimuli
        ACT_N1(x(idx),y(idx),z(idx),sub)  = nanmedian(vec(beta(x(idx),y(idx),z(idx),N(qs))));
        ACT_NN1(x(idx),y(idx),z(idx),sub) = nanmedian(vec(beta(x(idx),y(idx),z(idx),NN(qs))));
        
        % Identifying areas where ACTIVITY for two letter Native stimuli is greater than Non-native stimuli
        ACT_N2(x(idx),y(idx),z(idx),sub)  =  nanmedian(vec(beta(x(idx),y(idx),z(idx),N(qd))));
        ACT_NN2(x(idx),y(idx),z(idx),sub) =  nanmedian(vec(beta(x(idx),y(idx),z(idx),NN(qd))));
        
        % Single voxel model
        robsm = squeeze(beta(x(idx),y(idx),z(idx),:));
        qstim = N; robs2  = robsm(qstim(qd)); r1 = robsm(qstim(qs)); robs1 = r1(letterparts);
        X = [robs1 ones(length(robs2),1)]; b = regress(robs2,X); rpred2 = X*b;
        Rvmod_N(x(idx),y(idx),z(idx),sub) = nancorrcoef(robs2,rpred2);
        
        qstim = NN; robs2  = robsm(qstim(qd)); r1 = robsm(qstim(qs)); robs1 = r1(letterparts);
        X = [robs1 ones(length(robs2),1)]; b = regress(robs2,X); rpred2 = X*b;
        Rvmod_NN(x(idx),y(idx),z(idx),sub) = nancorrcoef(robs2,rpred2);
        
        
        if size(act,1) > 10
            % PART BASED MODEL
            act_N  = act(:,N(qd))';
            act_NN = act(:,NN(qd))';
            clear pred_N pred_NN w_N w_NN
            for img = 1:24
                N1 = N(qs);  NN1 = NN(qs);
                
                % Part based activity for individual letters
                act_N1 = act(:,N1(letterparts(img,1)));             act_N2 = act(:,N1(letterparts(img,2)));
                act_NN1 = act(:,NN1(letterparts(img,1)));           act_NN2 = act(:,NN1(letterparts(img,2)));
                mat_N = [act_N1, act_N2, ones(numel(act_N1),1)];    mat_NN = [act_NN1, act_NN2, ones(numel(act_NN1),1)];
                w_N(img,:) = regress(act_N(img,:)',mat_N);          w_NN(img,:) = regress(act_NN(img,:)',mat_NN);
                pred_N(img,:) = mat_N*w_N(img,:)';                  pred_NN(img,:) = mat_NN*w_NN(img,:)';
            end
            
            % Population model predictions
            Rpmod_N(x(idx),y(idx),z(idx),sub)  = nancorrcoef(act_N(:,vox),pred_N(:,vox));
            Rpmod_NN(x(idx),y(idx),z(idx),sub) = nancorrcoef(act_NN(:,vox),pred_NN(:,vox));
            
            disn1(x(idx),y(idx),z(idx),:,gid)  = nansum([squeeze(disn1(x(idx),y(idx),z(idx),:,gid)), pdist(act(:,N(qs))','spearman')'],2);
            disn2(x(idx),y(idx),z(idx),:,gid)  = nansum([squeeze(disn2(x(idx),y(idx),z(idx),:,gid)), pdist(act(:,N(qd))','spearman')'],2);
            disnn1(x(idx),y(idx),z(idx),:,gid) = nansum([squeeze(disnn1(x(idx),y(idx),z(idx),:,gid)), pdist(act(:,NN(qs))','spearman')'],2);
            disnn2(x(idx),y(idx),z(idx),:,gid) = nansum([squeeze(disnn2(x(idx),y(idx),z(idx),:,gid)), pdist(act(:,NN(qd))','spearman')'],2);
        end
    end
    disp(sub)
end

for i = 1:35
    mask = reshape(L2_str.mergedevtbetash{sub},[53 63 52]);   mask(~isnan(mask)) = 1;
    MASK(:,:,:,sub) = mask;
end
mask = nanmean(MASK,4); mask(mask > 0) = 1;


for gid = 1:2
    if gid == 1 % Telugu group
        Nd1 = d1psyTT;  NNd1 = d1psyTM;
        Nd2 = d2psyTT;  NNd2 = d2psyTM;
        N = qt; NN = qm; nsubj = numel(find(ismal == 0));
    else  % Malayalam group
        NNd1 = d1psyMT;  Nd1 = d1psyMM;
        NNd2 = d2psyMT;  Nd2 = d2psyMM;
        N = qm; NN = qt;  nsubj = numel(find(ismal == 1));
    end
    
    [x,y,z] = ind2sub([53 63 52],find(mask >0));
    for idx = 1:numel(x)
        % correlating BEHAVIOUR dissimilarity with voxel
        [R_beh_N1(x(idx),y(idx),z(idx),gid), p_beh_N1(x(idx),y(idx),z(idx),gid)]   = nancorrcoef(Nd1,  disn1(x(idx),y(idx),z(idx),:,gid)/nsubj);
        [R_beh_NN1(x(idx),y(idx),z(idx),gid),p_beh_NN1(x(idx),y(idx),z(idx),gid)]  = nancorrcoef(NNd1, disnn1(x(idx),y(idx),z(idx),:,gid)/nsubj);
        
        [R_beh_N2(x(idx),y(idx),z(idx),gid), p_beh_N2(x(idx),y(idx),z(idx),gid)]   = nancorrcoef(Nd2,  disn2(x(idx),y(idx),z(idx),:,gid)/nsubj);
        [R_beh_NN2(x(idx),y(idx),z(idx),gid),p_beh_NN2(x(idx),y(idx),z(idx),gid)]  = nancorrcoef(NNd2, disnn2(x(idx),y(idx),z(idx),:,gid)/nsubj);
    end
    disp(gid)
end

save normalised27
%%   Extracting and saving the data
R_beh_N1 = nanmean(R_beh_N1,4);   p_beh_N1  = nanmean(p_beh_N1,4);
R_beh_N2 = nanmean(R_beh_N2,4);   p_beh_N2  = nanmean(p_beh_N2,4);
R_beh_NN1 = nanmean(R_beh_NN1,4); p_beh_NN1 = nanmean(p_beh_NN1,4);
R_beh_NN2 = nanmean(R_beh_NN2,4); p_beh_NN2 = nanmean(p_beh_NN2,4);

ACT_diff  = ACT_N - ACT_NN;
pMOD_diff = Rpmod_N - Rpmod_NN;
vMOD_diff = Rvmod_N - Rvmod_NN;

X = NaN(53,63,52);
p_ACT_N = X;    p_ACT_NN = X;    p_ACT_diff = X;
ppR_mod_N = X;  ppR_mod_NN = X;  ppMOD_diff = X;
pvR_mod_N = X;  pvR_mod_NN = X;  pvMOD_diff = X;

for x = 1:53
    for y = 1:63
        for z = 1:52
            if ~isnan(nanmean(squeeze(ACT_N(x,y,z,:))) + nanmean(squeeze(Rpmod_N(x,y,z,:))))
                [p_ACT_N(x,y,z)]     = signrank(squeeze(ACT_N(x,y,z,:)));
                [p_ACT_NN(x,y,z)]    = signrank(squeeze(ACT_NN(x,y,z,:)));
                [p_ACT_diff(x,y,z)]  = signrank(squeeze(ACT_diff(x,y,z,:)));
                
                [ppR_mod_N(x,y,z)]   = signrank(squeeze(Rpmod_N(x,y,z,:)));
                [ppR_mod_NN(x,y,z)]  = signrank(squeeze(Rpmod_NN(x,y,z,:)));
                [ppMOD_diff(x,y,z)]  = signrank(squeeze(pMOD_diff(x,y,z,:)));
                
                [pvR_mod_N(x,y,z)]   = signrank(squeeze(Rvmod_N(x,y,z,:)));
                [pvR_mod_NN(x,y,z)]  = signrank(squeeze(Rvmod_NN(x,y,z,:)));
                [pvMOD_diff(x,y,z)]  = signrank(squeeze(vMOD_diff(x,y,z,:)));
            end
        end
    end
    disp(x)
end
%%
FDR = 1;
if ~FDR
    % Uncorrected
    thresh = .05;   % 1 for raw
    pth_ACT_N     = thresh;
    pth_ACT_NN    = thresh;
    pth_ACT_diff  = thresh;
    
    pth_R_pmod_N   = thresh;
    pth_R_pmod_NN  = thresh;
    pth_pMOD_diff = thresh;
    
    pth_R_vmod_N   = thresh;
    pth_R_vmod_NN  = thresh;
    pth_vMOD_diff = thresh;
    
    pth_R_beh_N2  = thresh;
    pth_R_beh_NN2 = thresh;
    Folder = 'UNC';
else
    % FDR corrected
    thresh = .05;
    pth_ACT_N     = FDRthreshold(p_ACT_N,thresh,mask);
    pth_ACT_NN    = FDRthreshold(p_ACT_NN,thresh,mask);
    pth_ACT_diff  = FDRthreshold(p_ACT_diff,thresh,mask);
    
    pth_R_pmod_N   = FDRthreshold(ppR_mod_N,thresh,mask);
    pth_R_pmod_NN  = FDRthreshold(ppR_mod_NN,thresh,mask);
    pth_pMOD_diff  = FDRthreshold(ppMOD_diff,thresh,mask);
    
    pth_R_vmod_N   = FDRthreshold(pvR_mod_N,thresh,mask);
    pth_R_vmod_NN  = FDRthreshold(pvR_mod_NN,thresh,mask);
    pth_vMOD_diff  = FDRthreshold(pvMOD_diff,thresh,mask);
    
    pth_R_beh_N2  = FDRthreshold(p_beh_N2,thresh,mask);
    pth_R_beh_NN2 = FDRthreshold(p_beh_NN2,thresh,mask);
    Folder = 'FDR';
end

mACT_N = nanmean(ACT_N,4);          mACT_N(p_ACT_N>pth_ACT_N | mask == 0) = nan;
mACT_NN = nanmean(ACT_NN,4);        mACT_NN(p_ACT_NN>pth_ACT_NN | mask == 0) = nan;
mACT_diff = nanmean(ACT_diff,4);    mACT_diff(p_ACT_diff>pth_ACT_diff | mask == 0) = nan;
% mACT_diff = (nanmean(ACT_diff,4)*sqrt(size(ACT_diff,4)))./(nanstd(ACT_diff,[],4)+.1);    mACT_diff(mask == 0) = nan; % Tvalue


mR_pmod_N = nanmean(Rpmod_N,4);      mR_pmod_N(ppR_mod_N>pth_R_pmod_N | mask == 0) = nan;
mR_pmod_NN = nanmean(Rpmod_NN,4);    mR_pmod_NN(ppR_mod_NN>pth_R_pmod_NN | mask == 0) = nan;
mpMOD_diff = nanmean(pMOD_diff,4);   mpMOD_diff(ppMOD_diff>pth_pMOD_diff | mask == 0) = nan;
% mpMOD_diff = (nanmean(pMOD_diff,4)*sqrt(size(pMOD_diff,4)))./(nanstd(pMOD_diff,[],4)+.01);    mpMOD_diff(mask == 0) = nan; % Tvalue
% xtemp_P = zeros(size(mpMOD_diff)); xtemp_N = zeros(size(mpMOD_diff)); 
% xtemp_N(mpMOD_diff <0) = -log(1+abs(mpMOD_diff(mpMOD_diff < 0))); xtemp_P(mpMOD_diff > 0) = log(1+abs(mpMOD_diff(mpMOD_diff > 0)));
% xtemp = xtemp_N + xtemp_P;

mR_vmod_N = nanmean(Rvmod_N,4);      mR_vmod_N(pvR_mod_N>pth_R_vmod_N | mask == 0) = nan;
mR_vmod_NN = nanmean(Rvmod_NN,4);    mR_vmod_NN(pvR_mod_NN>pth_R_vmod_NN | mask == 0) = nan;
mvMOD_diff = nanmean(vMOD_diff,4);    mvMOD_diff(pvMOD_diff>pth_vMOD_diff | mask == 0) = nan;

mR_beh_N2 = nanmean(R_beh_N2,4);    mR_beh_N2(p_beh_N2>pth_R_beh_N2 | mask == 0) = nan;
mR_beh_NN2 = nanmean(R_beh_NN2,4);  mR_beh_NN2(p_beh_NN2>pth_R_beh_NN2 | mask == 0) = nan;

% saving different matrix as a mask
V = spm_vol('..\preprocessing\20170329_READING_SUB05\glm\wdnlocmerged\spmT_0001.nii');
V = rmfield(V,'pinfo');

V.fname = ['searchlight\' Folder '\ACT_N.nii'];      spm_write_vol(V,mACT_N);
V.fname = ['searchlight\' Folder '\ACT_NN.nii'];     spm_write_vol(V,mACT_NN);
V.fname = ['searchlight\' Folder '\ACT_diff.nii'];   spm_write_vol(V,mACT_diff);
V.fname = ['searchlight\' Folder '\pMOD_N.nii'];      spm_write_vol(V,mR_pmod_N);
V.fname = ['searchlight\' Folder '\pMOD_NN.nii'];     spm_write_vol(V,mR_pmod_NN);
V.fname = ['searchlight\' Folder '\pMOD_diff.nii'];   spm_write_vol(V,mpMOD_diff);
V.fname = ['searchlight\' Folder '\vMOD_N.nii'];      spm_write_vol(V,mR_vmod_N);
V.fname = ['searchlight\' Folder '\vMOD_NN.nii'];     spm_write_vol(V,mR_vmod_NN);
V.fname = ['searchlight\' Folder '\vMOD_diff.nii'];   spm_write_vol(V,mvMOD_diff);
V.fname = ['searchlight\' Folder '\BEH_N.nii'];      spm_write_vol(V,mR_beh_N2);
V.fname = ['searchlight\' Folder '\BEH_NN.nii'];     spm_write_vol(V,mR_beh_NN2);

return

%% thresholding data

TACT_diff  = nanmedian((ACT_N(:,:,:,ismal == 0)  - ACT_NN(:,:,:,ismal == 0)),4);
TACT_diff1 = nanmedian((ACT_N1(:,:,:,ismal == 0) - ACT_NN1(:,:,:,ismal == 0)),4);
TACT_diff2 = nanmedian((ACT_N2(:,:,:,ismal == 0) - ACT_NN2(:,:,:,ismal == 0)),4);
Tbeh_diff1 = nanmedian((R_beh_N1(:,:,:,ismal == 0) - R_beh_NN1(:,:,:,ismal == 0)),4);
Tbeh_diff2 = nanmedian((R_beh_N2(:,:,:,ismal == 0) - R_beh_NN2(:,:,:,ismal == 0)),4);
TMOD24_diff =  nanmedian((R_24N(:,:,:,ismal == 0) - R_24NN(:,:,:,ismal == 0)),4);
TMOD_diff  = nanmedian((Rpmod_N(:,:,:,ismal == 0) - Rpmod_NN(:,:,:,ismal == 0)),4);

MACT_diff  = nanmedian((ACT_N(:,:,:,ismal == 1)  - ACT_NN(:,:,:,ismal == 1)),4);
MACT_diff1 = nanmedian((ACT_N1(:,:,:,ismal == 1) - ACT_NN1(:,:,:,ismal == 1)),4);
MACT_diff2 = nanmedian((ACT_N2(:,:,:,ismal == 1) - ACT_NN2(:,:,:,ismal == 1)),4);
Mbeh_diff1 = nanmedian((R_beh_N1(:,:,:,ismal == 1) - R_beh_NN1(:,:,:,ismal == 1)),4);
Mbeh_diff2 = nanmedian((R_beh_N2(:,:,:,ismal == 1) - R_beh_NN2(:,:,:,ismal == 1)),4);
MMOD24_diff =  nanmedian((R_24N(:,:,:,ismal == 1) - R_24NN(:,:,:,ismal == 1)),4);
MMOD_diff  = nanmedian((Rpmod_N(:,:,:,ismal == 1) - Rpmod_NN(:,:,:,ismal == 1)),4);

% saving different matrix as a mask
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = spm_vol(['..\preprocessing\20170329_READING_SUB05\glm\wdnlocmerged\spmT_0001.nii']);
V = rmfield(V,'pinfo');

% ACTIVITY data
V.fname = 'searchlight\TACT_N.nii'; spm_write_vol(V,nanmedian(ACT_N(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TACT_NN.nii';spm_write_vol(V,nanmedian(ACT_NN(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TACT_N1.nii'; spm_write_vol(V,nanmedian(ACT_N1(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TACT_NN1.nii'; spm_write_vol(V,nanmedian(ACT_NN1(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TACT_N2.nii'; spm_write_vol(V,nanmedian(ACT_N2(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TACT_NN2.nii'; spm_write_vol(V,nanmedian(ACT_NN2(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TACT_diff.nii'; spm_write_vol(V,TACT_diff);
V.fname = 'searchlight\TACT_diff1.nii'; spm_write_vol(V,TACT_diff1);
V.fname = 'searchlight\TACT_diff2.nii'; spm_write_vol(V,TACT_diff2);

V.fname = 'searchlight\MACT_N.nii'; spm_write_vol(V,nanmedian(ACT_N(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MACT_NN.nii';spm_write_vol(V,nanmedian(ACT_NN(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MACT_N1.nii'; spm_write_vol(V,nanmedian(ACT_N1(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MACT_NN1.nii'; spm_write_vol(V,nanmedian(ACT_NN1(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MACT_N2.nii'; spm_write_vol(V,nanmedian(ACT_N2(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MACT_NN2.nii'; spm_write_vol(V,nanmedian(ACT_NN2(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MACT_diff.nii'; spm_write_vol(V,MACT_diff);
V.fname = 'searchlight\MACT_diff1.nii'; spm_write_vol(V,MACT_diff1);
V.fname = 'searchlight\MACT_diff2.nii'; spm_write_vol(V,MACT_diff2);

% MODEL DATA
V.fname = 'searchlight\TR_mod_N.nii'; spm_write_vol(V,nanmedian(Rpmod_N(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TR_mod_NN.nii'; spm_write_vol(V,nanmedian(Rpmod_NN(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TR_24N.nii'; spm_write_vol(V,nanmedian(nanmedian(R_24N(:,:,:,ismal == 0,:),4),5));
V.fname = 'searchlight\TR_24NN.nii'; spm_write_vol(V,nanmedian(nanmedian(R_24NN(:,:,:,ismal == 0,:),4),5));
V.fname = 'searchlight\TMOD_diff.nii'; spm_write_vol(V,TMOD_diff);
V.fname = 'searchlight\TMOD24_diff.nii'; spm_write_vol(V,TMOD24_diff);

V.fname = 'searchlight\MR_mod_N.nii'; spm_write_vol(V,nanmedian(Rpmod_N(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MR_mod_NN.nii'; spm_write_vol(V,nanmedian(Rpmod_NN(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MR_24N.nii'; spm_write_vol(V,nanmedian(nanmedian(R_24N(:,:,:,ismal == 1,:),4),5));
V.fname = 'searchlight\MR_24NN.nii'; spm_write_vol(V,nanmedian(nanmedian(R_24NN(:,:,:,ismal == 1,:),4),5));
V.fname = 'searchlight\MMOD_diff.nii'; spm_write_vol(V,MMOD_diff);
V.fname = 'searchlight\MMOD24_diff.nii'; spm_write_vol(V,MMOD24_diff);

% BEHAVIOUR DATA

V.fname = 'searchlight\TBEH_N1.nii'; spm_write_vol(V,nanmedian(R_beh_N1(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TBEH_NN1.nii'; spm_write_vol(V,nanmedian(R_beh_NN1(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TBEH_N2.nii'; spm_write_vol(V,nanmedian(R_beh_N2(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TBEH_NN2.nii'; spm_write_vol(V,nanmedian(R_beh_NN2(:,:,:,ismal == 0),4));
V.fname = 'searchlight\TBEH_diff1.nii'; spm_write_vol(V,Tbeh_diff1);
V.fname = 'searchlight\TBEH_diff2.nii'; spm_write_vol(V,Tbeh_diff2);

V.fname = 'searchlight\MBEH_N1.nii'; spm_write_vol(V,nanmedian(R_beh_N1(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MBEH_NN1.nii'; spm_write_vol(V,nanmedian(R_beh_NN1(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MBEH_N2.nii'; spm_write_vol(V,nanmedian(R_beh_N2(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MBEH_NN2.nii'; spm_write_vol(V,nanmedian(R_beh_NN2(:,:,:,ismal == 1),4));
V.fname = 'searchlight\MBEH_diff1.nii'; spm_write_vol(V,Mbeh_diff1);
V.fname = 'searchlight\MBEH_diff2.nii'; spm_write_vol(V,Mbeh_diff2);

