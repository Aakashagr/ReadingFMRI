% Correlation with reading fluency
allclearL2
% load L2fmri_READINGw
ismal = L2_str.ismal;
fluency = L2_str.fluency;
qm = 1:34; qt = 35:68;
[idx, ROIname, aROIs] = getvoxind(L2_str);
%%
for rep = 1:100
    Subjsample = randi(numel(ismal),[numel(ismal),1])';
    ISMAL = ismal(Subjsample); FLUENCY = fluency(Subjsample);
    
    for roi = 1:5
        cnt = 0; clear tmean
        for sub = Subjsample
            cnt = cnt + 1;
            switch roi
                case 1
                    % for EVC
                    aidx = []; for aR = 1:8, aidx = [aidx; L2_str.aROI.(aROIs{aR}){sub}];end
                    ind = ismember(L2_str.fROI.EVC.ids{sub}, aidx);
                    tmean(cnt) = nanmedian(L2_str.fROI.EVC.tvalues{sub}(ind));
                    
                case 2
                    % for EVC
                    aidx = []; for aR = 9:12, aidx = [aidx; L2_str.aROI.(aROIs{aR}){sub}];end
                    ind = ismember(L2_str.fROI.EVC.ids{sub}, aidx);
                    tmean(cnt) = nanmedian(L2_str.fROI.EVC.tvalues{sub}(ind));
                    
                case 3
                    % for LOC
                    aidx = []; for aR = 13:18, aidx = [aidx; L2_str.aROI.(aROIs{aR}){sub}];end
                    ind = ismember(L2_str.fROI.LOC.ids{sub}, aidx);
                    tmean(cnt) = nanmedian(L2_str.fROI.LOC.tvalues{sub}(ind));
                    
                case 4
                    % for VWFA
                    tmean(cnt) = nanmedian(L2_str.fROI.VWFA.tvalues{sub}); % FOR VWFA
                    
                case 5
                    % for WSCR
                    aidx = []; for aR = 19:22, aidx = [aidx; L2_str.aROI.(aROIs{aR}){sub}];end
                    ind = ismember(L2_str.fROI.WSCR.ids{sub}, aidx);
                    tmean(cnt) = nanmedian(L2_str.fROI.WSCR.tvalues{sub}(ind));
            end
        end
        
        FLUENCY(FLUENCY > 150) = NaN;
        FLUENCY(tmean < 0) = NaN;
        ISMAL(isnan(FLUENCY)) = NaN;
        [R(roi,rep) P(roi,rep)] = nancorrcoef([zscore(tmean(ISMAL==0)) zscore(tmean(ISMAL==1))],[zscore(FLUENCY(ISMAL == 0)) zscore(FLUENCY(ISMAL == 1))]);
    end
end
barweb(nanmean(R,2),nanstd(R,[],2))
%%
for i = 1:35
    tmean(i,1) = nanmedian(L2_str.fROI.VWFA.tvalues{i}(1:20)); % In the main figure, this constraint of 20 voxels is not considered.
end
fluency = fluency';
tmean(tmean < 0) = NaN; fluency(fluency > 150) = NaN;
Tsub = setdiff(find(ismal == 0),[find(isnan(fluency)); find(isnan(tmean))]);
Msub = setdiff(find(ismal == 1),[find(isnan(fluency)); find(isnan(tmean))]);

Tfluency = fluency(Tsub); Mfluency = fluency(Msub); zTfluency = zscore(Tfluency); zMfluency = zscore(Mfluency);
Ttmean = tmean(Tsub); Mtmean = tmean(Msub); zTtmean = zscore(Ttmean); zMtmean = zscore(Mtmean);
corrplot([zTtmean ;zMtmean], [zTfluency ; zMfluency]); ylim([-2 3])

