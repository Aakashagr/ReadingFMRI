% On real data
allclearL2
% Load data structure
% load L2fmri_READINGw

ismal = L2_str.ismal; % Is subject Malayalam: 1 (for Malayalam) and 0 (for Telugu)
qm = 1:34; qt = 35:68; % Index of Malayalam and Telugu stimuli
qs = 1:10; qd = 11:34; % Index of single and double letter stimuli

ids{1} = find(spm_read_vols(spm_vol('roidef_figure5/E_loc.nii'))); 
ids{2} = find(spm_read_vols(spm_vol('roidef_figure5/O_loc.nii'))); 
ids{3} = find(spm_read_vols(spm_vol('roidef_figure5/E_vwfa.nii')));
ids{4} = find(spm_read_vols(spm_vol('roidef_figure5/O_vwfa.nii')));

qsubm = find(ismal == 1);
qsubt = find(ismal == 0);
odd_sub = [qsubt(1:2:end); qsubm(1:2:end)];
even_sub = [qsubt(2:2:end); qsubm(2:2:end)];
%% Modelling activation of Bigrams using activations of single letters
for roi = 1%:2 % 3 = LOC
    % Selecting top N voxels to model and initialising the variables
    topnvoxels = 40; N = NaN(numel(ismal),topnvoxels);
    cvn = N; pvn = N; cvnn = N; pvnn = N; bn = NaN(numel(ismal),topnvoxels,3); bnn = bn;
    cpvn = N; ppvn = N; cpvnn = N; ppvnn = N;
    
    clear bpn bpnn
    for sub = 1:35
        if ismal(sub); N = qm; NN = qt; else, N = qt; NN = qm; end
        
        if roi ==1 && ismember(sub,odd_sub) 
            voxid = ids{1};
        elseif roi ==1 && ismember(sub,even_sub) 
            voxid = ids{2};
        elseif roi ==2 && ismember(sub,odd_sub) 
            voxid = ids{3};
        elseif roi ==2 && ismember(sub,even_sub) 
            voxid = ids{4};            
        end
        
        
        nvox = min(numel(voxid),topnvoxels);                % Total number of voxels considered in the analysis
        betas = L2_str.mergedevtbeta{sub}(voxid(1:nvox),:);
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

figure; barweb([nanmean(R_n,2) nanmean(R_nn,2)], [nansem(R_n,2) nansem(R_nn,2)]); ylim([0 0.5])
set(gca,'Xticklabel',{'LOC','VWFA'})
legend('Known', 'Unknown'); ylabel('Correlation Coefficient');

for i = 1:2
    P(i) = signrank(R_n(i,:),R_nn(i,:));
end

figure; barweb([nanmean(R_n(1,ismal==0),2) nanmean(R_nn(1,ismal==0),2); nanmean(R_n(1,ismal==1),2) nanmean(R_nn(1,ismal==1),2)], ...
               [nansem(R_n(1,ismal==0),2) nansem(R_nn(1,ismal==0),2); nansem(R_n(1,ismal==1),2) nansem(R_nn(1,ismal==1),2)]); ylim([0 0.5])
set(gca,'Xticklabel',{'Telugu readers','Malayalam readers'})
legend('Known', 'Unknown'); ylabel('Correlation Coefficient');
signrank(R_n(1,ismal == 1),R_nn(1,ismal == 1))


% figure; statcomparemean([bpn(:,1)- bpn(:,2)]./[bpn(:,1)+ bpn(:,2)], [bpnn(:,1)- bpnn(:,2)]./[bpnn(:,1)+ bpnn(:,2)]);
% figure; subplot(121); statcomparemean(bpn(:,:,1),bpn(:,:,2));
% subplot(122); statcomparemean(bpnn(:,:,1),bpnn(:,:,2));







