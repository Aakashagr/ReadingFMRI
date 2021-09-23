
allclear
load ..\L2_telmal_bigram.mat;

L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.
mat = zeros(size(L2_str.RT)); mat(L2_str.RT>5) = 1; % identifying harder search pairs
temp = sum(sum(mat,2),3);
% figure; stem(temp) % identifying mistakes.
id = find(temp > 0 & temp < 2);

% removing higher RT outlier
for i = 1:length(id)
    [c, v] =ind2sub([size(L2_str.RT,2),2],find(mat(id(i),:,:)==1));
    L2_str.RT(id(i),c,v) = NaN;
end

RT = L2_str.RT; 
qtel = [1:300]; qmal = [301:600]; 

imgpairs = L2_str.img_pairs; 
ismal = L2_str.subjinfo.ismalayalam; 
imgpairs = nchoosek([1:25],2); 
allparts = L2_str.partsTD; 
q1 = find(allparts(:,2)==allparts(:,4)); % first letter different 
q2 = find(allparts(:,1)==allparts(:,3)); % second letter different 
% q1 = find(allparts(:,2)==allparts(:,4) & (allparts(:,1)~=allparts(:,2)) & (allparts(:,3)~=allparts(:,4))); % first letter different
% q2 = find(allparts(:,1)==allparts(:,3) & (allparts(:,1)~=allparts(:,2)) & (allparts(:,3)~=allparts(:,4))); % second letter different

%% Telugu bigrams, Telugu subjects
RT1 = RT(qtel(q1),find(~ismal),:);  
RT2 = RT(qtel(q2),find(~ismal),:); 

for letterid = 1:2
    for pid = 1:size(RT1,1)
        for subjid = 1:size(RT1,2)
            for repid = 1:size(RT1,3)
                if(letterid==1)
                    RTall(letterid,pid,subjid,repid) = RT1(pid,subjid,repid); 
                else
                    RTall(letterid,pid,subjid,repid) = RT2(pid,subjid,repid); 
                end
                pairall(letterid,pid,subjid,repid) = pid; 
                subjall(letterid,pid,subjid,repid) = subjid; 
                letterall(letterid,pid,subjid,repid) = letterid; 
            end
        end
    end
end

panovatt = anovan(vec(RTall),{vec(letterall),vec(pairall),vec(subjall)},'model','full','display','off'); 

figure; 
subplot(221); 
data = [nanmean(nanmean(RT1,3),2) nanmean(nanmean(RT2,3),2)]; 
bar(nanmean(data,1)); hold on; errorbar([1:2],nanmean(data,1),nansem(data,1),'.'); 
title(sprintf('TBTS, pletter = %2.2g',panovatt(1)));  
ylim([0 3.5]); 

%% Telugu bigrams, Mallu subjects
RT1 = RT(qtel(q1),find(ismal),:);  
RT2 = RT(qtel(q2),find(ismal),:); 

for letterid = 1:2
    for pid = 1:size(RT1,1)
        for subjid = 1:size(RT1,2)
            for repid = 1:size(RT1,3)
                if(letterid==1)
                    RTall(letterid,pid,subjid,repid) = RT1(pid,subjid,repid); 
                else
                    RTall(letterid,pid,subjid,repid) = RT2(pid,subjid,repid); 
                end
                pairall(letterid,pid,subjid,repid) = pid; 
                subjall(letterid,pid,subjid,repid) = subjid; 
                letterall(letterid,pid,subjid,repid) = letterid; 
            end
        end
    end
end

panovatm = anovan(vec(RTall),{vec(letterall),vec(pairall),vec(subjall)},'model','full','display','off');
subplot(222); 
data = [nanmean(nanmean(RT1,3),2) nanmean(nanmean(RT2,3),2)]; 
bar(nanmean(data,1)); hold on; errorbar([1:2],nanmean(data,1),nansem(data,1),'.'); 
title(sprintf('TBMS, pletter = %2.2g',panovatm(1)));  
ylim([0 3.5]); 
%% Mallu bigrams, Mallu subjects
RT1 = RT(qmal(q1),find(ismal),:);  
RT2 = RT(qmal(q2),find(ismal),:); 

for letterid = 1:2
    for pid = 1:size(RT1,1)
        for subjid = 1:size(RT1,2)
            for repid = 1:size(RT1,3)
                if(letterid==1)
                    RTall(letterid,pid,subjid,repid) = RT1(pid,subjid,repid); 
                else
                    RTall(letterid,pid,subjid,repid) = RT2(pid,subjid,repid); 
                end
                pairall(letterid,pid,subjid,repid) = pid; 
                subjall(letterid,pid,subjid,repid) = subjid; 
                letterall(letterid,pid,subjid,repid) = letterid; 
            end
        end
    end
end

panovamm = anovan(vec(RTall),{vec(letterall),vec(pairall),vec(subjall)},'model','full','display','off');
subplot(223); 
data = [nanmean(nanmean(RT1,3),2) nanmean(nanmean(RT2,3),2)]; 
bar(nanmean(data,1)); hold on; errorbar([1:2],nanmean(data,1),nansem(data,1),'.'); 
title(sprintf('MBMS, pletter = %2.2g',panovamm(1)));  
ylim([0 3.5]); 
%% Mallu bigrams, Telugu subjects
RT1 = RT(qmal(q1),find(~ismal),:); 
RT2 = RT(qmal(q2),find(~ismal),:); 

for letterid = 1:2
    for pid = 1:size(RT1,1)
        for subjid = 1:size(RT1,2)
            for repid = 1:size(RT1,3)
                if(letterid==1)
                    RTall(letterid,pid,subjid,repid) = RT1(pid,subjid,repid); 
                else
                    RTall(letterid,pid,subjid,repid) = RT2(pid,subjid,repid); 
                end
                pairall(letterid,pid,subjid,repid) = pid; 
                subjall(letterid,pid,subjid,repid) = subjid; 
                letterall(letterid,pid,subjid,repid) = letterid; 
            end
        end
    end
end

panovamt = anovan(vec(RTall),{vec(letterall),vec(pairall),vec(subjall)},'model','full','display','off');

subplot(224); 
data = [nanmean(nanmean(RT1,3),2) nanmean(nanmean(RT2,3),2)]; 
bar(nanmean(data,1)); hold on; errorbar([1:2],nanmean(data,1),nansem(data,1),'.'); 
title(sprintf('MBTS, pletter = %2.2g',panovamt(1)));  
ylim([0 3.5]); 