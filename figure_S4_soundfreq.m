allclear
load L2_sounds

ismalayalam = L2_str.subjinfo.ismalayalam;
images = L2_str.images;

dtel = zscore(L2_str.dissimilarity(:,ismalayalam == 0));
dmal = zscore(L2_str.dissimilarity(:,ismalayalam == 1));

%% Panel A
figure; corrplot(nanmean(dtel,2), nanmean(dmal,2),[],1); xlabel('Telugu subjects'); ylabel('Malayalam subjects')

%% saving data for reading experiment
sound_rm_mal = 9;  sound_rm_tel = 8; % Remove Uncommon sound pairs

id_mal = find(L2_str.sound_pairs(:,1) == sound_rm_mal | L2_str.sound_pairs(:,2) == sound_rm_mal);
id_tel = find(L2_str.sound_pairs(:,1) == sound_rm_tel | L2_str.sound_pairs(:,2) == sound_rm_tel);

sound_tel = mean(dtel,2); sound_tel(id_tel) = [];
sound_mal = mean(dmal,2); sound_mal(id_mal) = [];

%% panel B. Is sound error correlated with bigram frequency differences?
load('tel_bigram.mat'); 
load('mal_bigram.mat'); 

tel_bigram = mean(tel_bigram,2); tel_bigram = tel_bigram/sum(tel_bigram)*100;
mal_bigram = mean(mal_bigram,2); mal_bigram = mal_bigram/sum(mal_bigram)*100;

diff1 = tel_bigram - mal_bigram;
diff2 = sound_tel - sound_mal;

figure; 
xx = find(abs(diff1) > .335);
corrplot(diff1(xx),diff2(xx),'filtered pairs'); hold on;
plot(diff1(setdiff(1:630,xx)),diff2(setdiff(1:630,xx)),'.r')
xlabel('% Bigram differences, (Telugu - Malayalam)');
ylabel('Sound dissimilarity diff, (Telugu - Malayalam)');

%% Panel C and D
% Using perceived sound similarity
allclear
load('L2_letters.mat');
load('sound_distance.mat');

% % subjects to remove
L2_str.RT(:,[29 28 12],:) = [];
L2_str.subjinfo.ismalayalam([29 28 12]) = [];   
RT_mal = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:,:);
RT_tel = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:,:);

mRT_mal = nanmean(RT_mal,3);
mRT_tel = nanmean(RT_tel,3);

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2); 
RT_tel_mean = nanmean(mRT_tel,2);
mal = 1:630; tel = 631:1260;

visual_mal = 1./RT_mal_mean(mal); 
visual_tel = 1./RT_tel_mean(tel); 

figure; corrplot(sound_mal,visual_mal,'Malayalam Subjects'); 
ylabel('Visual dissimilarity'); xlabel('Audio dissimilarity')
figure; corrplot(sound_tel,visual_tel,'Telugu Subjects'); 
ylabel('Visual dissimilarity'); xlabel('Audio dissimilarity')

