%% Figure 1: All the letter stimuli shown in the expt.
allclear
load L2_letters
imgs = L2_str.images;

% Changing image contrast
for q = 1:length(imgs)
    imgs{q} = imcomplement(imgs{q});
end

manyimagesc(imgs,L2_str.sound,4,18);

