%% to look at granulometry, start by just getting an image in here
load('/Volumes/ryan_ims/sms_images/petro_im_instances/smg_7.mat');
petrim = smg_7;
green_im = im2double(imread(fullfile(petrim.main_path, petrim.im_subpaths{petrim.wavelengths == 530}, [petrim.sample_name, petrim.default_ext])));
im_mask = im2double(imread(fullfile(petrim.main_path, petrim.background_mask_path, [petrim.sample_name, petrim.default_ext])));

%% do some adaptive histogram adjustment
green_adj = adapthisteq(green_im,"NumTiles",[10 10]);
green_adj = imadjust(green_adj);

%% define an arbitrary set of radii, and use them to open the image
radii = linspace(10, 1000, 25);
% and we'll actually resize the image to simulate enlarging the radius at
% each step. 
resizes = radii(1)./radii;

% because ooids are dark here, need to invert the image to make them
% relatively light dots
invert_im = 1-green_adj;

% loop through and do stuff
opened_vol = zeros(round(size(invert_im,1)/10), round(size(invert_im,2)/10),numel(resizes));
for i = 1:numel(resizes)
    opened = imopen(imresize(invert_im.*im_mask, resizes(i)), strel("disk", radii(1)));
    
    opened_vol(:,:,i) = imresize(opened, [size(opened_vol, 1), size(opened_vol, 2)]);
end

% and we want to add a layer as the first element of the volume that
% represents opening with a radius of 0. 
opened = imopen(invert_im.*im_mask, strel("disk", 0));
opened_vol = cat(3,imresize(opened, [size(opened_vol, 1), size(opened_vol, 2)]), opened_vol);

%% The sort of detected grains are hidden in the differences between two consecutive openings
raw_grains = opened_vol(:,:,1:end-1) - opened_vol(:,:,2:end);

%% binarize the grain images with a threshold and do watershed

grain_binary = zeros(size(raw_grains));
watershed_grains = zeros(size(raw_grains));
for i = 1:size(raw_grains,3)
    grain_binary(:,:,i) = imbinarize(raw_grains(:,:,i));
    
    D = bwdist(~grain_binary(:,:,i));
    D = -D;
    D(~grain_binary(:,:,i)) = Inf;
        
    watershed_labels = watershed(D);
    watershed_labels(~grain_binary(:,:,i)) = 0;

    watershed_grains(:,:,i) = watershed_labels;
end

%% now we can do region props to get the centroids of all the grains

% get a cell array of all centroid positions
grain_centroids = cell(1, size(watershed_grains,3));
for i = 1:size(watershed_grains,3)
    these_cents = regionprops(watershed_grains(:,:,i),'centroid');
    extracted_centers = extractfield(these_cents,'Centroid');
    grain_centroids{i} = reshape(extracted_centers, 2,[])';
end

% to then select through grains, we want a mask that represents the at
% which radius each pixel takes its max
[~, max_mask] = max(raw_grains, [], 3);

% based on the cell array, construct grain images with perfectly circular
% grains at each centroid
[X,Y] = meshgrid(1:size(watershed_grains,2), 1:size(watershed_grains,1));
for i = 1:size(watershed_grains,3)
    these_cents = grain_centroids{i};
    
    for j = 1:size(these_cents,1)
        grain_mask = hypot(X - these_cents(j,1), Y - these_cents(j,2)) <= radii(i)/radii(1);
        if i ~= mode(max_mask(grain_mask))
            to_flip = find(watershed_grains(:,:,i) == j);
            [subinds1, subinds2] = ind2sub([size(watershed_grains,1), size(watershed_grains,2)], to_flip);
            watershed_grains(subinds1,subinds2,i) = 0;
        end
    end
end


%%
h = figure; 
ax = gca;
as.NextPlot = 'replaceChildren';

M(numel(resizes)-1) = struct('cdata', [], 'colormap', []);

h.Visible = 'off';

for i = 1:numel(resizes)-1
    c = imfuse(green_adj, imresize(opened_vol(:,:,i) - opened_vol(:,:,i+1), [size(green_adj, 1), size(green_adj, 2)]), "falsecolor");
    imshow(c)
    drawnow
    M(i) = getframe;
end

mean_rads = (radii(2:end) + radii(1:end-1))./2;

for i = 1:numel(M)
    imwrite(M(i).cdata, ['/Volumes/ryan_ims/sms_images/granulometry_experiment/' num2str(round(mean_rads(i))) '.tif'])
end

%% playing around a little with delta masks to get stream abundance vs. time
% first read in the datasets
data_10 = h5read('/Users/ryan/Downloads/agu_data.h5', '/manualChannelMasksqv10');
data_15 = h5read('/Users/ryan/Downloads/agu_data.h5', '/manualChannelMasksqv15');
data_30 = h5read('/Users/ryan/Downloads/agu_data.h5', '/manualChannelMasksqv30');

%% weird format, turn into 3d matrices
% brute force loop through to do it. Ugly, but works well enough

imcube_10 = zeros(size(data_10,1), size(data_10,2), size(data_10,3));
for i = 1:size(data_10, 1)
    for j = 1:size(data_10,2)
        for k = 1:size(data_10,3)
            if strcmp(data_10{i,j,k}, 'FALSE')
                imcube_10(i,j,k) = false;
            elseif strcmp(data_10{i,j,k}, 'TRUE')
                imcube_10(i,j,k) = true;
            end
        end
    end
end

imcube_15 = zeros(size(data_15,1), size(data_15,2), size(data_15,3));
for i = 1:size(data_15, 1)
    for j = 1:size(data_15,2)
        for k = 1:size(data_15,3)
            if strcmp(data_15{i,j,k}, 'FALSE')
                imcube_15(i,j,k) = false;
            elseif strcmp(data_15{i,j,k}, 'TRUE')
                imcube_15(i,j,k) = true;
            end
        end
    end
end

imcube_30 = zeros(size(data_30,1), size(data_30,2), size(data_30,3));
for i = 1:size(data_30, 1)
    for j = 1:size(data_30,2)
        for k = 1:size(data_30,3)
            if strcmp(data_30{i,j,k}, 'FALSE')
                imcube_30(i,j,k) = false;
            elseif strcmp(data_30{i,j,k}, 'TRUE')
                imcube_30(i,j,k) = true;
            end
        end
    end
end

%% get skeletonized binaries and number of branch points for each time slice

skel_10 = zeros(size(data_10,1), size(data_10,2), size(data_10,3));
bp_10 = zeros(size(skel_10, 3),1);
for i = 1:size(skel_10, 3)
    skel_10(:,:,i) = bwmorph(imcube_10(:,:,i), 'skel', Inf);
    bp_10(i) = sum(bwmorph(skel_10(:,:,i), 'branchpoints'), 'all');
end

skel_15 = zeros(size(data_15,1), size(data_15,2), size(data_15,3));
bp_15 = zeros(size(skel_15, 3),1);
for i = 1:size(skel_15, 3)
    skel_15(:,:,i) = bwmorph(imcube_15(:,:,i), 'skel', Inf);
    bp_15(i) = sum(bwmorph(skel_15(:,:,i), 'branchpoints'), 'all');
end

skel_30 = zeros(size(data_30,1), size(data_30,2), size(data_30,3));
bp_30 = zeros(size(skel_30, 3),1);
for i = 1:size(skel_30, 3)
    skel_30(:,:,i) = bwmorph(imcube_30(:,:,i), 'skel', Inf);
    bp_30(i) = sum(bwmorph(skel_30(:,:,i), 'branchpoints'), 'all');
end

%%
ids_10 = h5readatt('/Users/ryan/Downloads/agu_data.h5', '/manualChannelMasksqv10', 'IDs');
ids_15 = h5readatt('/Users/ryan/Downloads/agu_data.h5', '/manualChannelMasksqv15', 'IDs');
ids_30 = h5readatt('/Users/ryan/Downloads/agu_data.h5', '/manualChannelMasksqv30', 'IDs');

raw_meta1 = readtable('/Users/ryan/Downloads/tdb12_metadata.csv');
raw_meta2 = readtable('/Users/ryan/Downloads/tdb19_metadata.csv');

inds_10 = ismember(raw_meta1.linkID, ids_10); 
inds_15 = ismember(raw_meta2.linkID, ids_15); 
inds_30 = ismember(raw_meta2.linkID, ids_30); 

meta_10 = raw_meta1(inds_10,2:end);
meta_15 = raw_meta2(inds_15,2:end);
meta_30 = raw_meta2(inds_30,2:end);

bp_table_10 = array2table(bp_10, 'VariableNames', {'n_branch_points'});
bp_table_15 = array2table(bp_15, 'VariableNames', {'n_branch_points'});
bp_table_30 = array2table(bp_30, 'VariableNames', {'n_branch_points'});

final_table_10 = [bp_table_10, meta_10];
final_table_15 = [bp_table_15, meta_15];
final_table_30 = [bp_table_30, meta_30];

writetable(final_table_10, '/Users/ryan/Downloads/branchpoints_10.csv')
writetable(final_table_15, '/Users/ryan/Downloads/branchpoints_15.csv')
writetable(final_table_30, '/Users/ryan/Downloads/branchpoints_30.csv')