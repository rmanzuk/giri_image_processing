% wrapper to go through and assess textural metrics. Various applications
% of functions. Will need to be cleaned up once we know what is effective.
%
% R. A. Manzuk 10/20/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set everything up
% path to rgb images for identifying patches + directory
rgb_path = '/Users/ryan/Dropbox (Princeton)/reef_survey_project/nevada_qgis/rgb_jpgs';
rgb_dir = dir(fullfile(rgb_path,'*.jpg'));

% path to all of the full res tifs
path_470 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/470nm'; 
path_505 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/505nm'; 
path_530 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/530nm'; 
path_590 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/590nm'; 
path_625 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/625nm'; 
path_760 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/760nm'; 
path_940 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/940nm'; 
path_365 = '/Volumes/Extreme SSD/sms_images/original_tifs/365nm'; 

% in my dataset, we'll be clicking the patches on low-res jpegs for
% efficiency, but doing analyses on full-res tifs. So we need to know their
% scaling ratio
jpg_test = imread(fullfile(rgb_path,rgb_dir(1).name));
test_dir_470 = dir(fullfile(path_470, '*.tif'));
tif_test = im2double(imread(fullfile(path_470,test_dir_470(10).name)));
scale_ratio = size(jpg_test,1)/size(tif_test,1);
%% select patches from the rgb images
patch_types = {'whole rock', 'dolomitized sediments', 'nondolomitized sediments', 'archaeocyath', 'ooids'};

patch_table_path = '/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_texture.xlsx';
% use line below to save a new patch table
%[patches_table] = patch_selector(rgb_path,'*.jpg',patch_types,'savePath', patch_table_path,'scaleRatio',scale_ratio);
% use the line below if patch table already exists and want to add to it
[patches_table] = patch_selector(rgb_path,'*.jpg',patch_types,'savePath', patch_table_path,'scaleRatio',scale_ratio, 'existingTable',readtable(patch_table_path));


%% okay, we're ready to go through and grab a bunch of textural information
% make sure the spreadsheets are in here
patches_table = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_texture.xlsx');
image_key = readtable('/Volumes/Extreme SSD/sms_images/image_logs/masterlog.xlsx');

% and we'll go through each image once, so we need a uniqe set of names
% from the patches table
unique_names = unique(patches_table.image_name);

% make up a table for basic statistics
% start with the stats we'll collect for each channel
stats = {'glcm_contrast', 'glcm_correlation', 'glcm_energy' 'glcm_homogeneity', 'mean', 'variance', 'rms', 'entropy'};
wavelengths = cellstr(string(num2str(unique(image_key.im_wavelength))));

% loop through and give a heading for each possible stat with wavelength
stat_headings = cell(1,numel(stats)*numel(wavelengths));
for i = 1:numel(wavelengths)
    for j = 1:numel(stats)
        stat_headings{(i-1)*numel(stats) + j} = [wavelengths{i} '_' stats{j}];
    end
end

% make an empty table with the image names, patch names, and stat headings.
% To update as we go
im_stuff_headings = patches_table.Properties.VariableNames;
final_stat_table = cell2table(cell(0,numel(stat_headings) + numel(im_stuff_headings)),'VariableNames',[im_stuff_headings, stat_headings]);


% we need to know the size of a downsampled image of the same resolution as
% the one we use for textons...for later
resize_test = imresize(tif_test,0.1);

% set up a wait bar
f =  waitbar(0, 'running texture analysis');
% now we can loop through each image and get some info
for i = 1:numel(unique_names)
    % find the patches in the table
    these_patches = patches_table(strcmp(patches_table.image_name,unique_names{i}),:);
    % and we need the image names
    these_images = image_key(strcmp(image_key.actual_sample,unique_names{i}),:);
    % list out the wavelenths we'll use
    im_wavelengths = these_images.im_wavelength;
    % we can make up a table that will hold the stats for this image
    this_stat_table = cell2table(cell(0,numel(stat_headings)),'VariableNames',stat_headings);
    % and fill it with nans for now
    this_stat_table{1:size(these_patches,1),:} = missing;

    % now sub loop through all of the wavelengths and gather data from the
    % images
    for j = 1:numel(im_wavelengths)
        % if statements to choose directory
        if im_wavelengths(j) == 365
            search_folder = path_365;
        elseif im_wavelengths(j) == 470
            search_folder = path_470;
        elseif im_wavelengths(j) == 505
            search_folder = path_505;
        elseif im_wavelengths(j) == 530
            search_folder = path_530;
        elseif im_wavelengths(j) == 590
            search_folder = path_590;
        elseif im_wavelengths(j) == 625
            search_folder = path_625;
        elseif im_wavelengths(j) == 760
            search_folder = path_760;
        elseif im_wavelengths(j) == 940
            search_folder = path_940;
        end
        % so now we can read the image
        [~,image_name] = fileparts(these_images.file_name(j));
        this_image = im2double(imread(fullfile(search_folder,[image_name '.tif'])));

        % one last loop to nest to go into the patches :/
        for k = 1:size(these_patches,1)
            % index the patch
            % but first need to check the patch indices to see if they make
            % sense
            if these_patches.top_row(k) < 1
                these_patches.top_row(k) = 1;
            end
            if these_patches.bot_row(k) > size(this_image,1)
                these_patches.bot_row(k) = size(this_image,1);
            end
            if these_patches.left_col(k) < 1
                these_patches.left_col(k) = 1;
            end
            if these_patches.right_col(k) > size(this_image,2)
                these_patches.right_col(k) = size(this_image,2);
            end
  
            this_patch = this_image(these_patches.top_row(k):these_patches.bot_row(k),these_patches.left_col(k):these_patches.right_col(k));

            % and now we're clear to set grab the data. 
            % first to glcm stuff
            glcms = graycomatrix(this_patch);
            glcm_props = graycoprops(glcms);

            % calculate local binary pattern features
            lbp_features = extractLBPFeatures(this_patch);

            % take the mean, variance, and some other stuff
            im_mean = mean(this_patch(:));
            im_var = var(this_patch(:));
            im_rms = rms(this_patch(:));
            im_ent = entropy(this_patch);

            % phew, we did it. Now we need to get everything into tables.
            % first handle the basic stats, figure out which columns we are
            % in based upon wavelength
            these_columns =  find(contains(stat_headings,num2str(im_wavelengths(j))));
            this_stat_table(k,these_columns(1)) = {glcm_props.Contrast};
            this_stat_table(k,these_columns(2)) = {glcm_props.Correlation};
            this_stat_table(k,these_columns(3)) = {glcm_props.Energy};
            this_stat_table(k,these_columns(4)) = {glcm_props.Homogeneity};
            this_stat_table(k,these_columns(5)) = {im_mean};
            this_stat_table(k,these_columns(6)) = {im_var};
            this_stat_table(k,these_columns(7)) = {im_rms};
            this_stat_table(k,these_columns(8)) = {im_ent};    
        end
        
    end
    % now that we've gotten through all of the patches and wavelenghts, put
    % the patches table
    image_stat_table = [these_patches, this_stat_table];
    final_stat_table = [final_stat_table; image_stat_table];

    waitbar(i/numel(unique_names), f, 'running texture analysis');
end



%% set up filter bank(s) for textons
% we'll scale the images for textons. So define that scaling
scaling = 0.15;

% first textural thing we'll do is texton feature mapping
% start with Leung-Malik filter bank from online source: https://www.robots.ox.ac.uk/~vgg/research/texclass/filters.html
lm_small_filter_bank = makeLMfilters(18);
% and make a larger version too, twice the size
lm_large_filter_bank = imresize(lm_small_filter_bank,[99,99]);

% also a schmid filter bank from the same source as above
schmid_filter_bank = makeSfilters;

%% get the texton library by looking at the whole rock patch of each image

% make sure the spreadsheets are in here
patches_table = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_texture.xlsx');
image_key = readtable('/Volumes/Extreme SSD/sms_images/image_logs/masterlog.xlsx');

% and we'll go through each image once, so we need a uniqe set of names
% from the patches table
unique_names = unique(patches_table.image_name);

% we need to know the size of a downsampled image of the same resolution as
% the one we use for textons...for later
resize_test = imresize(tif_test,scaling);

% set up empty vectors to receive the filter responses
lm_small_responses = [];
lm_large_responses = [];
schmid_responses = [];

% set up a wait bar
f =  waitbar(0, 'building texton library');
% now we can loop through each image and get some info
for i = 1:numel(unique_names)

    % find the patches in the table
    these_patches = patches_table(strcmp(patches_table.image_name,unique_names{i}),:);
    % and we need the image names

    these_images = image_key(strcmp(image_key.actual_sample,unique_names{i}),:);

    % list out the wavelenths we have
    im_wavelengths = these_images.im_wavelength;

    % we're only gonna use the green whole rock patch for now
    % so index the green
    green_ind = im_wavelengths == 530;
    
    % so now we can read the image
    [~,image_name] = fileparts(these_images.file_name(green_ind));
    this_image = im2double(imread(fullfile(path_530,[image_name '.tif'])));

    % figure out which patch is for whole rock
    whole_rock_ind = strcmp(these_patches.patch_type,'whole rock');

    % but first need to check the patch indices to see if they make
    % sense
    % and only try if there was a whole rock patch
    if sum(whole_rock_ind) > 0 
        if these_patches.top_row(whole_rock_ind) < 1
            these_patches.top_row(whole_rock_ind) = 1;
        end
        if these_patches.bot_row(whole_rock_ind) > size(this_image,1)
            these_patches.bot_row(whole_rock_ind) = size(this_image,1);
        end
        if these_patches.left_col(whole_rock_ind) < 1
            these_patches.left_col(whole_rock_ind) = 1;
        end
        if these_patches.right_col(whole_rock_ind) > size(this_image,2)
            these_patches.right_col(whole_rock_ind) = size(this_image,2);
        end
        % now we're ready to index the patch
        this_patch = this_image(these_patches.top_row(whole_rock_ind):these_patches.bot_row(whole_rock_ind),these_patches.left_col(whole_rock_ind):these_patches.right_col(whole_rock_ind));

        % time to do filter bank texton stuff
        % normalize the patch
        patch_normalized = normalize_image(this_patch,'zeroCenter');

        % do it with the three different filter banks
        lm_small_responses = [lm_small_responses; texton_jemwa(patch_normalized,scaling,lm_small_filter_bank,100)];
        lm_large_responses = [lm_large_responses; texton_jemwa(patch_normalized,scaling,lm_large_filter_bank,100)];
        schmid_responses = [schmid_responses; texton_jemwa(patch_normalized,scaling,schmid_filter_bank,100)];

        % update the waitbar
        waitbar(i/numel(unique_names), f, 'building texton library')
    end
end

%% need to kmeans the texton libraries.
[~,lm_small_centroids] = kmeans(lm_small_responses,30,'MaxIter',1000);
[~,lm_large_centroids] = kmeans(lm_large_responses,30,'MaxIter',1000);
[~,schmid_centroids] = kmeans(schmid_responses,30,'MaxIter',1000);

%% and now we have to go through and use the library to assign each pixel in the dataset to a texton feature

% make sure the spreadsheets are in here
patches_table = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_texture.xlsx');
image_key = readtable('/Volumes/Extreme SSD/sms_images/image_logs/masterlog.xlsx');

% and we'll go through each image once, so we need a uniqe set of names
% from the patches table
unique_names = unique(patches_table.image_name);

% we need to know the size of a downsampled image of the same resolution as
% the one we use for textons...for later
resize_test = imresize(tif_test,scaling);

% set up empty vectors to receive the total number of assignments
lm_small_counts = [];
lm_large_counts = [];
schmid_counts = [];

% set up a wait bar
f =  waitbar(0, 'assigning textons');
% now we can loop through each image and get some info
for i = 1:numel(unique_names)

    % find the patches in the table
    these_patches = patches_table(strcmp(patches_table.image_name,unique_names{i}),:);
    % and we need the image names

    these_images = image_key(strcmp(image_key.actual_sample,unique_names{i}),:);

    % list out the wavelenths we have
    im_wavelengths = these_images.im_wavelength;

    % we're only gonna use the green whole rock patch for now
    % so index the green
    green_ind = im_wavelengths == 530;
    
    % so now we can read the image
    [~,image_name] = fileparts(these_images.file_name(green_ind));
    this_image = im2double(imread(fullfile(path_530,[image_name '.tif'])));

    % figure out which patch is for whole rock
    whole_rock_ind = strcmp(these_patches.patch_type,'whole rock');

    % but first need to check the patch indices to see if they make
    % sense
    % and only try if there was a whole rock patch
    if sum(whole_rock_ind) > 0 
        if these_patches.top_row(whole_rock_ind) < 1
            these_patches.top_row(whole_rock_ind) = 1;
        end
        if these_patches.bot_row(whole_rock_ind) > size(this_image,1)
            these_patches.bot_row(whole_rock_ind) = size(this_image,1);
        end
        if these_patches.left_col(whole_rock_ind) < 1
            these_patches.left_col(whole_rock_ind) = 1;
        end
        if these_patches.right_col(whole_rock_ind) > size(this_image,2)
            these_patches.right_col(whole_rock_ind) = size(this_image,2);
        end
        % now we're ready to index the patch
        this_patch = this_image(these_patches.top_row(whole_rock_ind):these_patches.bot_row(whole_rock_ind),these_patches.left_col(whole_rock_ind):these_patches.right_col(whole_rock_ind));

        % time to do filter bank texton stuff
        % normalize the patch
        patch_normalized = normalize_image(this_patch,'zeroCenter');

        % assign the textons for each of the filter banks
        texton_image_lm_small = texton_apply_clusters(patch_normalized, scaling, lm_small_filter_bank, lm_small_centroids);
        texton_image_lm_large = texton_apply_clusters(patch_normalized, scaling, lm_large_filter_bank, lm_large_centroids);
        texton_image_schmid = texton_apply_clusters(patch_normalized, scaling, schmid_filter_bank, schmid_centroids);
        
        % and write the texton images as outputs
        imwrite(uint8(texton_image_lm_small),fullfile('/Volumes/Extreme SSD/sms_images/texton_ims/lm_small',[image_name '.tif']));
        imwrite(uint8(texton_image_lm_large),fullfile('/Volumes/Extreme SSD/sms_images/texton_ims/lm_large',[image_name '.tif']));
        imwrite(uint8(texton_image_schmid),fullfile('/Volumes/Extreme SSD/sms_images/texton_ims/schmid',[image_name '.tif']));
    
        % and put the counts in the waiting matrices
        lm_small_counts = [lm_small_counts; texton_counts(texton_image_lm_small,30)'];
        lm_large_counts = [lm_large_counts; texton_counts(texton_image_lm_large,30)'];
        schmid_counts = [schmid_counts; texton_counts(texton_image_schmid,30)'];

        % update the waitbar
        waitbar(i/numel(unique_names), f, 'assigning textons')
    end
end

close(f);

%% from the saved texton images we can extract some statistics as well
% make sure the spreadsheets are in here
patches_table = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_texture.xlsx');
image_key = readtable('/Volumes/Extreme SSD/sms_images/image_logs/masterlog.xlsx');

% and we'll go through each image once, so we need a uniqe set of names
% from the patches table
unique_names = unique(patches_table.image_name);

% id the folders with all of the texton images
lm_small_dir = '/Volumes/Extreme SSD/sms_images/texton_ims/lm_small';
lm_large_dir = '/Volumes/Extreme SSD/sms_images/texton_ims/lm_large';
schmid_dir = '/Volumes/Extreme SSD/sms_images/texton_ims/schmid';

% set up a bunch of empty stuff
all_lm_small_counts = [];
all_lm_large_counts = [];
all_schmid_counts = [];

all_lm_small_area = {};
all_lm_large_area = {};
all_schmid_area = {};

all_lm_small_orientation = {};
all_lm_large_orientation = {};
all_schmid_area_orientation = {};

all_lm_small_eccentricity = {};
all_lm_large_eccentricity = {};
all_schmid_area_eccentricity = {};

all_lm_small_lacuniarity = {};
all_lm_large_lacuniarity = {};
all_schmid_area_lacuniarity = {};

max_box_sizes = [];

% set up a wait bar
f =  waitbar(0, 'getting texton stats');
% now we can loop through each image and get some info
for i = 1:numel(unique_names)

    % we need to go from the sample name to the image name
    these_images = image_key(strcmp(image_key.actual_sample,unique_names{i}),:);

    % list out the wavelenths we have
    im_wavelengths = these_images.im_wavelength;

    % we're only gonna use the green whole rock patch for now
    % so index the green
    green_ind = im_wavelengths == 530;
    
    % so now we can read the texton images
    [~,image_name] = fileparts(these_images.file_name(green_ind));
    lm_small_image = imread(fullfile(lm_small_dir,[image_name '.tif']));
    lm_large_image = imread(fullfile(lm_large_dir,[image_name '.tif']));
    schmid_image = imread(fullfile(schmid_dir,[image_name '.tif']));

    % and get the stats for all of them
    [all_lm_small_counts(:,i), lm_small_area, lm_small_orientation, lm_small_eccentricity, lm_small_lacuniarity,max_box_sizes(i)] = texton_patch_stats(lm_small_image,30);
    [all_lm_large_counts(:,i), lm_large_area, lm_large_orientation, lm_large_eccentricity, lm_large_lacuniarity] = texton_patch_stats(lm_large_image,30);
    [all_schmid_counts(:,i), schmid_area, schmid_orientation, schmid_eccentricity, schmid_lacuniarity] = texton_patch_stats(schmid_image,30);

    % set everything for saving
    all_lm_small_area{i} = lm_small_area;
    all_lm_large_area{i} = lm_large_area;
    all_schmid_area{i} = schmid_area;

    all_lm_small_orientation{i} = lm_small_orientation;
    all_lm_large_orientation{i} = lm_large_orientation;
    all_schmid_area_orientation{i} = schmid_orientation;

    all_lm_small_eccentricity{i} = lm_small_eccentricity;
    all_lm_large_eccentricity{i} = lm_large_eccentricity;
    all_schmid_area_eccentricity{i} = schmid_eccentricity;

    all_lm_small_lacuniarity{i} = lm_small_lacuniarity;
    all_lm_large_lacuniarity{i} = lm_large_lacuniarity;
    all_schmid_area_lacuniarity{i} = schmid_lacuniarity;

    % update the waitbar
    waitbar(i/numel(unique_names), f, 'getting texton stats')

end

close(f);

%% experiments with SMG 60 to detect skeletal material

% set the scaling
scaling = 0.25;

% read in the indices
shell_inds = imread('/Volumes/Extreme SSD/sms_images/masks/smg_60_shell.tif');
 
% resize the mask and take out the shell/not shell masks
inds_resized = imresize(shell_inds,scaling,"nearest");
shell_mask = inds_resized == 4;
not_shell_mask = inds_resized == 3;
bound_mask = inds_resized == 2;
bound_mask = ~bound_mask;

% and make a filter bank that's the schmid filter bank and the blobby parts
% of the lm filter bank
%blobby_filter_bank = cat(3,lm_small_filter_bank(:,:,[37:end]),schmid_filter_bank);

% read in and normalize images
blue_im = im2double(imread('/Volumes/Extreme SSD/sms_images/color_corrected_tifs/470nm/cap89070_60001.tif'));
cyan_im = im2double(imread('/Volumes/Extreme SSD/sms_images/color_corrected_tifs/505nm/cap89127_60002.tif'));
green_im = im2double(imread('/Volumes/Extreme SSD/sms_images/color_corrected_tifs/530nm/cap89184_60003.tif'));
yellow_im = im2double(imread('/Volumes/Extreme SSD/sms_images/color_corrected_tifs/590nm/cap89241_60004.tif'));
red_im = im2double(imread('/Volumes/Extreme SSD/sms_images/color_corrected_tifs/625nm/cap89298_60005.tif'));
rededge_im = im2double(imread('/Volumes/Extreme SSD/sms_images/color_corrected_tifs/760nm/cap89355_60006.tif'));
ir_im = im2double(imread('/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/760nm/cap89355_60006.tif'));
uv_im = im2double(imread('/Volumes/Extreme SSD/sms_images/original_tifs/365nm/cap89469_60008.tif'));

blue_normalized = normalize_image(blue_im,'zeroCenter');
cyan_normalized = normalize_image(cyan_im,'zeroCenter');
green_normalized = normalize_image(green_im,'zeroCenter');
yellow_normalized = normalize_image(yellow_im,'zeroCenter');
red_normalized = normalize_image(red_im,'zeroCenter');
rededge_normalized = normalize_image(rededge_im,'zeroCenter');
ir_normalized = normalize_image(ir_im,'zeroCenter');
uv_im_normalized = normalize_image(uv_im,'zeroCenter');

% make the filter responses on the yellow image
[responses_lm_small_yellow] = texton_jemwa(yellow_normalized, scaling, lm_small_filter_bank);
[responses_lm_large_yellow] = texton_jemwa(yellow_normalized, scaling, lm_large_filter_bank);
[responses_schmid_yellow] = texton_jemwa(yellow_normalized, scaling, schmid_filter_bank);
[responses_lm_small_blue] = texton_jemwa(blue_normalized, scaling, lm_small_filter_bank);
[responses_lm_large_blue] = texton_jemwa(blue_normalized, scaling, lm_large_filter_bank);
[responses_schmid_blue] = texton_jemwa(blue_normalized, scaling, schmid_filter_bank);

% now we can set up a column matrix of filter responses and color values
all_dataset = [reshape(responses_lm_small_yellow,[],size(lm_small_filter_bank,3)), ...
    reshape(responses_lm_large_yellow,[],size(lm_large_filter_bank,3)), ...
    reshape(responses_schmid_yellow,[],size(schmid_filter_bank,3)), ...
    reshape(responses_lm_small_blue,[],size(lm_small_filter_bank,3)), ...
    reshape(responses_lm_large_blue,[],size(lm_large_filter_bank,3)), ...
    reshape(responses_schmid_blue,[],size(schmid_filter_bank,3)), ...
    reshape(imresize(blue_im,scaling),[],1), reshape(imresize(cyan_im,scaling),[],1),...
    reshape(imresize(green_im,scaling),[],1), reshape(imresize(yellow_im,scaling),[],1),...
    reshape(imresize(red_im,scaling),[],1), reshape(imresize(rededge_im,scaling),[],1),...
    reshape(imresize(ir_im,scaling),[],1), reshape(imresize(uv_im,scaling),[],1)];

% turn the masks into columns for indexing
shell_col = reshape(shell_mask,[],1);
not_shell_col = reshape(not_shell_mask,[],1);
bound_col = reshape(bound_mask,[],1);

% get a random sampler for the not_shell training data to have equal
% numbers for each category
sampler = randsample(1:sum(not_shell_col),sum(shell_col));

% grab the masked portions of the dataset
all_not_shell = all_dataset(not_shell_col,:);
all_shell = all_dataset(shell_col,:);
all_masked = all_dataset(bound_col,:);

% index as training data
training_dataset = [all_shell;all_not_shell(sampler,:)];
training_dataset_inds = [ones(sum(shell_col),1);2*ones(sum(shell_col),1)];

% and train the model
shell_model = fitcsvm(training_dataset,training_dataset_inds);

% predict the rest of teh image
all_classifications = predict(shell_model,all_dataset);

% reshape into an image
shell_segmentation = reshape(all_classifications, size(responses_lm_small,1), size(responses_lm_small,2)).*bound_mask;


%% dig a bit deeper into the color/texture dataset with a pca
% first normalize with zero center, standard deviation of 1
all_data_normalized = normalize(all_masked);

% and run a pca on the normalized data
[coeff,score,latent] = pca(all_data_normalized);

% do tsne too
%tsne_data = tsne(all_data_normalized);

% adjust the shell/not shell column masks for teh boundary mask
shell_w_boundary = shell_col(bound_col);
not_shell_w_boundary = not_shell_col(bound_col);

%% plot some pca stuff

% get the percent explained by each variable
pct_explained = latent./sum(latent).*100;

% get the top 10 contributers to the first five principal components
[~,largest_inds] = maxk(abs(coeff(:,1:5)),10);

% make bar plots of loadings for those max contributers
figure()
colormap(flipud(cmocean('rain')))
for i = 1:5
    subplot(2,3,i)
    bar(coeff(largest_inds(:,i),i))
    these_names = num2cell(largest_inds(:,i));
    set(gca,'xticklabel',these_names)
    title(['PC ' num2str(i)])
end

% join the filter banks into 1
joint_filter_bank = cat(3,lm_small_filter_bank,lm_small_filter_bank,schmid_filter_bank);

% make figures for the pcs most important filters
figure()
colormap(flipud(cmocean('rain')))
tcl1 = tiledlayout(5,2);
for i = 1:10
    nexttile
    imagesc(joint_filter_bank(:,:,largest_inds(i,1)))
    axis image
end
title(tcl1,'PC 1')

figure()
colormap(flipud(cmocean('rain')))
tcl1 = tiledlayout(5,2);
for i = 1:10
    nexttile
    imagesc(joint_filter_bank(:,:,largest_inds(i,2)))
    axis image
end
title(tcl1,'PC 2')

figure()
colormap(flipud(cmocean('rain')))
tcl1 = tiledlayout(5,2);
for i = 1:10
    nexttile
    imagesc(joint_filter_bank(:,:,largest_inds(i,3)))
    axis image
end
title(tcl1,'PC 3')

figure()
colormap(flipud(cmocean('rain')))
tcl1 = tiledlayout(5,2);
for i = 1:10
    nexttile
    imagesc(joint_filter_bank(:,:,largest_inds(i,4)))
    axis image
end
title(tcl1,'PC 4')

figure()
colormap(flipud(cmocean('rain')))
tcl1 = tiledlayout(5,2);
for i = 1:10
    nexttile
    imagesc(joint_filter_bank(:,:,largest_inds(i,5)))
    axis image
end
title(tcl1,'PC 5')

%% ichnofabric index images to develop anisotropy metric

% make the directory, images from bottjer and droser 1991
ichnofab_path = ('/Volumes/Extreme SSD/sms_images/icnofabric_inds');

% make the filter bank. Let's use a combo schmid and lm filter bank, where
% the lm filter bank has a 18 orientations (separated by 10 degrees)
n_orientations = 18;
lm_filter_bank = makeLMfilters(n_orientations);
schmid_filter_bank = makeSfilters;
all_filter_bank = cat(3,lm_filter_bank,schmid_filter_bank);
n_filters = size(all_filter_bank,3);

% we want to mark down which filters in the bank are oriented, and at which
% angle
orient_inds = [true(n_orientations*6,1); false(n_filters-n_orientations*6,1)];
orient_angles = repmat(linspace(0,180-180/n_orientations,n_orientations)',[6,1]);

% function to verify anisotropy metric
%calibrate_anisotropy(ichnofab_path, '*.tif', all_filter_bank, orient_inds);

%% apply anisotropy metric to image dataset
% make sure the spreadsheets are in here
patches_table = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_texture.xlsx');
image_key = readtable('/Volumes/Extreme SSD/sms_images/image_logs/masterlog.xlsx');

% and we'll go through each image once, so we need a uniqe set of names
% from the patches table
unique_names = unique(patches_table.image_name);

% make an empty table with the image names, patch names, and stat headings.
% To update as we go
im_stuff_headings = patches_table.Properties.VariableNames;
stat_headings = {'pct_var_explained', 'filt_1_ind', 'filt_2_ind', 'filt_3_ind', 'filt_4_ind', 'filt_5_ind', 'filt_6_ind'};
final_stat_table = cell2table(cell(0,numel(stat_headings) + numel(im_stuff_headings)),'VariableNames',[im_stuff_headings, stat_headings]);

% set our scaling ratio for filtering
scaling = 0.15;

% set up a wait bar
f =  waitbar(0, 'running anisotropy analysis');
% now we can loop through each image and get some info
for i = 1:numel(unique_names)
    % find the patches in the table
    these_patches = patches_table(strcmp(patches_table.image_name,unique_names{i}),:);
    % and we need the image names
    these_images = image_key(strcmp(image_key.actual_sample,unique_names{i}),:);
    % list out the wavelenths we'll use
    im_wavelengths = these_images.im_wavelength;
    % we can make up a table that will hold the stats for this image
    this_stat_table = cell2table(cell(0,numel(stat_headings)),'VariableNames',stat_headings);
    % and fill it with nans for now
    this_stat_table{1:size(these_patches,1),:} = missing;

    % for this, we'll use a brightness image, so we need to read in red,
    % green and blue
    red_ind = im_wavelengths == 625;
    green_ind = im_wavelengths == 530;
    blue_ind = im_wavelengths == 470;

    [~,red_name] = fileparts(these_images.file_name(red_ind));
    [~,green_name] = fileparts(these_images.file_name(green_ind));
    [~,blue_name] = fileparts(these_images.file_name(blue_ind));

    red_image = im2double(imread(fullfile(path_625,[red_name '.tif'])));
    green_image = im2double(imread(fullfile(path_530,[green_name '.tif'])));
    blue_image = im2double(imread(fullfile(path_470,[blue_name '.tif'])));
    % assemble the grayscale image
    grayscale_image = rgb2gray(cat(3,red_image,green_image,blue_image));

    % we're going to get the filter responses for this whole image before
    % we bother diving into the patches
    % first window normalize. 
    median_filted = medfilt2(imresize(grayscale_image,scaling,'bicubic'));
    smoothed = imgaussfilt(median_filted,3);
    for_filtering = normalize_image(smoothed, 'windows', 7);

    % get the filter_responses
    normalized_responses = texton_jemwa(for_filtering, 1, gabor_filter_bank);
    
    % now loop into the patches
    for j = 1:size(these_patches,1)
         % index the patch
        % but first need to check the patch indices to see if they make
        % sense
        if these_patches.top_row(j) < 1
            these_patches.top_row(j) = 1;
        end
        if these_patches.bot_row(j)*scaling > size(normalized_responses,1)
            these_patches.bot_row(j) = size(normalized_responses,1)/scaling;
        end
        if these_patches.left_col(j) < 1
            these_patches.left_col(j) = 1;
        end
        if these_patches.right_col(j)*scaling > size(normalized_responses,2)
            these_patches.right_col(j) = size(normalized_responses,2);
        end
        
        % set the indices
        top_ind = ceil(these_patches.top_row(j)*scaling);
        bot_ind = floor(these_patches.bot_row(j)*scaling);
        left_ind = ceil(these_patches.left_col(j)*scaling);
        right_ind = floor(these_patches.right_col(j)*scaling);

        % grab the patch
        this_patch = normalized_responses(top_ind:bot_ind,left_ind:right_ind, :);

        % and make into column
        responses_col = reshape(this_patch,[],n_filters);

        % do the assessment of anisotropy
        [pct_explained, filt_inds] = assess_anisotropy(responses_col, orient_inds);

        % phew, we did it. Now we need to get everything into tables.
        this_stat_table(j,1) = {pct_explained};
        this_stat_table(j,2) = {filt_inds(1)};
        this_stat_table(j,3) = {filt_inds(2)};
        this_stat_table(j,4) = {filt_inds(3)};
        this_stat_table(j,5) = {filt_inds(4)};
        this_stat_table(j,6) = {filt_inds(5)};
        this_stat_table(j,7) = {filt_inds(6)};
    end
    % update the whole stat table
    image_stat_table = [these_patches, this_stat_table];
    final_stat_table = [final_stat_table; image_stat_table];
    waitbar(i/numel(unique_names), f, 'running anisotropy analysis');
end
%%
im_data = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/data_spreadsheets/patch_anisotropy.xlsx');
sample_locs = readtable('/Users/ryan/Dropbox (Princeton)/reef_survey_project/nevada_qgis/smg_sample_gps_points.xls');

% first go through and average patches with the same name
% which stats can we average
for_averaging = [7,14:77];
% what are their headings
all_headings = im_data.Properties.VariableNames;
averaged_headings = all_headings(for_averaging);
info_headings = all_headings(1:2);

% and combine with other headers to make table
averaged_stat_table = cell2table(cell(0,numel(info_headings)+numel(averaged_headings)),'VariableNames',[info_headings, averaged_headings]);

unique_names = unique(im_data.image_name);
for i = 1:numel(unique_names)
    these_patches = im_data(strcmp(unique_names(i),im_data.image_name),:);
    % get which patches are unique
    unique_patches = unique(these_patches.patch_type);
    % loop through the unique patches
    for j = 1:numel(unique_patches)
        to_average = these_patches(strcmp(unique_patches(j),these_patches.patch_type),:);
        % and we need patch sizes for the weighting of the mean
        if size(to_average,1) > 1
            patch_sizes = (to_average.bot_row-to_average.top_row).*(to_average.right_col-to_average.left_col);
            weights = patch_sizes./sum(patch_sizes);
            mean_stats = sum(to_average{:,for_averaging}.*weights);
        else
            mean_stats = to_average{:,for_averaging};
        end

        % put everything in a table
        table_to_place = cell2table(cell(0,numel(info_headings)+numel(averaged_headings)),'VariableNames',[info_headings, averaged_headings]);
        table_to_place.image_name{1} = string(unique_names{i});
        table_to_place.patch_type = string(table_to_place.patch_type);
        table_to_place.patch_type(1) = string(unique_patches{j});
        table_to_place(1,3:end) = num2cell(mean_stats);

        averaged_stat_table = [averaged_stat_table;table_to_place];
    end
end

data_sample_names = cellfun(@(x) strsplit(x,'_'),averaged_stat_table.image_name, 'UniformOutput',false);

data_sample_nums = {};
for i = 1:numel(data_sample_names)
    data_sample_nums{i} = data_sample_names{i}{2};
    if numel(data_sample_nums{i}) < 2
        data_sample_nums{i} = ['0' data_sample_nums{i}];
    end
end

loc_sample_names = cellfun(@(x) strsplit(x,'G'),sample_locs.name, 'UniformOutput',false);

loc_sample_nums = {};
for i = 1:numel(loc_sample_names)
    loc_sample_nums{i} = loc_sample_names{i}{2};
end

stat_headings = {'latitude', 'longitude', 'msl', 'hae', 'litho_code'};
final_stat_table = cell2table(cell(size(averaged_stat_table,1),numel(stat_headings)),'VariableNames',[stat_headings]);


for i = 1:size(averaged_stat_table,1)
    this_loc_ind = strcmpi(data_sample_nums{i},loc_sample_nums);
    latit = {sample_locs.lat(this_loc_ind')};
    final_stat_table{i,1} = latit;
    final_stat_table(i,2) = {num2cell(sample_locs.long(this_loc_ind'))};
    final_stat_table(i,3) = {num2cell(sample_locs.msl(this_loc_ind'))};
    final_stat_table(i,4) = {num2cell(sample_locs.hae(this_loc_ind'))};
    final_stat_table(i,5) = {num2cell(sample_locs.litho_code(this_loc_ind'))};
end

big_shabang = [averaged_stat_table,final_stat_table];

%%
archaeos = strcmp(big_shabang.patch_type, 'archaeocyath');
nondolo_seds = strcmp(big_shabang.patch_type, 'nondolomitized sediments');
dolo_seds = strcmp(big_shabang.patch_type, 'dolomitized sediments');
ooids = strcmp(big_shabang.patch_type, 'ooids'); 
whole_rocks = strcmp(big_shabang.patch_type, 'whole rock');

shabang_nondolo = big_shabang(nondolo_seds,:);