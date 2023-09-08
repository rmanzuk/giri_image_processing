% compilation of little bits of code used to construct instances of the
% object classes defined in this project. Could maybe make into formal
% functions within class definitions later
%
% R. A. Manzuk
%% compiling superpixel stats for all petro_images

input_dir = dir(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances', '*.mat'));

for i = 1:numel(input_dir)
    tic
    load(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances',input_dir(i).name));
    [~,samp_name] = fileparts(input_dir(i).name);
    eval([samp_name '.fill_geo;'])
    eval([samp_name '.superpixel_stats([1000, 5000, 10000, 25000, 50000]);'])
    eval(['save(fullfile(''/Volumes/ryan_ims/sms_images/petro_im_instances'', samp_name), samp_name);']);
    eval(['clear  ' samp_name]);
    toc
end

%% takeing petro_images and making petro_superpix

input_dir = dir(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances', '*.mat'));

supers_50000 = petro_superpix;

supers_50000.n_superpixels = 50000;
 
[~,supers_50000.sample_names] = cellfun(@fileparts, {input_dir(:).name}, 'UniformOutput', false);

supers_50000.superpixel_ind_path = '/Volumes/ryan_ims/sms_images/superpixel_inds/50000';


for i = 1:numel(input_dir)
    tic
    load(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances',input_dir(i).name))    
    [~,samp_name] = fileparts(input_dir(i).name);
%     eval([samp_name '.fill_geo;'])
    eval(['this_pet = ' samp_name ';'])
%     eval(['save(fullfile(''/Volumes/ryan_ims/sms_images/petro_im_instances'', samp_name), samp_name);']);
    eval(['clear  ' samp_name]);

    supers_50000.utm_xyz(i,:) = this_pet.utm_xyz;
    supers_50000.superpix_stats{i} = this_pet.superpix_stats{this_pet.n_superpixels == supers_50000.n_superpixels}; 
    supers_50000.class_labels{i} = this_pet.class_labels{this_pet.n_superpixels == supers_50000.n_superpixels};

    supers_50000.im_strikes(i) = this_pet.im_strike;
    supers_50000.im_dips(i) = this_pet.im_dip;
    supers_50000.local_strikes(i) = this_pet.local_strike;
    supers_50000.local_dips(i) = this_pet.local_dip;
    supers_50000.field_lithos(i) = this_pet.field_litho;
    supers_50000.im_bedding_angles(i) = this_pet.im_bedding_angle;
    toc
end

supers_50000.make_stat_titles;

%%
input_dir = dir(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances', '*.mat'));

for i = 1:numel(input_dir)
    tic
    load(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances',input_dir(i).name));
    [~,samp_name] = fileparts(input_dir(i).name);
    eval([samp_name '.fill_geo;'])
    eval([samp_name '.im_bedding_angle;'])


    eval(['save(fullfile(''/Volumes/ryan_ims/sms_images/petro_im_instances'', samp_name), samp_name);']);
    eval(['clear  ' samp_name]);
    toc
end

%% take all clicked geochem inputs and make them into one cell array for superpix
input_dir = dir(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances', '*.mat'));
input_dir(strncmp({input_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'

all_geochem_locs = cell(0,3);

for i = 1:numel(input_dir)
    tic
    load(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances',input_dir(i).name));
    [~,samp_name] = fileparts(input_dir(i).name);
    eval(['this_pet = ' samp_name ';'])

    these_geochem_locs = cell(size(this_pet.geochem_locs,1),3);

    these_geochem_locs(:,1) = {samp_name};
    these_geochem_locs(:,2:3) = this_pet.geochem_locs(:,:);

    all_geochem_locs = [all_geochem_locs; these_geochem_locs];

    eval(['clear  ' samp_name]);
    toc
end

% and load and save in the superpix instances
load('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_1000.mat')
supers_1000.all_geochem_locs = all_geochem_locs;
save('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_1000.mat', 'supers_1000', '-v7.3')
clear supers_1000

load('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_5000.mat')
supers_5000.all_geochem_locs = all_geochem_locs;
save('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_5000.mat', 'supers_5000', '-v7.3')
clear supers_5000

load('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_10000.mat')
supers_10000.all_geochem_locs = all_geochem_locs;
save('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_10000.mat', 'supers_10000', '-v7.3')
clear supers_10000

load('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_25000.mat')
supers_25000.all_geochem_locs = all_geochem_locs;
save('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_25000.mat', 'supers_25000', '-v7.3')
clear supers_25000

load('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_50000.mat')
supers_50000.all_geochem_locs = all_geochem_locs;
save('/Volumes/ryan_ims/sms_images/petro_superpix_instances/supers_50000.mat', 'supers_50000', '-v7.3')
clear supers_50000