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

supers_5000 = petro_superpix;

supers_5000.n_superpixels = 5000;
 
[~,supers_5000.sample_names] = cellfun(@fileparts, {input_dir(:).name}, 'UniformOutput', false);

supers_5000.superpixel_ind_path = '/Volumes/ryan_ims/sms_images/superpixel_inds/5000';


for i = 1:numel(input_dir)
    tic
    load(fullfile('/Volumes/ryan_ims/sms_images/petro_im_instances',input_dir(i).name))    
    [~,samp_name] = fileparts(input_dir(i).name);
    eval([samp_name '.fill_geo;'])
    eval(['this_pet = ' samp_name ';'])
    eval(['save(fullfile(''/Volumes/ryan_ims/sms_images/petro_im_instances'', samp_name), samp_name);']);
    eval(['clear  ' samp_name]);

    supers_5000.utm_xyz(i,:) = this_pet.utm_xyz;
    supers_5000.superpix_stats{i} = this_pet.superpix_stats{this_pet.n_superpixels == supers_5000.n_superpixels}; 

    supers_5000.im_strikes(i) = this_pet.im_strike;
    supers_5000.im_dips(i) = this_pet.im_dip;
    supers_5000.local_strikes(i) = this_pet.local_strike;
    supers_5000.local_dips(i) = this_pet.local_dip;
    supers_5000.field_lithos(i) = this_pet.field_litho;
    toc
end

supers_5000.make_stat_titles;