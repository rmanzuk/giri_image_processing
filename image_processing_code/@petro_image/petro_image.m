classdef petro_image < handle
    % straight up properties the user can declare
    properties
        % the path to the master directory with everything
        main_path {mustBeText} = '/Volumes/ryan_ims/sms_images';
        % name of the sample
        sample_name {mustBeText} = '';
        % wavelenghts for the images, in order that the paths occur
        wavelengths(1, :) {mustBeNumeric} = [365, 470, 505, 530, 590, 625, 760, 940];
        % default extension for image files
        default_ext {mustBeText} = '.tif';
        % file with the geology/geography info csv sheet    
        geo_info_file {mustBeText} = 'geo_info.csv'
        % cell array with the sub paths to the image folders corresponding
        % to wavelengths
        im_subpaths cell {mustBeText} = {'/color_corrected_tifs/365nm', '/color_corrected_tifs/470nm', '/color_corrected_tifs/505nm', '/color_corrected_tifs/530nm', '/color_corrected_tifs/590nm', '/color_corrected_tifs/625nm', '/color_corrected_tifs/760nm', '/color_corrected_tifs/940nm'};
        % the strike of the image plane in the real world
        im_strike {mustBeNumeric}
        % dip of the image plane in the real world
        im_dip {mustBeNumeric}
        % the field lithology id of the sample
        field_litho
        % pointing of the arrow of strike in the image (NSEW style)
        strike_im_heading {mustBeNumeric}
        % the strike of local bedding planes in the field
        local_strike {mustBeNumeric}
        % dip of the local bedding planes in the field
        local_dip {mustBeNumeric}
        % latitude and longitude location of the image
        utm_xyz (1,3)
        % subpath to superpixel label matrix
        superpixel_subpath {mustBeText} = '/superpixel_inds';
        % how many superpixels we use
        n_superpixels {mustBeInteger}
        % stats for all of the superpixels 
        superpix_stats cell
        % which types of classes are present, identified by the user
        classes_present cell {mustBeText}
        % text-based labels of superpixels same length as number of indices
        % of superpixels
        class_labels cell
        % region coordinates with class labels that can be json formatted
        class_regions struct
        % the filter bank used for responses when getting stats
        filter_bank = makeLMfilters(6);
    end

    % properties that we can calculate from the others
    properties (Dependent)  
        num_channels
    end

    % all of the functions we can perform with this class
    methods
        % can calculate the number of channels based upon how many images
        % are in the object
        function value = get.num_channels(obj)
            value = numel(obj.im_subpaths);
        end

        % get any 3 channel image with an input vector of wavelengths
        function three_channel = get_3channel(obj,wave_vec)
            % option to just input the object and no vector. default output
            % is an rgb
            if nargin == 1
                im1 = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 625}, [obj.sample_name, obj.default_ext])));
                im2 = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 530}, [obj.sample_name, obj.default_ext])));
                im3 = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 470}, [obj.sample_name, obj.default_ext])));
                % concatenate together final product
                three_channel = cat(3,im1,im2,im3);
            else
                % quick check and make sure that the input wavelenthgs are
                % present in the object wavelengths
                if sum(ismember(wave_vec,obj.wavelengths)) == 3
                    im1 = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == wave_vec(1)}, [obj.sample_name, obj.default_ext])));
                    im2 = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == wave_vec(2)}, [obj.sample_name, obj.default_ext])));
                    im3 = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == wave_vec(3)}, [obj.sample_name, obj.default_ext])));
                    % concatenate together final product
                    three_channel = cat(3,im1,im2,im3);
                elseif numel(wave_vec) ~= 3
                     % if we don't have 3 inputs, gotta say something
                    error('improper number of wavelengths requested.')
                elseif ~all(ismember(wave_vec,obj.wavelengths))
                    % if we don't have all matches in the wavelents 
                    error('One or more wavelengths requested not present in this object.')
                end
            end
        end

        % get superpixel training data (function in separate file)
        superpixel_trainer(obj, n_supers);

        % assemble superpixel stats
        superpixel_stats(obj, n_supers);

        % fill in geology or geography info from csv file
        function obj = fill_geo(obj)
            % make up the composite path to the file
            geo_filename = fullfile(obj.main_path,obj.geo_info_file);
            % auto detect the import options for the csv
            opts = detectImportOptions(geo_filename);
            % and import
            raw_data = readtable(geo_filename,opts);
            % identify the corresponding sample
            this_sample = strcmpi(raw_data.sample,obj.sample_name);
            % now we can simply fill in
            obj.im_strike = raw_data.im_strike(this_sample);
            obj.im_dip = raw_data.im_dip(this_sample);
            obj.field_litho = raw_data.field_litho(this_sample);
            obj.local_strike = raw_data.local_strike(this_sample);
            obj.local_dip = raw_data.local_dip(this_sample);
            % lat and long are a bit more complicated, we'll want to switch
            % to utm for ease of use later
            % get the utm zone from the coordinates
            position = [raw_data.lat(this_sample),raw_data.long(this_sample)];
            utm_zone = utmzone(position);
            % get the geoid of the zon and construct its projection structure
            [ellipsoid,estr] = utmgeoid(utm_zone);
            utmstruct = defaultm('utm');
            utmstruct.zone = utm_zone;
            utmstruct.geoid = ellipsoid;
            utmstruct = defaultm(utmstruct);
            % and just do the conversion, put into the object right away
            [obj.utm_xyz(1), obj.utm_xyz(2)] = mfwdtran(utmstruct,raw_data.lat(this_sample),raw_data.long(this_sample));
            % last to put in is msl
            obj.utm_xyz(3) = raw_data.msl(this_sample);
        end

        % make and save superpixels 
        function [label_mat] =  get_superpix(obj, number_desired, rgb_im)
            % first make a list of all files (cell array) that should exist
            % with this number of superpixels
            f_names = cell(1,numel(number_desired));
            for i = 1:numel(number_desired)
                f_names{i} = fullfile(obj.main_path, obj.superpixel_subpath, num2str(number_desired(i)), [obj.sample_name, obj.default_ext]);
            end
            % if any of the desired superpixel files don't exist, we have
            % to go through and make them
            needed_superpix = ~isfile(f_names);
            if any(needed_superpix)
                % the third argument is optional input of an rgb image, so
                % make one if we don't have that argument
                if nargin < 3
                    % read in the images and tell the user
                    disp('reading in red, green, and blue images');
                    red_im = imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 625}, [obj.sample_name, obj.default_ext]));
                    green_im = imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 530}, [obj.sample_name, obj.default_ext]));
                    blue_im = imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 470}, [obj.sample_name, obj.default_ext]));
        
                    % concatenate into rgb and clear individual images
                    rgb_im = im2double(cat(3,red_im,green_im,blue_im));
                    clear red_im green_im blue_im
                end
                % now we need to loop through the missing files, and make
                % and store the superpixels
                needed_inds = find(needed_superpix);
                for i = 1:numel(needed_inds)
                    % we don't have superpixels, so we should make them
                    n_supers = number_desired(needed_inds(i));
                    disp(['Creating superpixelation for ', obj.sample_name ' with ' num2str(n_supers) ' super pixels.']);
                    label_mat = superpixels(rgb_im, n_supers);

                    % and update the list of n_superpixels stored in the
                    % object 
                    obj.n_superpixels = [obj.n_superpixels, n_supers];

                    % check if we even have a path
                    if isempty(obj.superpixel_subpath)
                        % if not make a simple, informative directory, and
                        % subdirectory for this number of superpixels
                        obj.superpixel_subpath = '/superpixel_inds';
                    end
                    
                    % make the subfolder for this particular n superpixels
                    % if needed
                    if ~isfolder(fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers)))
                        mkdir(fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers)));
                    end

                    % do the save
                    imwrite(uint16(label_mat),fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers), [obj.sample_name, obj.default_ext]));
                    % and inform the user
                    disp(['Superpixel indices saved in ', fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers), [obj.sample_name, obj.default_ext])]);
                end
            else 
                % otherwise we already have superpixels and should tell the
                % user
                disp(['file(s) already exist for the requested superpixelation']);
            end
        end
        
        % validate the number of wavelengths in the instance
        function validate_wavelengths(obj)
            if numel(obj.wavelengths) ~= numel(obj.im_subpaths) 
              warning('Not a corresponding wavelength for each channel');
            else
                disp('wavelenghts match up with number of channels, good to go!')
            end
        end
    end
end