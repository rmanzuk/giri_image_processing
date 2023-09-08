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
        % the average strike of bedding for the whold outcrop
        regional_strike = 38.1;
         % the average dip of bedding for the whold outcrop
        regional_dip = 34.7;
        % subpath to superpixel label matrix
        superpixel_subpath {mustBeText} = '/superpixel_inds';
        % subpath to tif background mask of this image
        background_mask_path {mustBeText} = '/masks';
        % boundary coordinates of the mask to be stored easily. NOTE: This
        % value does not denote if any of the boundary sections are holes
        % in an element of the mask. That must be handled elsewhere.
        background_mask_bounds
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
        % locations of geochemistry drill points with labels
        geochem_locs cell
        % the filter bank used for responses when getting stats
        filter_bank = makeLMfilters(6);
        % file name of sheet with carbon and oxygen data
        carb_ox_file {mustBeText} = '/geochem/carbon_oxygen_master.csv';
        % file name of sheet with element concentration data
        element_file {mustBeText} = '/geochem/element_concentrations_master.csv';
        % cell array with the indices of superpixels associated with each
        % geochemical click
        geochem_superpix_inds
    end

    % properties that we can calculate from the others
    properties (Dependent)  
        num_channels
        % angle between the image plane and the local bedding plane
        im_bedding_angle {mustBeNumeric}
        % cell array with geochemical
        geochem_data {mustBeCell}
        % the file name for point count data will depend on the sample name
        point_count_file {mustBeText}
        % cell array with point count data
        point_count_data {mustBeCell}
    end

    % all of the functions we can perform with this class
    methods
        %%
        % can calculate the number of channels based upon how many images
        % are in the object
        function value = get.num_channels(obj)
            value = numel(obj.im_subpaths);
        end
        
        %% 
        % make the file name for point counting based upon sample name
        function value = get.point_count_file(obj)
            value = ['/point_counts/' obj.sample_name '_Point_Counting_regular_300_coordinants.csv'];
        end
        %% 
        % from the image orientation stats and local or regional bedding,
        % can calcualate the angle between the image and the bedding plane
        function value = get.im_bedding_angle(obj)
            % first check that we have enough information to proceed
            has_im_orientation = all([~isnan(obj.im_strike), ~isnan(obj.im_dip), ~isempty(obj.im_strike), ~isempty(obj.im_dip)]);
            has_local_bedding = all([~isnan(obj.local_strike), ~isnan(obj.local_dip), ~isempty(obj.local_strike), ~isempty(obj.local_dip)]);
            has_regional_bedding = all([~isnan(obj.regional_strike), ~isnan(obj.regional_dip), ~isempty(obj.regional_strike), ~isempty(obj.regional_dip)]);
            % if we don't have an image orientation, we know not to
            % proceed, and just set the value to nan
            if ~has_im_orientation
                value = NaN;
            elseif ~has_local_bedding && ~has_regional_bedding
                % also if we have neither local or regional bedding, we
                % can't do anything
                value = NaN;
            else
                % here we're good to go, so just need to decide if we're
                % using local or regional bedding (local is the priority)   
                if has_local_bedding
                    % and we need the strike to be a dip direction to
                    % calculate normal
                    bedding_ddir = obj.local_strike - 90;
                    bedding_dip = obj.local_dip;
                else 
                    % and we need the strike to be a dip direction to
                    % calculate normal
                    bedding_ddir = obj.regional_strike - 90;
                    bedding_dip = obj.regional_dip;
                end
                
                % with bedding figured out, we also need the image strike
                % as a dip direction
                image_ddir = obj.im_strike - 90;

                % extract plane normals from strikes and dips
                im_normal_x = -sind(obj.im_dip)*sind(image_ddir);
                im_normal_y = sind(obj.im_dip)*cosd(image_ddir);
                im_normal_z = -cosd(obj.im_dip);
                im_normal = [im_normal_x,im_normal_y,im_normal_z];
                bedding_normal_x = -sind(bedding_dip)*sind(bedding_ddir);
                bedding_normal_y = sind(bedding_dip)*cosd(bedding_ddir);
                bedding_normal_z = -cosd(bedding_dip);
                bedding_normal = [bedding_normal_x,bedding_normal_y,bedding_normal_z];

                % calculate the angle at which bedding dips into the image,
                % via the normals
                bedding_wrt_im = atan2d(norm(cross(im_normal,bedding_normal)),dot(im_normal,bedding_normal));

                % calculate the sign of the plane angles 
                % the normal of the plane fit to the two surface normal vectors is just
                % their cross product
                plane_angle_cross = cross(im_normal,bedding_normal);
                % and fill the final signed value, based on the dot product
                value = bedding_wrt_im * sign(dot(plane_angle_cross, plane_angle_cross));
            end
        end

        %%
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
        %%
        % get superpixel training data (function in separate file)
        superpixel_trainer(obj, n_supers);

        % assemble superpixel stats
        superpixel_stats(obj, n_supers);
        %%
        % create a polygon from the mask tif to store it in a small way
        function obj = fill_background_bounds(obj)
            % start by loading in the mask
            mask_logical = imread(fullfile(obj.main_path,obj.background_mask_path,[obj.sample_name obj.default_ext]));
            
            % and get all of the region boundaries 
            region_bounds = bwboundaries(mask_logical);

            % the masks may have some spurious small scale features, so we
            % can get rid of those based upon the number of elements of the
            % boundary. The number of elements threshold is arbitrary
            bound_numel = cell2mat(cellfun(@numel,region_bounds,'UniformOutput',false));
            not_trash = bound_numel > 300;
            
            % and properly fill in
            obj.background_mask_bounds = region_bounds(not_trash);
        end
        %%
        % given either a background mask path or the boundaries, label all
        % background superpixels for the desired superpixel level
        function obj = label_background(obj, number_desired, background_source)
            
            % assemble the list of expected file names for the superpixel
            % label mats 
            f_names = cell(1,numel(number_desired));
            for i = 1:numel(number_desired)
                f_names{i} = fullfile(obj.main_path, obj.superpixel_subpath, num2str(number_desired(i)), [obj.sample_name, obj.default_ext]);
            end

            % if any of the desired superpixel files don't exist, we have
            % to go through and make them
            needed_superpix = ~isfile(f_names);
            if sum(needed_superpix) > 0
                disp('creating missing superpixel label images');
                obj.get_superpix(number_desired(needed_superpix));
            end

            % now we need the background mask, which we can load depending
            % on the source specified by the user
            if strcmpi(background_source, 'mask')
                % here the user wants to use the mask image saved at the
                % path specified in the object
                % so load the image
                disp('loading mask image')
                mask_logical = imread(fullfile(obj.main_path,obj.background_mask_path,[obj.sample_name obj.default_ext]));
            elseif strcmpi(background_source, 'boundaries')
                % here the user wants to use the polygon stored in the
                % images

                % check we even have coordinates
                if isempty(obj.background_mask_bounds)
                    % if not, make them
                    obj.fill_background_bounds
                end

                % time to make the mask
                disp('making mask from boundary coordinates')

                % gotta worry about holes, but only if there is more than one
                % boundary
                if numel(obj.background_mask_bounds) > 1
                    % need to decide if any of these are holes by looking for overlap
                    % first construct a polyshape vector
                    polyvec = repmat(polyshape,1,numel(obj.background_mask_bounds));
                    warning('off','all')
                    for j = 1:numel(obj.background_mask_bounds)
                        polyvec(j) = polyshape([obj.background_mask_bounds{j}(:,1), obj.background_mask_bounds{j}(:,2)]);
                    end
                    warning('on','all');
        
                    % and test for overlap
                    is_overlapping = overlaps(polyvec);
                    
                    % the diagonal of the overlap matrix is automatically true, but
                    % we want it to be false. so get the linear indices of the
                    % diagonal and set them to false
                    diag_inds = linspace(1,numel(is_overlapping),size(is_overlapping,1));
                    is_overlapping(diag_inds) = false;
            
                    % if we do have overlap, need to do something about it
                    if sum(is_overlapping,'all') > 0
                        holes = false(1,numel(obj.background_mask_bounds));
                        % because the overlapping matrix is symmetric, we only need to
                        % iterate over one dimension
                        hole_pairs = nchoosek(find(any(is_overlapping)),2);
                        for j = 1:size(hole_pairs,1)
                            % we'll use areas to decide which one is the hole
                            area1 = area(polyvec(hole_pairs(j,1)));
                            area2 = area(polyvec(hole_pairs(j,2)));
                            [~,hole_ind] =  min([area1, area2]);
                            % and set the proper one as the hole in the logical
                            holes(hole_pairs(j,hole_ind)) = true;
                        end
                    else
                        % in this case no holes
                        holes = false(1,numel(obj.background_mask_bounds));
                    end
                else
                    % in this case there is only one element, can't be a
                    % hole, so we just need a 1x1 false flag.
                    holes = false(1,1);
                end

                % now we have holes sorted out, so we can make the mask
                % from the boundaries
                % start by setting up separate array for good boundaries
                norm_bounds = obj.background_mask_bounds(~holes);

                % and we'll make a mask stack to compress down into a
                % composite
                % to get the size, we need the size of an image. could do
                % in a faster way in the futre
                test_im = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{1}, [obj.sample_name, obj.default_ext])));
                mask_stack = zeros(size(test_im,1), size(test_im,2), numel(norm_bounds));

                % iterate through and fill up the mask stack
                for i = 1:numel(norm_bounds)
                    % when making masks, sometimes there were weird strips
                    % near the edges of a few pixels in thickness that get
                    % masked out, so can correct for those here
                    these_x = norm_bounds{i}(:,2);
                    these_x(these_x <= 4) = 1;
                    these_x(these_x >= size(test_im,2)-3) = size(test_im,2);
                    these_y = norm_bounds{i}(:,1);
                    these_y(these_y <= 4) = 1;
                    these_y(these_y >= size(test_im,1)-3) = size(test_im,1);

                    mask_stack(:,:,i) = poly2mask(these_x, these_y, size(test_im,1), size(test_im,2));
                end
                
                % and the final mask is just the where any of the mask
                % stack is true in the 3rd dimension
                mask_logical = any(mask_stack,3);

                % and if we had any holes earlier, do the stack thing to
                % get a hole logical and get rid of it from the mask
                if sum(holes) > 0 
                    % make a separate hole bounds cell array
                    hole_bounds = obj.background_mask_bounds(holes);
                    hole_stack = zeros(size(test_im,1), size(test_im,2), numel(hole_bounds));

                    % iterate through and fill up the mask stack
                    for i = 1:numel(hole_bounds)
                        hole_stack(:,:,i) = poly2mask(hole_bounds{i}(:,2), hole_bounds{i}(:,1), size(test_im,1), size(test_im,2));
                    end
                    
                    % collapse the logical
                    hole_logical = any(hole_stack,3);

                    % and just turn the holes to false in the mask
                    mask_logical(hole_logical) = false;
                end
            else
                error('Invalid background source string for matching.')
            end
            
            disp('labeling background superpixels.')
            % okay, we have a mask logical, now just have to mark the right
            % superpixels, which we will do with column versions of the
            % mask and superpixel indices
            mask_col = mask_logical(:);
            % do this iteratively through all the numbers of superpixels
            % requested
            for i = 1:numel(number_desired)
                % load the label matrix and make it a column
                label_mat = imread(fullfile(obj.main_path, obj.superpixel_subpath, num2str(number_desired(i)), [obj.sample_name, obj.default_ext]));
                label_col = label_mat(:);

                % the superpixels to mask are the unique entries in the
                % superpixel column where the mask is false
                as_background = unique(label_col(~mask_col));

                % and just do the labeling. start by making sure there is
                % an array at the right spot to fill
                this_super_ind = find(obj.n_superpixels == number_desired(i));
                if numel(obj.class_labels) <= this_super_ind || isempty(obj.class_labels{this_super_ind})
                    obj.class_labels{this_super_ind} = cell(max(label_col),1);
                end
                obj.class_labels{this_super_ind}(as_background) = {'background'};
            end
        end
        %% 
        % ginput geochemical drill spots
        function obj = input_geochem_locs(obj)
            % maybe there already some geochem locs here. If so, we need to
            % ask the user if they want to clear and start fresh, or
            % update.
            updating = false;
            if ~isempty(obj.geochem_locs)
                prompt = 'geochemistry locations already exist for this sample. \nwould you like to clear and start new (c), or update the existing list(u)?\n';
                user_desire = input(prompt,'s');

                % just compare the users input with possibilites
                if strcmpi(user_desire, 'c')
                    obj.geochem_locs = {};
                elseif strcmpi(user_desire,'u')
                    updating = true;
                else
                    % here the input was wrong, so tell the user to try
                    % again
                    disp('not a valid response. Try the function again.')
                    return
                end
            end

            % we need an rgb image to display and ask for ginputs
            rgb_im = obj.get_3channel;
            % and we want to ask the user which phase labels they have
            prompt = 'please input the phase labels present in this sample. \nlist duplicates with no additional qualifiers. \n';
            raw_phases = input(prompt,"s");
            % now we're ready to ginput. We'll store everything in a 2
            % column cell array where the first column is the labels and
            % the second column is a vector of the coordinates
            loc_cell = cell(numel(raw_phases),2);
            % and we want to track any duplicates, so we can add a
            % qualifyer
            [unique_phases, first_inds] = unique(raw_phases);
            duplicate_inds = setdiff(1:numel(raw_phases), first_inds);

            % and we may have more than 2 of a particular phase, so we need
            % a counter for each phase to update in that case
            counters = ones(1,numel(unique_phases));

            % and if we're updating a previous set, the new set may include
            % duplicates of existing phases, so account for that
            if updating
                % here we'll just go through each of the new phases and see
                % if we need to do anything
                for i = 1:numel(unique_phases)
                    where_exists = strfind(obj.geochem_locs(:,1),raw_phases(i));
                    % if there are no duplicates, we're good, otherwise we
                    % need to do something
                    if sum(cellfun(@isempty,where_exists)) ~= numel(where_exists)
                        % first, we need to indicate that this index
                        % actually is a duplicate
                        duplicate_inds = [duplicate_inds, first_inds(i)];
                        
                        % and then we need to update the counter to be the
                        % max of the existing qualifiers
                        % so extract all the necessary qualifiers
                        is_this_phase = ~cellfun(@isempty,where_exists);
                        % cells with 2 characters are going to have a
                        % number we want to look at
                        has_2 = cellfun(@(x) numel(x)==2, obj.geochem_locs(is_this_phase,1));
                        % and extract the previous qualifiers
                        is_phase_cell = obj.geochem_locs(is_this_phase,1);
                        prev_quals = cellfun(@(x) x(2), is_phase_cell(has_2,1));
                        
                        % and update the counters, if prev_quals is empty,
                        % that means there were no duplicates previously,
                        % so this is the first duplicate
                        if isempty(prev_quals)
                            counters(i) = 1;
                        else
                            counters(i) = max(str2num(prev_quals));
                        end
                    end
                end
            end

          
            % populate label column in a for loop so we can handle
            % particularities
            for i = 1:numel(raw_phases)
                % and add the qualifier if needed
                if ismember(i, duplicate_inds)
                    % find the counter index of this character
                    counter_ind = strfind(unique_phases,raw_phases(i));
                    % and the qualifier is just 1 plus the existing counter
                    % value
                    qualifier = num2str(counters(counter_ind) + 1);
                    % and then just update the counter and fill the cell
                    counters(counter_ind) = counters(counter_ind) + 1;
                    loc_cell{i,1} = [raw_phases(i) qualifier];
                else
                    loc_cell{i,1} = raw_phases(i);
                end
                
            end

            % now we ginput, can do this on a loop with nice display
            % prompts
            f = figure();
            imshow(rgb_im)
            for i = 1:size(loc_cell,1)
                % tell the user which phase to click
                disp(['Please click the location of spot ' loc_cell{i,1} '.']);
                % and get the ginput
                [loc_col, loc_row] = ginput(1);
                % and then just put in the cell array
                loc_cell{i,2} = [round(loc_row), round(loc_col)];
            end
            close(f)
            
            % and we're good, just put into the object at the right spot
            % depending on if we're updating or not
            if updating
                obj.geochem_locs = cat(1,obj.geochem_locs, loc_cell);
            else
                obj.geochem_locs = loc_cell;
            end
        end
        
        %%
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
        %%
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
        %%
        % load in carbon and oxygen isotope data from master spreadsheet
        function value = get.geochem_data(obj)
            % load in the master spreadsheets
            raw_carbox = readtable(fullfile(obj.main_path,obj.carb_ox_file));
            % issue with some of the element concentrations reading in as
            % string, so fix import options
            opts = detectImportOptions(fullfile(obj.main_path,obj.element_file));
            new_formats = cell(1, numel(opts.VariableTypes)-1);
            new_formats(:) = {'double'};
            opts.VariableTypes(2:end) = new_formats;
            raw_elements = readtable(fullfile(obj.main_path,obj.element_file), opts);
            % split the strings in the sample name column so we can find
            % this sample's data
            split_names = cellfun(@(x) strsplit(x,'_'), raw_carbox.export_id, 'UniformOutput',false);
            % and get the inds for non split names (likely samples not
            % standards)
            is_sample = cellfun(@numel, split_names) > 1;
            % okay, finally ready to see what the part of the sample name
            % is before the underscore
            just_samples = vertcat(split_names{is_sample});
            % and get a logical for those that match this sample
            sample_num = strsplit(obj.sample_name, '_');
            is_this_sample = strcmpi(just_samples(:,1), num2str(sample_num{2}));
            % now we're ready to put together the cell arry with geochem
            % data
            % first row of cells will be the headers of the columns
            data_cell = cell(2,width(raw_carbox)+25);
            % I'm cheating here, preprescribing the number of data columns
            % I'll be taking from the element concentrations data
            field_names = fieldnames(raw_carbox);
            field_names{1} = 'phase_code';
            
            data_cell(1,1:7) = field_names(1:width(raw_carbox))';
            % and fill in everything else. Could make this neater and more
            % flexible in the future.
            data_cell(2,1) = {just_samples(is_this_sample,2)};
            sample_data = raw_carbox(is_sample,:);
            data_cell(2,2) = {sample_data.line(is_this_sample)};
            data_cell(2,3) = {sample_data.delta_13C_mean(is_this_sample)};
            data_cell(2,4) = {sample_data.delta_13C_std(is_this_sample)};
            data_cell(2,5) = {sample_data.delta_18O_mean(is_this_sample)};
            data_cell(2,6) = {sample_data.delta_18O_std(is_this_sample)};
            data_cell(2,7) = {datetime(sample_data.date_run(is_this_sample),"Format","MM/dd/uu")};

            % now that we're through carbon and oxygen isotopes, need to
            % fill in element concentrations

            % split the strings in the sample name column so we can find
            % this sample's data
            split_names = cellfun(@(x) strsplit(x,'_'), raw_elements.name, 'UniformOutput',false);
            % sometimes have a third element of the split names due to a
            % note in the spreadsheet, only want the first two
            for i = 1:numel(split_names)
                if numel(split_names{i}) ~= 2
                    split_names{i} = split_names{i}(1:2);
                end
            end
            % get rid of the cell array wrapping the split names
            split_names = vertcat(split_names{:});
            % and get a logical for those that match this sample
            sample_num = strsplit(obj.sample_name, '_');
            is_this_sample = strcmpi(split_names(:,1), num2str(sample_num{2}));
            % now we're ready to put together the cell arry with geochem
            % data
            % add the headers to the data cell
            field_names = fieldnames(raw_elements);
            data_cell(1,8:end) = field_names(12:36)';
            % and fill in everything else. Could make this neater and more
            % flexible in the future.
            sample_data = raw_elements(is_this_sample, 2:end);
            % final step, because we want to match the phase order of the
            % existing carbon data, we need to get a sorting vector based
            % upon that
            existing_phases = data_cell{2,1};
            elements_phases = split_names(is_this_sample,2);
            [~,sorting_idx] = ismember(elements_phases,existing_phases);
            for i = 1:25
                data_cell(2,i+7) = {table2array(sample_data(sorting_idx,i+10))};
            end

            % and set the value
            value = data_cell;
        end


        %%
        % load in point count data from spreadsheet
        function value = get.point_count_data(obj)
            % load in the master spreadsheet. It'll be formatted weird, so
            % supress the warnings
            warning('off','all')
            raw_point_counts = readtable(fullfile(obj.main_path,obj.point_count_file));
            warning('on','all')
            
            % now we're ready to put together the cell arry with point
            % count data
            % first row of cells will be the headers of the columns, but
            % we'll make the xy positions combined
            data_cell = cell(2,width(raw_point_counts)-1);
            field_names = fieldnames(raw_point_counts);
            field_names{2} = 'xy_pos';
            data_cell(1,:) = field_names(1:width(raw_point_counts)-1)';
            % and fill in everything else. Could make this neater and more
            % flexible in the future.
            data_cell(2,1) = {raw_point_counts.Class};
            data_cell(2,2) = {[raw_point_counts.PosX, raw_point_counts.PosY]};

            % and set the value
            value = data_cell;
        end

        %% 
        % display the point counts overlain on the image
        function overlay_point_counts(obj)
            % first need an rgb image for the sample
            rgb_im = obj.get_3channel;

            % and get the mask
            im_mask = imread(fullfile(obj.main_path, obj.background_mask_path, [obj.sample_name obj.default_ext]));
            
            % and then the unique entries in the class portion of the point
            % count data
            class_col = strcmpi(obj.point_count_data(1,:), 'Class');
            unique_classes = unique(obj.point_count_data{2,class_col});

            % we may need more than 1 symbol to display points uniquely if
            % we have more classes than colors, so get the color order for
            % keeping track
            color_num = size(get(gca,'colororder'), 1);
            symbols = {'o', 'square', 'diamond', '^', 'v'};
            
            % cool, just display the image and iterate through the columns
            % to display everything
            figure()
            imshow(rgb_im .* im_mask)
            % note point counts were made on 5000 pixel wide images scale
            % images, so need a scaling factor
            mult_fac = size(rgb_im,2)/5000;
            hold on
            for i = 1:numel(unique_classes)
                these_points = strcmpi(obj.point_count_data{2,class_col}, unique_classes{i});
                % scatter keeping track of symbols
                symbol_num = ceil(i/color_num);
                scatter(obj.point_count_data{2,~class_col}(these_points,1)*mult_fac, obj.point_count_data{2,~class_col}(these_points,2)*mult_fac, 20, 'filled', symbols{symbol_num}, 'DisplayName', unique_classes{i}, 'MarkerEdgeColor',[0 0 0]);
            end
            legend
            hold off
        end
        %% 
        % overlay geochemical data on the sample
        function overlay_geochem(obj,chem_variable)
            % first need an rgb image for the sample
            rgb_im = obj.get_3channel;

            % and get the mask
            im_mask = imread(fullfile(obj.main_path, obj.background_mask_path, [obj.sample_name obj.default_ext]));
            
            % because carbon data and their localities are input
            % separately, we need to match the clicked localities to the
            % isotope data
            % so run the string comparison between the two lists,
            % iteratively is easy enough. Maybe could be cooler some other
            % way
            sample_col = strcmpi(obj.geochem_data(1,:), 'phase_code');
            point_match = zeros(numel(obj.geochem_locs(:,1)),1);
            for i = 1:numel(obj.geochem_locs(:,1))
                match_logical = strcmpi(obj.geochem_locs{i,1}, obj.geochem_data{2,sample_col});
                % maybe the '1' is missing on the point position label, so
                % check
                if sum(match_logical) == 0
                    match_logical2 = strcmpi([obj.geochem_locs{i,1}, '1'], obj.geochem_data{2,sample_col});
                    point_match(i) = find(match_logical2);
                else
                    point_match(i) = find(match_logical);
                end
            end
            
            % so at this point, the indices in point_match are the indices
            % of the geochem data where they would match the point locations

            % now we just need to sort out which geochemical variable we're
            % plotting
            if strcmpi(chem_variable,'carb')
                data_col = strcmpi(obj.geochem_data(1,:), 'delta_13C_mean');
                color_label = '\delta^{13}C';
            elseif strcmpi(chem_variable,'ox')
                data_col = strcmpi(obj.geochem_data(1,:), 'delta_18O_mean');
                color_label = '\delta^{18}O';
            elseif sum(strcmpi(chem_variable,obj.geochem_data(1,:))) == 1
                data_col = strcmpi(chem_variable,obj.geochem_data(1,:));
                color_label = chem_variable;
            else
                error('not a valid geochemical variable.')
            end

            % now we're ready to plot
            figure()
            imshow(rgb_im .* im_mask)
            hold on
            % and we just need to extract the geochem locs into a usable
            % array
            point_locs = vertcat(obj.geochem_locs{:,2});
            % may have to deal with nans
            is_data = ~isnan(obj.geochem_data{2,data_col});
            scatter(point_locs(is_data,2), point_locs(is_data,1), 50, obj.geochem_data{2,data_col}(is_data), 'filled','MarkerEdgeColor',[0 0 0]);
            % set a diverging colormap for this plot
            burd_raw = [33, 102, 172; 67, 147, 195; 146, 197, 222; 209, 229, 240; ...
                247, 247, 247; 253, 219, 199; 244, 165, 130; 214, 96, 77; 178, 24, 43]./255;
            cmap = imresize(burd_raw,[256, 3], 'bilinear');
            colormap(cmap)
            h = colorbar;
            ylabel(h, color_label)
            hold off
        end
      
        %%
        % validate the number of wavelengths in the instance
        function validate_wavelengths(obj)
            if numel(obj.wavelengths) ~= numel(obj.im_subpaths) 
              warning('Not a corresponding wavelength for each channel');
            else
                disp('wavelenghts match up with number of channels, good to go!')
            end
        end
        %%
        % from the ginput drill spots, automatically assign the superpixel
        % index associated with that spot for each number of superpixels
        function set_geochem_super_inds(obj)
            % only do any thing if there are locs to begin with 
            if numel(obj.geochem_locs) > 0 
                % start by just setting up the cell array to receive the
                % indices. basically jsut copy the loc one. 
                clicked_inds = cell(size(obj.geochem_locs));
                clicked_inds(:,1) = obj.geochem_locs(:,1);
                
                % we just need to loop through all of the n_superpix in the
                % object. No reason to do this faster because we have to load
                % through the loading (or creation) of superpix, which would be
                % the rate-limiting step
                loc_array = vertcat(obj.geochem_locs{:,2});
                inds_array = zeros(size(clicked_inds,1), numel(obj.n_superpixels));
                for i = 1:numel(obj.n_superpixels)
                    n_supers = obj.n_superpixels(i);
                    % now we have to check if there are superpixels already made
                    superpix_fname = fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers), [obj.sample_name, obj.default_ext]);
                    if isfile(superpix_fname)
                        disp('loading in superpixel indices');
                        % if they are there we can just load them
                        label_mat = imread(superpix_fname);
                    else
                        label_mat = obj.get_superpix(n_supers, rgb_im);
                    end
    
                    % easy enough with the label mat loaded, so just extract
                    % the superpixel indices at the clicked locations
                    linear_locs = sub2ind(size(label_mat),loc_array(:,1),loc_array(:,2));
                    inds_array(:,i) = label_mat(linear_locs);
                end
    
                % all done, just put in the cell array
                clicked_inds(:,2) = num2cell(inds_array,2);
                obj.geochem_superpix_inds = clicked_inds;
            end
        end
    end
end