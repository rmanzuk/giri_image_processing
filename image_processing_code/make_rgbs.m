function make_rgbs(red_folder, green_folder, blue_folder, destination_folder, im_in_ext, im_out_ext, scaling, varargin)
%Applies gray or lcc/white card image adjustments to a set of sample images
%taken at the same wavelength
%
% INPUT 
% red_folder: String with the full path to the red channel images
% green_folder: String with the full path to the green channel images
% blue_folder: String with the full path to the blue channel images
% destination_folder: String with the full path to the folder where you
% would like the final rgbs to end up
% im_in_ext: String with the extension for the input image files
% im_out_ext: String with the extension for the desired output rgb images
% scaling: fraction by which you would like the images resized. Use 1 if
% you would like rgbs to be same resolution as inputs.
% name_map: (optional) n_sample x 2 cell array with the file names of the
% input images and the associated names you would like for the rgb images
% in the 2nd column. Function will otherwise base naming off of the input
% red image.
% mask_path: (optional) character vector with the full path to the folder
% which contains masks to be applied to the images. Likely to get rid of
% background and specular patches
%
% OUTPUT
% there are no outputs, but the rgb images will be in the indicated
% output folder
%
% Written by R.A. Manzuk
% 10/12/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_in_ext(1) ~= '*'
        im_in_ext = ['*' im_in_ext];
    end

    % start with false logicals for optional inputs
    name_mapped = false;
    masking = false;
    
    % gotta parse optional arguments
    if numel(varargin) > 0 
        % first see what the types of the variables are
        var_classes = cell(size(varargin));
        for i = 1:numel(varargin)
            var_classes{i} = class(varargin{i});
        end
        % a name_map would be if any of these classes were a cell
        if any(strcmp(var_classes, 'cell'))
            map_ind = strcmp(var_classes, 'cell');
            name_map = varargin{map_ind};
            name_mapped = true;
        end
        % a mask path woud be a character vector
        if any(strcmp(var_classes, 'char'))
            mask_ind = strcmp(var_classes, 'char');
            mask_path = varargin{mask_ind};
            masking = true;
        end
    end

    % make directories
    red_dir = dir(fullfile(red_folder, im_in_ext));
    red_dir(strncmp({red_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
    green_dir = dir(fullfile(green_folder, im_in_ext));
    green_dir(strncmp({green_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
    blue_dir = dir(fullfile(blue_folder, im_in_ext));
    blue_dir(strncmp({blue_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'

    % optional directory for masks if exists
    if masking 
        mask_dir = dir(fullfile(mask_path, im_in_ext));
        mask_dir(strncmp({mask_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
    end
    
    % and then go through all of the images, read them, concatenate them
    % and save them

    % Make a waitbar to keep track of how long it's taking
    f = waitbar(0, 'Making RGB images');
    for i = 1:numel(red_dir)
       
        this_red = im2double(imread(fullfile(red_folder, red_dir(i).name)));
        this_green = im2double(imread(fullfile(green_folder, green_dir(i).name)));
        this_blue = im2double(imread(fullfile(blue_folder, blue_dir(i).name)));
        
        this_rgb = cat(3, this_red, this_green, this_blue);

        % if we're masking, gotta apply the mask
        if masking
            this_mask = double(imread(fullfile(mask_path,mask_dir(i).name)));
            this_rgb = this_rgb.*this_mask;
        end

        if scaling < 1
            this_rgb = imresize(this_rgb,scaling);
        end


        if name_mapped
            % in this case the user provided a naming key. we'll search off
            % the red image
            [~,red_name,~] = fileparts(red_dir(i).name); 
            [~,search_names,~] = cellfun(@fileparts,name_map(:,1),'UniformOutput',false);
            map_ind = ismember(search_names,red_name);
            
            % now we can make the output name
            output_name = [name_map{map_ind,2} im_out_ext];
        else
            % if the user didn't provide a key, we'll just use the name of
            % the red image
            [~,red_name,~] = fileparts(red_dir(i).name); 
            output_name = [red_name im_out_ext];
        end

        imwrite(this_rgb, [destination_folder '/' output_name]);
        waitbar(i/numel(red_dir),f, 'Making RGB images');
    end

    % close the waitbar
    close(f);

end