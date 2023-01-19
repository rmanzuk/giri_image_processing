function [gray_vals,wavelengths] = assess_colorbalance(colorcard_dir, im_ext)
% Assesses colorbalance of all images in a set based upon gray boxes in the
% input color card images. For this function, be sure to put the
% wavelenghts of each image in the image names.
%
% INPUT 
% colorcard_dir: String with the path to the directory 
% im_ext: String with the extension for the images
%
% OUTPUT
% gray_vals: number wavelengths x number gray boxes matrix with the average
% value within that box for each channel.
% wavelenghts: cell array with the character vectors with the wavelengts of
% each channel as gleaned from the file names.
% 
% Written by R.A. Manzuk
% 10/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end

    % get all of the images into a directory structure
    input_files = dir(fullfile(colorcard_dir, im_ext));
    input_files(strncmp({input_files.name}, '.', 1)) = []; %remove files in dir starting with '.'

    % the user should have named their images nicely corresponding to
    % wavelength, so let's get the names
    im_names = {input_files(:).name};
    [~,wavelengths,~] = cellfun(@fileparts,im_names,'UniformOutput',false);

    % and we'll work with an rgb to start, so let's get the indices for
    % those images
    red_ind = contains(im_names,'625');
    green_ind = contains(im_names,'530');
    blue_ind = contains(im_names,'470');

    % read in those images
    red_im = im2double(imread(fullfile(colorcard_dir,input_files(red_ind).name)));
    green_im = im2double(imread(fullfile(colorcard_dir,input_files(green_ind).name)));
    blue_im = im2double(imread(fullfile(colorcard_dir,input_files(blue_ind).name)));
    
    % make an rgb and show it
    rgb_im = cat(3,red_im,green_im,blue_im);
    figure()
    imshow(rgb_im)

    % and ask the user to input the centers of all gray boxes
    to_disp = 'please click the centers of all boxes that should be considered for gray values\nin ascending order of brighness. Press enter when done.\n';
    fprintf(to_disp);
    [col_coords,row_coords] = ginput();
    col_coords = round(col_coords);
    row_coords = round(row_coords);

    % we'll extract boxes of pixels from around the centers the box size
    % should be about 1/25 the size of the long dimension of the image.
    % Thus the dist from the center pixel should be 1/50.
    % Change this value if the color card or FoV ever changes
    dist_from_center = round(size(red_im,2)/50);

    clear rgb_im

    % now let's make a multispectral image of all inputs
    multispec_image = zeros(size(red_im,1), size(red_im,2), numel(input_files));
    for i = 1:numel(input_files)
        multispec_image(:,:,i) = im2double(imread(fullfile(colorcard_dir,input_files(i).name)));
    end

    % and then we should extract all the gray values from the means within
    % the squares from the centers we clicked earlier
    gray_vals = zeros(numel(input_files),numel(col_coords));
    for i = 1:numel(col_coords)
        % extract just the box in all channels
        this_box = multispec_image(row_coords(i)-dist_from_center:row_coords(i)+dist_from_center,col_coords(i)-dist_from_center:col_coords(i)+dist_from_center,:);
        % place the mean values in the gray_vals array
        gray_vals(:,i) = squeeze(mean(this_box,[1,2]));
    end
end